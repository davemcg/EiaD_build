### Function to create count data for non-GTEX samples -----

create_count_data_frames <- 
  function(recount3_library_url, 
           projects_vector, 
           count_file_name, 
           aggregated_count_file_name, 
           mapping_file_name, 
           metadata,
           empty_cache = TRUE) {
    
    if (empty_cache){
      print("Clearing recount3 cache")
      recount3_cache_rm()
    }
    #Indicating which recount library the function will pull from
    options(recount3_url = recount3_library_url)
    human_projects <- available_projects()
    #Create an empty list to gather the count data
    counts_data_list <- list()
    #Create an empty list to gather mapping information
    mapping_data_list <- list()
    # creates list of count data frames
    for (i in projects_vector) {
      success <- NULL
      attempt <- 1
      while( is.null(success) && attempt <= 3 ) {
        try({
          print(paste(i, 'attempt number:', attempt))
          attempt <- attempt + 1
          rse_gene = create_rse_manual(i, annotation = "gencode_v29")
          assays(rse_gene)$counts <- transform_counts(rse_gene)
          
          rse_gene_counts <- recount::getTPM(rse_gene, 
                                             length_var = 'bp_length',
                                             mapped_var = 'recount_qc.star.all_mapped_reads') %>% 
            data.frame()
          rse_gene_counts <- rownames_to_column(rse_gene_counts, var = "gene_id")
          
          row.names(rse_gene_counts) <- rse_gene_counts$gene_id
          rse_gene_counts <- rse_gene_counts[,-1, drop=FALSE]
          counts_data_list[[i]] <- rse_gene_counts
          #Input mapping data
          mapping_data <- colData(rse_gene) %>% as.data.frame()
          mapping_data <- mapping_data %>% mutate(across(everything(), as.character))
          mapping_data_list[[i]] <- mapping_data
          success <- 'YES'
          
        })
      }
    }
    
    # confirm gene name equality across each data frame
    gene_name_list <- list()
    for (i in 1:length(counts_data_list)){
      gene_name_list[[i]] <- row.names(counts_data_list[[i]])
    }
    gene_name_matrix <- gene_name_list %>% do.call(cbind, .)
    ## counting equality across each gene name row
    if (sum(!apply(gene_name_matrix, 1, function(x) all(x==x[1]))) > 0){
      stop("Gene names not aligned")
    }
    # cbind to make large matrix
    counts_data_frame <- do.call(cbind, counts_data_list)
    # remove name of study from column name
    colnames(counts_data_frame) <- c(str_extract(colnames(counts_data_frame), '\\SRR\\S+'))
    #Add rownames back
    counts_data_frame <- counts_data_frame %>% as_tibble(rownames = 'gene_id')
    
    #Create counts file
    system('mkdir -p gene_counts')
    write.csv(counts_data_frame, 
              file = paste0("gene_counts/", count_file_name, ".csv"), row.names = FALSE)
    #Create mapping file
    system('mkdir -p mapping_data')
    mapping_data_frame <- bind_rows(mapping_data_list)
    write.csv(mapping_data_frame, 
              file = paste0("mapping_data/", mapping_file_name, ".csv"), row.names = FALSE)
    
    #Aggregating the data
    #Checking to make sure only data in the metadata is included in analysis
    counts_data_frame_keep <- intersect(names(counts_data_frame), metadata$run_accession)
    counts_data_frame_final <- counts_data_frame %>% select(one_of(c("gene_id", counts_data_frame_keep)))
    
    #Creating empty list to store aggregated counts in
    aggregated_counts_data_list <- list()
    #Obtain study names so data can be subset during aggregation
    aggregation_subset_data <- 
      counts_data_frame_final %>% select(-"gene_id") %>% colnames() %>% as_tibble() %>% 
      rename("value" = "run_accession") %>% 
      left_join(metadata %>% select("study_accession", "run_accession"), by= "run_accession")
    
    for (i in aggregation_subset_data$study_accession %>% unique()){
      #Subset data to ensure memory is not exhausted
      run_accessions = metadata %>% filter(study_accession == i) %>% pull(run_accession)
      counts_subset <- counts_data_frame_final[, c("gene_id", run_accessions)]
      #Reformatting data for aggregation
      long_counts_subset <- counts_subset %>% pivot_longer(-gene_id)
      names(long_counts_subset) <- c("gene_id", "run_accession", "value")
      #Joining metadata for aggregation
      long_counts_meta <- long_counts_subset %>% left_join(metadata %>% select(sample_accession, run_accession),
                                                           by = c('run_accession'))
      dt_long_counts <- data.table(long_counts_meta)
      long_counts_srs <- dt_long_counts[, .(value=mean(value)), by=list(gene_id, sample_accession)]
      
      aggregated_counts_data_list[[i]] <- long_counts_srs
    }
    #Bind rows to create aggregated count data frame
    aggregated_counts_data <- aggregated_counts_data_list %>% bind_rows()
    #Create aggregated counts file
    write.csv(aggregated_counts_data,
              file = paste0("gene_counts/", aggregated_count_file_name, ".csv"), row.names = FALSE)
  }

### Function to create count data for GTEX samples -----
# Since the GTEX data in recount3 is labeled differently, we will use a different function to locate and download this data

create_gtex_count_data_frames <- function(projects_vector, 
                                          count_file_name, 
                                          aggregated_count_file_name, 
                                          mapping_file_name, 
                                          metadata,
                                          empty_cache = TRUE) {
  
  if (empty_cache){
    print("Clearing recount3 cache")
    recount3_cache_rm()
  }
  #Indicating which recount library the function will pull from
  options(recount3_url = "http://duffel.rail.bio/recount3/")
  human_projects <- available_projects()
  #Subset to only include GTEX data
  gtex_data <- subset(human_projects, file_source == "gtex" & project_type == "data_sources")
  #Create an empty list to gather the count data
  counts_data_list <- list()
  #Create an empty list to gather mapping information
  mapping_data_list <- list()
  
  # creates list of count data frames
  for (i in projects_vector) {
    
    rse_gene = create_rse(subset(gtex_data, project == i), annotation = "gencode_v29")
    # cut down to only in metadata
    rse_gene <- rse_gene[,gsub('-','.', colnames(rse_gene)) %in% metadata$run_accession]
    assays(rse_gene)$counts <- transform_counts(rse_gene)
    
    rse_gene_counts <- recount::getTPM(rse_gene, 
                                       length_var = 'bp_length',
                                       mapped_var = 'recount_qc.star.all_mapped_reads') %>% 
      data.frame()
    rse_gene_counts <- rownames_to_column(rse_gene_counts, var = "gene_id")
    
    row.names(rse_gene_counts) <- rse_gene_counts[,1]
    rse_gene_counts <- rse_gene_counts[,-1]
  
    counts_data_list[[i]] <- rse_gene_counts
    #Input mapping data
    mapping_data <- colData(rse_gene) %>% as.data.frame()
    mapping_data <- mapping_data %>% mutate(across(everything(), as.character))
    mapping_data_list[[i]] <- mapping_data
  }
  
  # confirm gene name equality across each data frame
  gene_name_list <- list()
  for (i in 1:length(counts_data_list)){
    gene_name_list[[i]] <- row.names(counts_data_list[[i]])
  }
  gene_name_matrix <- gene_name_list %>% do.call(cbind, .)
  ## counting equality across each gene name row
  if (sum(!apply(gene_name_matrix, 1, function(x) all(x==x[1]))) > 0){
    stop("Gene names not aligned")
  }
  
  # cbind to make large matrix
  counts_data_frame <- do.call(cbind, counts_data_list)
  # remove name of study from column name
  #Using if-else to account for differences in the naming of the bone marrow samples included
  colnames(counts_data_frame) <- 
    ifelse(!grepl("GTEX", colnames(counts_data_frame)), c(str_split(colnames(counts_data_frame), '.+?(?=K.562)') %>% map(2) %>% unlist()),
           c(str_split(colnames(counts_data_frame), '.+?(?=GTEX)') %>% map(2) %>% unlist()))
  
  #Add rownames back
  counts_data_frame <- counts_data_frame %>% as_tibble(rownames = 'gene_id')
  
  #Create counts file
  system('mkdir -p gene_counts')
  write.csv(counts_data_frame, 
            file = paste0("gene_counts/", count_file_name, ".csv"), row.names = FALSE)
  
  #Create mapping file
  system('mkdir -p mapping_data')
  mapping_data_frame <- bind_rows(mapping_data_list)
  write.csv(mapping_data_frame, 
            file = paste0("mapping_data/", mapping_file_name, ".csv"), row.names = FALSE)
  
  #Aggregating the data
  #Checking to make sure only data in the metadata is included in analysis
  counts_data_frame_keep <- intersect(names(counts_data_frame), metadata$run_accession)
  counts_data_frame_final <- counts_data_frame %>% select(one_of(c("gene_id", counts_data_frame_keep)))
  
  #Creating empty list to store aggregated counts in
  aggregated_counts_data_list <- list()
  #Obtain study names so data can be subset during aggregation
  long_counts <- counts_data_frame_final %>% pivot_longer(-gene_id)
  names(long_counts) <- c("gene_id", "run_accession", "value")
  #Joining metadata for aggregation
  long_counts_meta <- long_counts %>% left_join(metadata %>% select(sample_accession, run_accession),
                                                by = c('run_accession'))
  dt_long_counts <- data.table(long_counts_meta)
  long_counts_srs <- dt_long_counts[, value:=mean(value), by=c("sample_accession", "gene_id")]
  
  aggregated_counts_data_list[[i]] <- long_counts_srs
  
  #Bind rows to create aggregated count data frame
  aggregated_counts_data <- aggregated_counts_data_list %>% bind_rows()
  #Create aggregated counts file
  write.csv(aggregated_counts_data,
            file = paste0("gene_counts/", aggregated_count_file_name, ".csv"), row.names = FALSE)
}
