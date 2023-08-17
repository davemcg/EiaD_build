### Function to create count data for non-GTEX samples -----

create_count_data_frames <- 
  function(recount3_library_url, 
           projects_vector, 
           TPM_matrix_file, 
           long_TPM_matrix_file, 
           count_matrix_file,
           long_count_matrix_file,
           mapping_file_name, 
           metadata,
           empty_cache = TRUE,
           GTEX = FALSE) {
    
    if (empty_cache){
      print("Clearing recount3 cache")
      recount3_cache_rm()
    }
    #Indicating which recount library the function will pull from
    options(recount3_url = recount3_library_url)
    human_projects <- available_projects()
    if (GTEX){
      gtex_data <- subset(human_projects, file_source == "gtex" & project_type == "data_sources")
    }
    #Create an empty list to gather the TPM data
    TPM_data_list <- list()
    #Create an empty list to gather the count data
    counts_data_list <- list()
    #Create an empty list to gather mapping information
    mapping_data_list <- list()
    # creates list of count data frames
    for (i in projects_vector) {
      success <- NULL
      attempt <- 1
      # while loop to make three attempts to pull the recount3 data
      while( is.null(success) && attempt <= 3 ) {
        try({
          print(paste(i, 'attempt number:', attempt))
          
          attempt <- attempt + 1
          if (GTEX){
            rse_gene = create_rse(subset(gtex_data, project == i), annotation = "gencode_v29")
          } else {
            rse_gene = create_rse_manual(i, annotation = "gencode_v29")
          }
          assays(rse_gene)$counts <- transform_counts(rse_gene)
          assays(rse_gene)$raw_counts <- compute_read_counts(rse_gene)
          # extract TPM
          rse_gene_TPM <- recount::getTPM(rse_gene, 
                                          length_var = 'bp_length',
                                          mapped_var = 'recount_qc.star.all_mapped_reads') %>% 
            data.frame()
          rse_gene_TPM <- rownames_to_column(rse_gene_TPM, var = "gene_id")
          
          row.names(rse_gene_TPM) <- rse_gene_TPM$gene_id
          rse_gene_TPM <- rse_gene_TPM[,-1, drop=FALSE]
          TPM_data_list[[i]] <- rse_gene_TPM
          # extract raw counts
          rse_gene_counts <- assays(rse_gene)$raw_counts %>% data.frame()
          rse_gene_counts <- rownames_to_column(rse_gene_counts, var = "gene_id")
          row.names(rse_gene_counts) <- rse_gene_counts$gene_id
          rse_gene_counts <- rse_gene_counts[,-1, drop=FALSE]
          counts_data_list[[i]] <- rse_gene_counts
          #extract mapping data
          mapping_data <- colData(rse_gene) %>% as.data.frame()
          mapping_data <- mapping_data %>% mutate(across(everything(), as.character))
          mapping_data_list[[i]] <- mapping_data
          success <- 'YES'
          
        })
      }
    }
    ############################################
    # confirm gene name equality across each data frame
    gene_name_list <- list()
    for (i in 1:length(TPM_data_list)){
      gene_name_list[[i]] <- row.names(TPM_data_list[[i]])
    }
    gene_name_matrix <- gene_name_list %>% do.call(cbind, .)
    ## counting equality across each gene name row
    if (sum(!apply(gene_name_matrix, 1, function(x) all(x==x[1]))) > 0){
      stop("Gene names not aligned")
    }
    ##########################################
    
    ##########################################
    # cbind to make large matrix
    TPM_data_frame <- do.call(cbind, TPM_data_list)
    count_data_frame <- do.call(cbind, counts_data_list)
    # fix column names
    if (GTEX){
      #Using if-else to account for differences in the naming of the bone marrow samples included
      colnames(TPM_data_frame) <- 
        ifelse(!grepl("GTEX", colnames(TPM_data_frame)), c(str_split(colnames(TPM_data_frame), '.+?(?=K.562)') %>% map(2) %>% unlist()),
               c(str_split(colnames(TPM_data_frame), '.+?(?=GTEX)') %>% map(2) %>% unlist()))
      colnames(count_data_frame) <- 
        ifelse(!grepl("GTEX", colnames(count_data_frame)), c(str_split(colnames(count_data_frame), '.+?(?=K.562)') %>% map(2) %>% unlist()),
               c(str_split(colnames(count_data_frame), '.+?(?=GTEX)') %>% map(2) %>% unlist()))
    } else {
      # remove name of study from column name
      colnames(TPM_data_frame) <- colnames(TPM_data_frame) %>% gsub('.*\\.','',.)
      colnames(count_data_frame) <- colnames(count_data_frame) %>% gsub('.*\\.','',.)
    }
    
    #Add rownames back
    TPM_data_frame <- TPM_data_frame %>% as_tibble(rownames = 'gene_id')
    count_data_frame <- count_data_frame %>% as_tibble(rownames = 'gene_id')
    #########################################
    

    
    #########################################
    #filter matrix to only samples in our metadata
    #Checking to make sure only data in the metadata is included in analysis
    ## edit metadata run_accession to prepend an "X" if it starts with a digit
    metadata$run_accession[grep("^\\d", metadata$run_accession)] <- paste0("X", metadata$run_accession[grep("^\\d", metadata$run_accession)])
    TPM_data_frame_keep <- intersect(names(TPM_data_frame), metadata$run_accession)
    TPM_data_frame_final <- TPM_data_frame %>% select(one_of(c("gene_id", TPM_data_frame_keep)))
    
    count_data_frame_keep <- intersect(names(count_data_frame), metadata$run_accession)
    count_data_frame_final <- count_data_frame %>% select(one_of(c("gene_id", count_data_frame_keep)))
    #Create counts file
    system('mkdir -p gene_counts')
    write_csv(TPM_data_frame_final, 
              file = paste0("gene_counts/", TPM_matrix_file, ".csv.gz"))
    write_csv(count_data_frame_final, 
              file = paste0("gene_counts/", count_matrix_file, ".csv.gz"))
    #Create mapping file
    system('mkdir -p mapping_data')
    mapping_data_frame <- bind_rows(mapping_data_list)
    write_csv(mapping_data_frame, 
              file = paste0("mapping_data/", mapping_file_name, ".csv.gz"))
    
    ######################################
    # Make long
    #Obtain study names so data can be subset during aggregation
    studies <- metadata %>% filter(run_accession %in% colnames(TPM_data_frame_final)) %>% pull(study_accession) %>% unique()
    
    # make TPM long
    long_TPM_data_list <- list()
    for (i in studies){
      #Subset data to ensure memory is not exhausted
      run_accessions = metadata %>% filter(study_accession == i, run_accession %in% colnames(TPM_data_frame_final)) %>% pull(run_accession)
      counts_subset <- TPM_data_frame_final[, c("gene_id", run_accessions)]
      #Reformatting data for aggregation
      long_counts_subset <- counts_subset %>% pivot_longer(-gene_id)
      names(long_counts_subset) <- c("gene_id", "run_accession", "value")
      #Joining metadata for aggregation
      long_counts_meta <- long_counts_subset %>% left_join(metadata %>% as_tibble() %>% select(sample_accession, run_accession) %>% unique(),
                                                           by = c('run_accession'))
      dt_long_counts <- data.table(long_counts_meta)
      long_counts_srs <- dt_long_counts[, .(value=mean(value)), by=list(gene_id, sample_accession)]
      
      long_TPM_data_list[[i]] <- long_counts_srs
    }
    
    #Bind rows to create aggregated count data frame
    long_TPM <- long_TPM_data_list %>% bind_rows() %>% unique()
    #Create aggregated counts file
    write_csv(long_TPM,
              file = paste0("gene_counts/", long_TPM_matrix_file, ".csv.gz"))
    #############################################################################
    
    
    ########################################################################################
    # make count long
    long_count_data_list <- list()
    for (i in studies){
      #Subset data to ensure memory is not exhausted
      run_accessions = metadata %>% filter(study_accession == i, run_accession %in% colnames(count_data_frame_final)) %>% pull(run_accession)
      counts_subset <- count_data_frame_final[, c("gene_id", run_accessions)]
      #Reformatting data for aggregation
      long_counts_subset <- counts_subset %>% pivot_longer(-gene_id)
      names(long_counts_subset) <- c("gene_id", "run_accession", "value")
      #Joining metadata for aggregation
      long_counts_meta <- long_counts_subset %>% left_join(metadata %>% as_tibble() %>% select(sample_accession, run_accession) %>% unique(),
                                                           by = c('run_accession'))
      dt_long_counts <- data.table(long_counts_meta)
      long_counts_srs <- dt_long_counts[, .(value=mean(value)), by=list(gene_id, sample_accession)]
      
      long_count_data_list[[i]] <- long_counts_srs
    }
    
    
    #Bind rows to create aggregated count data frame
    long_count <- long_count_data_list %>% bind_rows() %>% unique()
    #Create aggregated counts file
    write_csv(long_count,
              file = paste0("gene_counts/", long_count_matrix_file, ".csv.gz"))
    ######################################################################################
  }
