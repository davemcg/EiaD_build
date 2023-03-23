load('../data/EiaD_pca_analysis.Rdata')

library(tidyverse)
library(matrixStats)

do_pca <- function(gene_by_sample, meta, ntop = 1000, scale = TRUE){
  # count values for gene_by_sample
  # metadata rows must match columns of gene_by_sample
  # quick checker for data mismatch
  if (!(colnames(gene_by_sample) == meta$sample_accession) %>% sum() == ncol(gene_by_sample)){
    stop("Check metadata and matrix samples name matching")
  }
  if (scale){
    gene_by_sample <- scale(gene_by_sample)
  }
  # remove RPL/RPS/MT genes from consideration 
  keep_genes <- grep('^RPS|^RPL|^MT-|^MPRS|^MPRL',row.names(gene_by_sample),value = TRUE, invert = TRUE)
  gene_by_sample <- gene_by_sample[keep_genes,]
  Pvars <- rowVars(log2(gene_by_sample+1))
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                        length(Pvars)))]
  PCA <- prcomp(t(log2(gene_by_sample[select, ]+1)), scale = FALSE)
  percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  
  dataGG = cbind(PCA$x, meta)
  gene_names = PCA$rotation %>% rownames()
  out <- list(PCA, dataGG, percentVar, select, gene_names)
  out
}

# all minus swaroop AMD
PCA_object <- do_pca(mat_all[,(meta_mat_all %>% filter(!grepl("AMD", Perturbation)) %>% pull(sample_accession))], meta_mat_all %>% filter(!grepl("AMD", Perturbation)) )
PCA <- PCA_object[[1]]


new_input_data <- data.table::fread('~/Desktop/test.csv.gz')
# turn into matrix
raw_matrix <- new_input_data[,-1] %>% as.matrix()
row_genes <- new_input_data[,1] %>% pull(1) %>% gsub('\\.\\d+','',.)

# gene name repair
### build a gene id table
gene_id_table <- row.names(mat_all) %>% enframe() %>% separate(value, c('gene_id','ensgene'), sep = ' \\(') %>% mutate(ensgene = gsub(')','',ensgene)) %>% select(-name)

overlap_with_ID <- row_genes[row_genes %in% gene_id_table$gene_id]
overlap_with_ens <- row_genes[row_genes %in% gene_id_table$ensgene]

if ((length(overlap_with_ID) < 1000) &
    (length(overlap_with_ens) < 1000)){
  stop("Failure to align gene names, please check your input matrix first column")
} else {
  # select column ID type to use
  if (length(overlap_with_ID) >
      length(overlap_with_ens)){
    column_val <- 1
  } else {
    column_val <- 2
  }
}



# scale
scaled <- scale(raw_matrix)
row.names(scaled) <- row_genes
# Only genes that match our internal genes used 
gene_universe <- gene_id_table %>% pull(column_val) 
scaled_cutdown <- scaled[gene_universe[gene_universe %in% row.names(scaled)], ]

# Code chunk to insert missing columns
not_included <- gene_universe[!gene_universe %in% intersect(row.names(scaled_cutdown), gene_universe)]
data <- matrix(nrow=length(not_included), ncol=ncol(scaled_cutdown))
row.names(data) <- not_included
colnames(data) <- colnames(scaled_cutdown)
scaled_cutdown <- rbind(scaled_cutdown, data)

# Replace NAs with 0
scaled_cutdown[is.na(scaled_cutdown)] <- 0

# Match order of gene_ids in kallisto data to columns from original PCA
## build PCA rotation gene object to ensure the input data is in the
## same order
rotation_gene_table <- PCA$rotation %>% row.names() %>% enframe() %>% separate(value, c('gene_id','ensgene'), sep = ' \\(') %>% mutate(ensgene = gsub(')','',ensgene)) %>% select(-name)
rotation_gene_vector <- rotation_gene_table %>% pull(column_val) # pull from the matched gene ID type
scaled_cutdown <- scaled_cutdown[rotation_gene_vector, ]

# rotate to get samples as rows
scaled_rotate <-  scaled_cutdown %>% t()

# project new data onto the PCA space
projected_input <- log2(scaled_rotate + 1) %*% PCA$rotation %>% as.data.frame()


########## merge all data
pca_projected_merge <- bind_rows(PCA_object[[2]] %>% mutate(Data = 'EiaD'),
                                 projected_input %>% mutate(Data = input$user_given_input_project_name))


# facet plot
pca_projected_merge %>% ggplot(aes(x=PC1,y=PC2)) + geom_point() + facet_wrap(~Data)
# color plot for outside
pca_projected_merge %>% filter(Cohort == 'Eye') %>% ggplot(aes(x=PC1,y=PC2, color = Data)) + geom_point() 
