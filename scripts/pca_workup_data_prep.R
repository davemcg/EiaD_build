# PCA
library(tidyverse)
library(metamoRph)
##################
# Load in file containing samples to remove from PCA analysis
excluded_samples <- scan("data/excluded_samples.txt", what = "character")
# Load in metadata
emeta <- data.table::fread('data/eyeIntegration23_meta_2023_09_01.csv.gz') %>% 
  filter(!sample_accession %in% excluded_samples) %>% as_tibble()
##################

###################
# Load in counts
# samples are columns
# genes are rows
mat <- vroom::vroom("counts/gene_counts.csv.gz") %>% data.frame()
mat <- mat[,(!colnames(mat) %in% excluded_samples)]# Subset to remove excluded samples
mat2 <- mat[,2:ncol(mat)] %>% as.matrix()
###################


######################
# Light processing
# set rownames and remove .digit at the end
row.names(mat2) <- mat[,1] %>% str_extract(., 'ENSG\\d+')
sample_mat <- mat2
##########################


###########################
# Align metadata to counts matrix
gene_count_meta <- colnames(sample_mat) %>% 
  enframe(value = 'sample_accession') %>% 
  left_join(emeta %>% as_tibble() %>% 
              # these fields create duplicated entries for the sort-of-odd GTEX metadata
              dplyr::select(-sample_title, -gtex_sra_run_accession, -run_accession) %>% unique())
##############################




#################################
# hard set samples for PCA
# up to 50 per tissue
# the hope is not have the change this much
# but to simply project new data onto it
# set.seed(2023-08-18)
# core_set <- gene_count_meta %>% 
#   filter(!grepl("AMD", Perturbation )) %>% 
#   group_by(Tissue) %>% 
#   sample_n(50, replace = TRUE) %>% 
#   unique() %>% 
#   pull(sample_accession)

# write(core_set, file  = '../data/core_pca_samples.txt')
core_set <- scan('../data/core_pca_samples.txt', what='character')
#####################################

#####################################
# add in human readable gene names
gene_annotations <-
  data.table::fread("counts/gene_anno.csv.gz")
row.names(sample_mat) <- row.names(sample_mat) %>% 
  enframe(value = 'gene_id') %>% 
  left_join(gene_annotations %>% select(gene_id, gene_name) %>% unique() %>% 
              mutate(gene_id = gsub("\\.\\d+$","",gene_id))) %>% 
  mutate(ID = paste0(gene_name, ' (', gene_id, ')')) %>% 
  pull(ID)
####################################




#########################################
# run pca and project the rest onto it
core_mat <- sample_mat[,core_set]
gene_count_meta_core <- gene_count_meta %>% data.frame()
row.names(gene_count_meta_core) <- gene_count_meta_core$sample_accession
core_mm <- metamoRph::run_pca(core_mat, gene_count_meta_core[core_set,], irlba_n = 50)

other_mat <- sample_mat[,colnames(sample_mat)[!colnames(sample_mat) %in% core_set]]
row.names(other_mat) <- str_extract(row.names(other_mat), "ENSG\\d+")
proj_mm <- metamoRph::metamoRph(other_mat, core_mm$PCA$rotation, core_mm$center_scale)

full_pca_mat <- rbind(core_mm$PCA$x,
                      proj_mm)[,1:50]
mm_model <- metamoRph::model_build(core_mm$PCA$x, core_mm$meta$Tissue)
# guesses <- metamoRph::model_apply(mm_model, proj_mm, gene_count_meta_core[colnames(other_mat),] %>% pull(Tissue))
# guesses %>% filter(sample_label != predict) %>% 
#   left_join(gene_count_meta, by = c("sample_id" = "sample_accession")) %>% 
#   select(sample_id:max_score, Source, Age, Sub_Tissue) %>% data.frame()
#########################################



# 
# # ##############################
# # # collapse technical replicates to the sample level
# # library(DESeq2)
# # dds <- DESeqDataSetFromMatrix(mat, gene_count_meta, design = ~ Tissue)
# # dds <- collapseReplicates(dds, groupby = colData(dds)$sample_accession)
# # 
# # # now extract the sample accession level counts
# # sample_mat <- assay(dds, 'counts')
# # # merge with TPMpb (pb is pseudoBulk)
# 
# # load in scRNA pseudobulk counts data from plae/scEiaD
# pb <- read_tsv('http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/4000-counts-universe-study_accession-scANVIprojection-15-5-20-50-0.1-CellType_predict-Homo_sapiens.Staged.pseudoCounts.tsv.gz') 
# pb_genes <- pb$Gene
# pb_mat <- pb[,-1] %>% as.matrix()
# row.names(pb_mat) <- pb_genes
# 
# mat_all <- cbind(sample_mat[intersect(row.names(sample_mat) , row.names(pb_mat)),],
#                 pb_mat[intersect(row.names(sample_mat) , row.names(pb_mat)),])
# 
# #####
# # proj scrna
# plae_mm <- metamoRph::metamoRph(pb_mat, core_mm$PCA$rotation, core_mm$center_scale)
# sc_guesses <- metamoRph::model_apply(mm_model, plae_mm, str_extract(row.names(plae_mm), "__.*__") %>% gsub("__","",.))
# sc_guesses %>% filter(sample_label != predict) %>% left_join(gene_count_meta, by = c("sample_id" = "sample_accession")) %>% select(sample_id:max_score, Source, Age, Sub_Tissue) %>% data.frame()
# #######################




########################################################

######################################################
# # recreate metadata
# ## cant use previous gene_count_meta as that was built around run_accession
# ## now we have sample level counts
# ## also sc data
# meta_mat_all <- colnames(mat_all) %>% 
#   enframe(value = 'sample_accession') %>% 
#   left_join(gene_count_meta %>% dplyr::select(-name)%>% dplyr::distinct(), by = 'sample_accession') %>% 
#   # now label the scRNA pb cell types
#   mutate(CellType = case_when(grepl('__', sample_accession) ~ str_extract(sample_accession, '__.*__') %>% gsub('__','',.)),
#          Age = case_when(grepl('Matur', sample_accession) ~ 'Adult',
#                          grepl("Dev\\.$", sample_accession) ~ 'Fetal',
#                          TRUE ~ Age),
#          Source = case_when(grepl('__', sample_accession) ~ 'scRNA',
#                             TRUE ~ Source))

#########################################
save(core_mm, mm_model, full_pca_mat, proj_mm, file = 'data/EiaD_metamoRph_2023.Rdata')
system("cp data/EiaD_metamoRph_2023.Rdata ~/git/eyeIntegration_app/inst/app/www/2023/")
#write_csv(mat_all %>% as_tibble(rownames = 'Gene'), file = 'data/EiaD_pca_analysis_mat_all.csv.gz')
#write_csv(meta_mat_all, file = 'data/EiaD_pca_analysis_meta_mat_all.csv.gz')
#########################################