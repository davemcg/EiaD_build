library(tidyverse)
library(Rtsne)
library(dbscan)

args = commandArgs(trailingOnly=TRUE)
metadata <- args[1]
tpm <- args[2]
gtf_info <- args[3]
output <- args[4]

load(metadata)
core_tight <- core_tight %>% select(-run_accession) %>% unique()
core_tight <- core_tight %>% mutate(Origin = case_when(Origin == 'Tissue' ~ 'Adult Tissue', TRUE ~ Origin))
tpm_gene <- read_csv(tpm)
load(gtf_info)
set.seed(23452345)
kept_samples <- core_tight %>%
  filter(Kept == 'Kept') %>% 
  pull(sample_accession)

TPM <- tpm_gene[tpm_gene %>% pull(ID) %in% (anno %>% filter(gene_type == 'protein_coding') %>% pull(gene_name) %>% unique()), kept_samples] 

# keep only 3000 most variable genes
variance <- apply(TPM, 1, var, na.rm=TRUE) 
TPM$variance <- variance
high_var_TPM <- TPM %>% arrange(variance) %>% tail(3000) %>% select(-variance)

# tsne loop
tsne_list <- list()
for (i in seq(5, 10, by = 5)){
  print(paste("Perplexity", i, "running."))
  dimRed_out <- Rtsne::Rtsne(as.matrix(log2(t(high_var_TPM)+1)),perplexity = i, check_duplicates = FALSE, theta=0.0)
  #dimRed_out <- dimRed(as.matrix(log2(t(high_var_TPM)+1)), n_neighbors = 20, min_dist = 0.001, pca = 50)
  #dimRed_out <- prcomp(as.matrix(log2(t(high_var_TPM)+1)))
  
  dimRed_plot <- data.frame(dimRed_out$Y)
  dimRed_plot$sample_accession <- colnames(high_var_TPM)
  
  
  # dbscan
  dbscan_cluster <- dbscan(dimRed_plot %>% dplyr::select(X1,X2), eps=2.5, minPts = 3)
  dimRed_plot$Cluster <- dbscan_cluster$cluster
  
  # create label points for each cluster
  cluster_centers <- dimRed_plot %>%
    left_join(.,core_tight) %>% group_by(Cluster) %>%
    summarise(C1=mean(X1),C2=mean(X2),Tissue=paste(unique(Tissue),collapse=','))
  
  # samples closest to center in each cluster
  center_samples <- dimRed_plot %>% left_join(.,core_tight)  %>%
    left_join(.,cluster_centers, by=c('Cluster')) %>%
    mutate(Distance = (X1-C1)+(X2-C2)) %>%
    group_by(Cluster) %>%
    dplyr::slice(which.min(Distance)) %>%
    dplyr::slice(1) %>% 
    .[['sample_accession']]
  
  # cluster stats
  cluster_stats <- dimRed_plot %>% left_join(.,core_tight)  %>%
    mutate(Origin=factor(Origin, levels=c('Adult Tissue', 'Fetal Tissue', 'Stem Cell Line', 'Cell Line', 'Organoid'))) %>%
    mutate(Cluster = as.factor(Cluster)) %>%
    group_by(Cluster) %>%
    summarise(Cluster_Tissues = paste(unique(Tissue), collapse=', '), Cluster_Counts = paste(n(), ' samples', sep=''))
  
  # set up for ggplot
  dimRed_plot_prep <- dimRed_plot %>% left_join(.,core_tight)  %>%
    mutate(Origin=factor(Origin, levels=c('Adult Tissue', 'Fetal Tissue', 'Stem Cell Line', 'Cell Line', 'Organoid'))) %>%
    mutate(Cluster = as.factor(Cluster)) %>%
    left_join(., cluster_stats, by=c('Cluster')) %>%
    mutate(Label = paste(Cluster, Cluster_Tissues,sep=': ')) %>%
    mutate(Label = ifelse(sample_accession %in% center_samples, Label, ""))
  dimRed_plot_prep$perplexity <- i
  tsne_list[[as.character(i)]] <- dimRed_plot_prep
}

all_tsne_plot_prepped <- bind_rows(tsne_list)
save(all_tsne_plot_prepped, file = output)
# 
# 
# # interactive
# # color by origin
# plot <- ggplot(all_tsne_plot_prepped %>% filter(perplexity == 10) %>% 
#                  mutate(Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label),
#                         Info = paste('SRA: ', 
#                                      sample_accession,
#                                      '\nStudy: ', 
#                                      study_title, '\n', 
#                                      gsub('\\|\\|', '\n', 
#                                           sample_attribute), 
#                                      sep ='')),
#                aes(x=X1,y=X2)) +
#   scale_shape_manual(values=c(0:2,5,6,15:20)) +
#   ggtitle(paste0('t-SNE')) +
#   geom_point(size=20, alpha=0.1, aes(colour=Tissue)) +
#   geom_point_interactive(size=5, alpha=0.6, aes(shape=Origin, tooltip = Info)) +
#   geom_label_repel(aes(label=Label), alpha=0.8, size=4, box.padding = unit(0.3, "lines")) + theme_minimal() + xlab('Dimension 1') + ylab('Dimension 2')
# girafe(ggobj = plot, width_svg = 20, height_svg = 14)
