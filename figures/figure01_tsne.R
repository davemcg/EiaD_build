library(tidyverse)
library(uwot)
library(dbscan)
library(ggrepel)
library(cowplot)
load('/Volumes/data/projects/nei/mcgaughey/auto_eyeIntegration/results/core_tight.Rdata')
core_tight <- core_tight %>% select(-run_accession) %>% unique()
core_tight <- core_tight %>% mutate(Origin = case_when(Origin == 'Tissue' ~ 'Adult Tissue', TRUE ~ Origin))
tpm_gene <- read_csv('/Volumes/data/projects/nei/mcgaughey/auto_eyeIntegration/results/smoothed_filtered_tpms_gene.csv')
load('/Volumes/data/projects/nei/mcgaughey/auto_eyeIntegration/results/gene_tx_gtf_info.Rdata')
set.seed(23452345)
eye_samples <- core_tight %>%
  filter(sample_accession %in% colnames(tpm_gene)) %>%
  filter(!Tissue %in% 'ENCODE Cell Line') %>%
  filter(study_accession!='SRP012682') %>%
  filter(!sample_accession %in% c('SRS523795','SRS360124','SRS360123')) %>%
  .[['sample_accession']]

kept_samples <- core_tight %>%
  filter(Kept == 'Kept') %>% 
  pull(sample_accession)


# only keep eye samples and protein coding genes
eye_TPM <- tpm_gene[tpm_gene %>% pull(ID) %in% (anno %>% filter(gene_type == 'protein_coding') %>% pull(gene_name) %>% unique()),unique(eye_samples)] 
# all that pass QC
TPM <- tpm_gene[tpm_gene %>% pull(ID) %in% (anno %>% filter(gene_type == 'protein_coding') %>% pull(gene_name) %>% unique()), kept_samples] 

plot_maker <- function(df, algorithm = 'pca', dbscan_eps = 2.5){
  # keep only 3000 most variable genes
  variance <- apply(df, 1, var, na.rm=TRUE) 
  df$variance <- variance
  high_var_TPM <- df %>% arrange(variance) %>% tail(3000) %>% select(-variance)
   
  variance <- apply(df, 1, var, na.rm=TRUE) 
  df$variance <- variance
  high_var_TPM <- df %>% arrange(variance) %>% tail(3000) %>% select(-variance)

  if (algorithm == 'pca'){
    dimRed_out <- prcomp(as.matrix(log2(t(high_var_TPM)+1)))
    dimRed_plot <- data.frame(dimRed_out$x)
    col1 = 'PC1'
    col2 = 'PC2'
  } else if (algorithm == 'tsne'){
    dimRed_out <- Rtsne::Rtsne(as.matrix(log2(t(high_var_TPM)+1)),perplexity = 50, check_duplicates = FALSE, theta=0.0)
    dimRed_plot <- data.frame(dimRed_out$Y)
    col1 = 'X1'
    col2 = 'X2'
  } else {
    stop('Algorithm?')
  }
  dimRed_plot$sample_accession <- colnames(high_var_TPM)
  
  
  # dbscan
  dbscan_cluster <- dbscan(dimRed_plot %>% dplyr::select(!!sym(col1), !!sym(col2)), eps=dbscan_eps, minPts = 3)
  dimRed_plot$Cluster <- dbscan_cluster$cluster
  
  # create label points for each cluster
  cluster_centers <- dimRed_plot %>%
    left_join(.,core_tight) %>% group_by(Cluster) %>%
    summarise(C1=mean(!!sym(col1)),C2=mean(!!sym(col2)),Tissue=paste(unique(Tissue),collapse=','))
  
  # samples closest to center in each cluster
  center_samples <- dimRed_plot %>% left_join(.,core_tight)  %>%
    left_join(.,cluster_centers, by=c('Cluster')) %>%
    mutate(Distance = (!!sym(col1)-C1)+(!!sym(col2)-C2)) %>%
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
  
  dimRed_plot_prep %>% filter(sample_accession != 'variance')
  colnames(dimRed_plot_prep)[1:2] <- c('Dimension 1', 'Dimension 2')
  dimRed_plot_prep
}

dimRed_plot_prep_eye_pca <- plot_maker(eye_TPM, 'pca', 12)
dimRed_plot_prep_eye_tsne <- plot_maker(eye_TPM, 'tsne', 2.5)

dimRed_plot_prep_pca <- plot_maker(TPM, 'pca', 9)
dimRed_plot_prep_tsne <- plot_maker(TPM, 'tsne', 2.5)

# eye PCA
eye_pca_plot <- ggplot(dimRed_plot_prep_eye_pca %>% 
         mutate(Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label)),
       aes(x=`Dimension 1`,y=`Dimension 2`)) +
  scale_shape_manual(values=c(0:20,35:50)) +
  ggtitle(paste0('Eye Tissue PCA')) +
  geom_point(size=16, alpha=0.1, aes(colour=Tissue)) +
  geom_point(size=6, alpha=0.6, position='jitter', aes(shape=Origin)) +
  geom_label_repel(aes(label=Label), alpha=0.8, size=4, box.padding = unit(0.3, "lines")) + theme_minimal() +
  xlab('Dimension 1') + ylab('Dimension 2')
# eye tSNE
eye_tsne_plot <- ggplot(dimRed_plot_prep_eye_tsne %>% 
                         mutate(Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label)),
                       aes(x=`Dimension 1`,y=`Dimension 2`)) +
  scale_shape_manual(values=c(0:20,35:50)) +
  ggtitle(paste0('Eye Tissue t-SNE')) +
  geom_point(size=16, alpha=0.1, aes(colour=Tissue)) +
  geom_point(size=6, alpha=0.6, position='jitter', aes(shape=Origin)) +
  geom_label_repel(aes(label=Label), alpha=0.8, size=4, box.padding = unit(0.3, "lines")) + theme_minimal()
plot_grid(eye_pca_plot, eye_tsne_plot)

# pca
pca_plot <- ggplot(dimRed_plot_prep_pca %>% 
                     mutate(Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label)),
                   aes(x=`Dimension 1`,y=`Dimension 2`)) +
  scale_shape_manual(values=c(0:20,35:50)) +
  ggtitle(paste0('All Tissue PCA')) +
  geom_point(size=16, alpha=0.1, aes(colour=Tissue)) +
  geom_point(size=6, alpha=0.6, position='jitter', aes(shape=Origin)) +
  geom_label_repel(aes(label=Label), alpha=0.8, size=4, box.padding = unit(0.3, "lines")) + theme_minimal() +
  xlab('Dimension 1') + ylab('Dimension 2')

# tSNE
tsne_plot <- ggplot(dimRed_plot_prep_tsne %>% 
                         mutate(Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label)),
                       aes(x=`Dimension 1`,y=`Dimension 2`)) +
  scale_shape_manual(values=c(0:20,35:50)) +
  ggtitle(paste0('All Tissue PCA')) +
  geom_point(size=16, alpha=0.1, aes(colour=Tissue)) +
  geom_point(size=6, alpha=0.6, position='jitter', aes(shape=Origin)) +
  geom_label_repel(aes(label=Label), alpha=0.8, size=4, box.padding = unit(0.3, "lines")) + theme_minimal() +
  xlab('Dimension 1') + ylab('Dimension 2')
