library(tidyverse)
library(uwot)
library(dbscan)
library(ggrepel)
library(ggiraph)
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

eye_TPM <- tpm_gene[tpm_gene %>% pull(ID) %in% (anno %>% filter(gene_type == 'protein_coding') %>% pull(gene_name) %>% unique()),unique(eye_samples)] # only keep eye samples and protein coding genes

#eye_TPM <- tpm_gene[,unique(eye_samples)] # only keep eye samples ALL "genes"


dimRed_out <- Rtsne::Rtsne(as.matrix(log2(t(eye_TPM)+1)),perplexity = 50, check_duplicates = FALSE, theta=0.0)
#dimRed_out <- dimRed(as.matrix(log2(t(eye_TPM)+1)), n_neighbors = 20, min_dist = 0.001, pca = 50)
#dimRed_out <- prcomp(as.matrix(log2(t(eye_TPM)+1)))

dimRed_plot <- data.frame(dimRed_out$Y)
dimRed_plot$sample_accession <- colnames(eye_TPM)


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

dimRed_plot_prep_eye <- dimRed_plot_prep

# no color, no label
ggplot(dimRed_plot_prep_eye %>% 
         mutate(Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label)),aes(x=X1,y=X2)) +
  scale_shape_manual(values=c(0:20,35:50)) +
  ggtitle(paste0('Eye tissue t-sne')) +
  #geom_point(size=16, alpha=0.1, aes(colour=Tissue)) +
  geom_point(size=6, alpha=0.6, aes(shape=Origin)) +
  #geom_label_repel(aes(label=Label), alpha=0.8, size=4, box.padding = unit(0.3, "lines")) + 
  theme_minimal() + xlab('Dimension 1') + ylab('Dimension 2')

# color by tissue
ggplot(dimRed_plot_prep_eye %>% 
         mutate(Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label)),aes(x=X1,y=X2)) +
  scale_shape_manual(values=c(0:20,35:50)) +
  ggtitle(paste0('Eye tissue t-sne')) +
  geom_point(size=16, alpha=0.1, aes(colour=Tissue)) +
  geom_point(size=6, alpha=0.6, position='jitter', aes(shape=Origin)) +
  geom_label_repel(aes(label=Label), alpha=0.8, size=4, box.padding = unit(0.3, "lines")) + theme_minimal() + xlab('Dimension 1') + ylab('Dimension 2')

# color by origin
ggplot(dimRed_plot_prep_eye %>% 
         mutate(Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label)),aes(x=X1,y=X2)) +
  scale_shape_manual(values=c(0:20,35:50)) +
  ggtitle(paste0('Eye tissue t-sne')) +
  geom_point(size=16, alpha=0.1, aes(colour=Origin)) +
  geom_point(size=6, alpha=0.6, aes(shape=Origin)) +
  geom_label_repel(aes(label=Label), alpha=0.8, size=4, box.padding = unit(0.3, "lines")) + theme_minimal() + xlab('Dimension 1') + ylab('Dimension 2')

# only organoid
ggplot(dimRed_plot_prep_eye %>% 
         mutate(Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label)),
       aes(x=X1,y=X2)) +
  scale_shape_manual(values=c(0:20,35:50)) +
  ggtitle(paste0('Eye tissue t-sne')) +
  geom_point(data = dimRed_plot_prep_eye %>% filter(Origin == 'Organoid'), size=16, alpha = 0.1, color = '#00BFC4') +
  geom_point(size=6, alpha=0.6, aes(shape = Origin)) +
  geom_label_repel(aes(label=Label), alpha=0.8, size=4, box.padding = unit(0.3, "lines")) + theme_minimal() + xlab('Dimension 1') + ylab('Dimension 2')


# interactive
# color by origin
plot <- ggplot(dimRed_plot_prep_eye %>% 
                 mutate(Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label),
                        Info = paste('SRA: ', 
                                     sample_accession,
                                     '\nStudy: ', 
                                     study_title, '\n', 
                                     gsub('\\|\\|', '\n', 
                                          sample_attribute), 
                                     sep ='')),
               aes(x=X1,y=X2)) +
  scale_shape_manual(values=c(0:2,5,6,15:20)) +
  ggtitle(paste0('Eye tissue t-sne')) +
  geom_point(size=20, alpha=0.1, aes(colour=Tissue)) +
  geom_point_interactive(size=5, alpha=0.6, position = 'jitter', aes(shape=Origin, tooltip = Info)) +
  geom_label_repel(aes(label=Label), alpha=0.8, size=4, box.padding = unit(0.3, "lines")) + theme_minimal() + xlab('Dimension 1') + ylab('Dimension 2')
girafe(ggobj = plot, width_svg = 10, height_svg = 7)
