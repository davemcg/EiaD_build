# PCA

#pca 
emeta <- data.table::fread('data/eyeIntegration22_meta_2022_11_04.01.csv') %>% as_tibble()
r3meta <- data.table::fread('mapping_data/mapping_data.csv.gz') %>% as_tibble()
# samples are columns
# genes are rows
mat <- vroom::vroom("gene_counts/gene_counts_matrix.csv.gz") %>% data.frame()
row.names(mat) <- mat[,1]
mat <- mat[,2:ncol(mat)]
mat <- as.matrix(mat)
######## remove 0 sum genes #################
mat <- mat[rowSums(mat)!=0,]
###############################################
#mat <- log2(mat + 1)





gene_count_meta <- colnames(mat) %>% 
  enframe(value = 'run_accession') %>% 
  left_join(emeta %>% as_tibble() %>% 
              mutate(run_accession = case_when(grepl('^\\d+', run_accession) ~ paste0('X',run_accession),
                                               TRUE ~ run_accession)) %>% 
              mutate(sample_accession = case_when(grepl('^\\d+', sample_accession) ~ paste0('X',sample_accession),
                                                  TRUE ~ sample_accession)) %>% 
              select(run_accession, sample_accession, Tissue, Source, Cohort, Source_details, Age,Sub_Tissue, study_accession, study_title) %>% unique())

# http://duffel.rail.bio/recount3/human/new_annotations/gene_sums/human.gene_sums.G029.gtf.gz
gene_annotations <- 
  rtracklayer::import("data/human.gene_sums.G029.gtf.gz") %>% 
  as.data.frame()


# add in gene names
row.names(mat) <- row.names(mat) %>% 
  enframe(value = 'gene_id') %>% 
  left_join(gene_annotations) %>% 
  mutate(ID = paste0(gene_name, ' (', gene_id, ')')) %>% 
  pull(ID)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(mat, gene_count_meta, design = ~ Tissue)
dds <- collapseReplicates(dds, groupby = colData(dds)$sample_accession)


# pull raw counts
raw_counts <- assay(dds, 'counts')
# get length info ('score' from recount3 gtf)
len_info <- row.names(raw_counts) %>% enframe() %>% mutate(gene_id = str_extract(value, 'ENSG\\d+.*') %>% gsub(')','',.)) %>% left_join(gene_annotations) 
# Divide each gene by transcript length
TPMlen <- apply( raw_counts, 2, function(x){ x / len_info$score } )
# Divide again the transcript length and apply 1e6 multiplier
TPM <- apply( TPMlen, 2, function(x) { x / sum(x) * 1E6} )



# only use protein coding genes

library(matrixStats)
ntop = 1000
Pvars <- rowVars(log2(TPM+1))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]
PCA <- prcomp(t(log2(TPM[select, ]+1)), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
gene_count_meta <- colnames(TPM) %>% 
  enframe(value = 'sample_accession') %>% 
  left_join(emeta %>% as_tibble() %>% 
              mutate(sample_accession = case_when(grepl('^\\d+', sample_accession) ~ paste0('X',sample_accession),
                                                  TRUE ~ sample_accession)) %>% 
              dplyr::select(sample_accession, Tissue, Source, Cohort, Source_details, Age,Sub_Tissue, study_accession, study_title) %>% unique())


dataGG = cbind(PCA$x, gene_count_meta)

# rotations for pc1 and pc2
rotations <- c(PCA$rotation[,1] %>% sort() %>% head(3) %>% names(),
               PCA$rotation[,1] %>% sort() %>% tail(3) %>% names(),
               PCA$rotation[,2] %>% sort() %>% head(3) %>% names(),
               PCA$rotation[,2] %>% sort() %>% tail(3) %>% names()) %>% 
  unique()

dataGG %>% as_tibble() %>% 
  mutate(Tissue = case_when(Tissue == 'Brain' ~ Tissue,
                            Cohort == 'Body' ~ 'Body',
                            TRUE ~ Tissue)) %>% 
  ggplot(., aes(PC1, PC2)) +
  geom_point(size=3, aes(color=Tissue, shape=Source)) +
  geom_segment(data = PCA$rotation[rotations,] %>% data.frame(), aes(x=0,y=0, xend = PC1*1000, yend = PC2*1000)) +
  ggrepel::geom_label_repel(data = PCA$rotation[rotations,] %>% 
                              as_tibble(rownames = 'Gene') %>% 
                              mutate(Gene = gsub(' \\(.*','',Gene)), 
                            aes(x=PC1*1000, y = PC2 * 1000, label = Gene)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(pals::glasbey(), pals::alphabet2()) %>% unname())


# facet by tissue
dataGG %>% 
  mutate(Tissue = case_when(Tissue == 'Brain' ~ 'Brain',
                            Cohort == 'Body' ~ 'Body',
                            TRUE ~ Tissue)) %>% 
  ggplot(., aes(PC1, PC2, color=Tissue, shape=Source)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(pals::glasbey(), pals::alphabet2()) %>% unname())+
  facet_wrap(~Tissue, scales = 'fixed')
# facetted, shape by  age
dataGG %>% 
  mutate(Tissue = case_when(Tissue == 'Brain' ~ 'Brain',
                            Cohort == 'Body' ~ 'Body',
                            TRUE ~ Tissue),
         Age = case_when(is.na(Age) ~ 'None',
                         TRUE ~ Age)) %>% 
  ggplot(., aes(PC1, PC2, color=Tissue, shape=Age)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(pals::glasbey(), pals::alphabet2()) %>% unname())+
  facet_wrap(~Tissue, scales = 'fixed')

# RPE, by source_details and Age
dataGG %>% 
  mutate(Age = case_when(is.na(Age) ~ 'None',
                         TRUE ~ Age)) %>% 
  filter(Tissue == 'RPE') %>% 
  ggplot(., aes(PC1, PC2, color=Source_details, shape=Age)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(pals::glasbey(), pals::alphabet2()) %>% unname())

# RPE, by study
dataGG %>% 
  mutate(Age = case_when(is.na(Age) ~ 'None',
                         TRUE ~ Age),
         study_title = str_wrap(study_title)) %>% 
  filter(study_accession == 'SRP273695') %>% 
  ggplot(., aes(PC1, PC2, color=study_title)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(pals::alphabet2(), pals::glasbey()) %>% unname()) +
  facet_wrap(~Tissue)

# mislabled retina/RPE????
# RPE, by study
dataGG %>% 
  mutate(Age = case_when(is.na(Age) ~ 'None',
                         TRUE ~ Age),
         study_title = str_wrap(study_title)) %>% 
  filter(study_accession == 'SRP273695') %>% 
  ggplot(., aes(PC1, PC2, color=Sub_Tissue, label = sample_accession)) +
  geom_point(size=3) + ggrepel::geom_text_repel() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(pals::alphabet2(), pals::glasbey()) %>% unname()) +
  facet_wrap(~Tissue,ncol = 1)


#####################
# PCA with only eye
######################

eye_samples <- gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat)) %>% pull(sample_accession)
mat_eye <- mat[,eye_samples]

ntop = 1000
Pvars <- rowVars(mat_eye)
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]
PCA <- prcomp(t(mat_eye[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
dataGG_eye = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                        PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                        Tissue = gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat)) %>% pull(Tissue),
                        Source = gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat)) %>% pull(Source),
                        Source_details = gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat)) %>% pull(Source_details),
                        Cohort =gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat)) %>% pull(Cohort),
                        Age = gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat)) %>% pull(Age),
                        Sub_Tissue = gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat)) %>% pull(Sub_Tissue))



dataGG_eye %>% 
  ggplot(., aes(PC1, PC2, color=Tissue, shape=Source)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(pals::glasbey(), pals::alphabet2()) %>% unname())




#####################
# PCA with only eye passing simple QC (not on first brush much different)
######################
# qc workup
r3meta <- data.table::fread('mapping_data/mapping_data.csv.gz') %>% as_tibble()
r3meta <- r3meta %>% mutate(sample_accession  = case_when(!is.na(sra.sample_acc) ~ sra.sample_acc, !is.na(sra.sample_acc.x) ~ sra.sample_acc.x)) 

r3_aggr <- r3meta %>% group_by(sample_accession) %>% summarise(intron_align_perc = mean(recount_qc.intron_sum_.), unique_mapping_perc = mean(recount_qc.star.uniquely_mapped_reads_.))

problem_samples <- r3_aggr %>% filter(intron_align_perc > 20 | unique_mapping_perc < 60) %>% pull(sample_accession)

eye_samples_qc <- gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat), !sample_accession %in% problem_samples) %>% pull(sample_accession)
mat_eye_qc <- mat[,eye_samples_qc]

ntop = 1000
Pvars <- rowVars(mat_eye)
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
                                                      length(Pvars)))]
PCA <- prcomp(t(mat_eye_qc[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
dataGG_eye_qc = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
                           PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
                           Tissue = gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat), !sample_accession %in% problem_samples) %>% pull(Tissue),
                           Source = gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat), !sample_accession %in% problem_samples) %>% pull(Source),
                           Source_details = gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat), !sample_accession %in% problem_samples) %>% pull(Source_details),
                           Cohort =gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat), !sample_accession %in% problem_samples) %>% pull(Cohort),
                           Age = gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat), !sample_accession %in% problem_samples) %>% pull(Age),
                           Sub_Tissue = gene_count_meta %>% filter(Cohort == 'Eye', sample_accession %in% colnames(mat), !sample_accession %in% problem_samples) %>% pull(Sub_Tissue))



dataGG_eye_qc %>% 
  ggplot(., aes(PC1, PC2, color=Tissue, shape=Source)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(pals::glasbey(), pals::alphabet2()) %>% unname())



