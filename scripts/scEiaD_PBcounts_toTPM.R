# script to turn pseudo bulk (pb) cell type matrix from plae/scEiaD into a TPM matrix
library(tidyverse)
# get pb counts matrix
pb <- read_tsv('http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/4000-counts-universe-study_accession-scANVIprojection-15-5-20-50-0.1-CellType-Homo_sapiens.pseudoCounts.tsv.gz') 

######################################
# TPM conversion needs the GTF to calculate the gene sizes
# wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz
library(GenomicFeatures) 
# build DB
txdb <- makeTxDbFromGFF("gencode.v35.annotation.gtf.gz",format="gtf")
# aggregate by gene
exons_gene <- exonsBy(txdb,by="gene")
# calc gene sizes
gene_sizes <- as.data.frame(sum(width(reduce(exons_gene))))
# gene name (ENSG) has \.\d+ which are not used in the PB matrix
# so we erase them to facilitate the join later
row.names(gene_sizes) <- gsub('\\.\\d+','',row.names(gene_sizes))
colnames(gene_sizes)[1] <- 'length'
###########################################


# only keep pb genes that also exist gtf (because our PB matrix has mouse/macaque gene names)
pbH <- pb %>% filter(Gene %in% row.names(gene_sizes))
# get length info you calculated in the block above with makeTxDbFromGFF/etc/gene_sizes
## doing this wacky join to ensure the sizes are the in the same order as the matrix
len_info <- pbH$Gene %>% enframe() %>% left_join(gene_sizes %>% as_tibble(rownames = 'value')) 
# remove gene name col and turn into matrix
pbM <- pbH[,2:ncol(pbH)] %>% as.matrix()
# divide each gene by transcript length
TPMpblen <- apply( pbM, 2, function(x){ x / len_info$length } )
# divide again the transcript length and apply 1e6 multiplier
TPMpb <- apply( TPMpblen, 2, function(x) { x / sum(x) * 1E6} )
row.names(TPMpb) <- pbH$Gene



####################################################
# meta
scmeta <- read_tsv('http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/metadata_filter.tsv.gz')
# janky bit to match celltypes to their most common tissue (a handful of fibroblasts in the 
# retina for example...but those are not retinal data per se)
sc_tissue_ct <- scmeta %>% 
  filter(!is.na(CellType)) %>% 
  group_by(CellType, Tissue) %>% 
  summarise(Count = n()) %>% 
  mutate(Ratio = Count/sum(Count)) %>% 
  slice_max(order_by = Ratio, n = 1) %>% 
  left_join(scmeta %>% filter(!is.na(CellType)) %>% dplyr::select(Organ, Tissue) %>% unique()) %>% 
  mutate(Organ = case_when(Organ != 'Eye' ~ 'Body',
                           TRUE ~ Organ),
         CellType_ugly = gsub("\\/|-| ",'.', CellType),
         Tissue = case_when(Tissue == 'RPE-Choroid' ~ 'RPE',
                            TRUE ~ Tissue))
##################################################


