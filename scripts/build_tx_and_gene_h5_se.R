suppressPackageStartupMessages({
  library(tidyverse)
  library(tximport)
  library(data.table)
  library(SummarizedExperiment)
  library(HDF5Array)
})

gtf_path <- '/data/mcgaugheyd/genomes/GRCh38/gencode.v47.primary_assembly.annotation.gtf.gz'
meta <- '/home/mcgaugheyd/git/EiaD_build/data/2026_05_28_EiaD.meta.tsv.gz'

# ---------------------------------------------------------
# 1. Set import paths and metadata
# ---------------------------------------------------------
nmeta <- data.table::fread(meta)

nmeta_filter <- nmeta %>%
  filter(Organism == 'Homo sapiens') %>% 
  group_by(across(!c(run_accession, 
                     gtex_sra_run_accession, 
                     gtex_accession,
                     sample_title))) |> 
  summarise(run_accession = paste0(unique(run_accession), collapse = ','),
            gtex_sra_run_accession = paste0(unique(gtex_sra_run_accession), collapse = ', '),
            gtex_accession =  paste0(unique(gtex_accession), collapse = ', '),
            sample_title  =  paste0(unique(sample_title), collapse = ', ')) |> 
  ungroup()

salmon_base <- "salmon_quant/"
salmon_dirs <- paste0(salmon_base, nmeta_filter$sample_accession, "/quant.sf") %>% unique()
# remove empty dirs
salmon_dirs <- salmon_dirs[file.exists(paste0(salmon_dirs))]
# trim down meta to match what exists
nmeta_filter <- nmeta_filter %>% filter(sample_accession %in% gsub(salmon_base, "", gsub("/quant.sf", "", salmon_dirs)))

# ---------------------------------------------------------
# 2. Make tx2gene
# --------------------------------------------------------

anno <- fread(salmon_dirs[1])
anno$gene_id <- sapply(anno$Name,function(x) strsplit(x,'\\|')[[1]][2])
anno$gene_name <- sapply(anno$Name,function(x) strsplit(x,'\\|')[[1]][6])
anno_tximport <- anno %>% 
  dplyr::select(target_id = Name, gene_id)

# ---------------------------------------------------------
# 3. Read in counts 
# ---------------------------------------------------------
txi <- tximport(salmon_dirs, 
	type = "salmon", tx2gene = anno_tximport, 
	dropInfReps = TRUE, 
	importer = function(x) as.data.frame(data.table::fread(x)))
txi_tx <- tximport(salmon_dirs, 
	type = "salmon", tx2gene = anno_tximport, 
	dropInfReps = TRUE, 
	importer = function(x) as.data.frame(data.table::fread(x)), txOut = TRUE)
# ---------------------------------------------------------


# ---------------------------------------------------------
# 4. Construct SummarizedExperiment objects
# ---------------------------------------------------------
gtf_genes   <- rtracklayer::import(gtf_path, feature.type = "gene")
gtf_tx <- rtracklayer::import(gtf_path, feature.type = "transcript")

# alignment
genes_ordered <- row.names(txi$counts) %>% enframe(value = 'gene_id') %>% 
	left_join(gtf_genes %>% as_tibble()) %>% select(-name)

tx <- sapply(row.names(txi_tx$counts),function(x) strsplit(x,'\\|')[[1]][1])
tx_ordered <- tx %>% enframe(value = 'transcript_id') %>%
    left_join(gtf_tx %>% as_tibble()) %>% select(-name)

se_gene <- SummarizedExperiment(
  assays = list(
    counts    = txi$counts,
    abundance = txi$abundance,
    length    = txi$length
  ),
  rowData = as.data.frame(genes_ordered),
  colData = as.data.frame(nmeta_filter)
)

se_gene_cpm <- SummarizedExperiment(
  assays = list(
    cpm    = edgeR::cpm(txi$counts)
  ),
  rowData = as.data.frame(genes_ordered),
  colData = as.data.frame(nmeta_filter)
)


se_tx <- SummarizedExperiment(
  assays = list(
    counts    = txi_tx$counts,
    abundance = txi_tx$abundance,
    length    = txi_tx$length
  ),
  rowData = as.data.frame(tx_ordered),
  colData = as.data.frame(nmeta_filter)
)


se_tx_cpm <- SummarizedExperiment(
  assays = list(
    cpm    = edgeR::cpm(txi_tx$counts)
  ),
  rowData = as.data.frame(tx_ordered),
  colData = as.data.frame(nmeta_filter)
)

saveRDS(se_gene, file = 'counts/se_gene_human.20261013.rds')
saveRDS(se_gene_cpm, file = 'counts/se_gene_cpm_human.20261013.rds')
saveRDS(se_tx, file = 'counts/se_tx_human.20261013.rds')
saveRDS(se_tx_cpm, file = 'counts/se_tx_cpm_human.20261013.rds')

saveHDF5SummarizedExperiment(se_gene, dir = 'counts/se_gene_human.20261013.hdf5', replace = TRUE)
saveHDF5SummarizedExperiment(se_gene_cpm, dir = 'counts/se_gene_cpm_human.20261013.hdf5', replace=TRUE)
saveHDF5SummarizedExperiment(se_tx, dir = 'counts/se_tx_human.20261013.hdf5', replace = TRUE)
saveHDF5SummarizedExperiment(se_tx_cpm, dir = 'counts/se_tx_cpm_human.20261013.hdf5', replace = TRUE)
