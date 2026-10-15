suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(SummarizedExperiment)
  library(HDF5Array)
  library(edgeR)
  library(rtracklayer)
  library(BiocParallel)
})

gtf_path <- '/data/mcgaugheyd/genomes/GRCh38/gencode.v47.primary_assembly.annotation.gtf.gz'
meta     <- '/home/mcgaugheyd/git/EiaD_build/data/2026_05_28_EiaD.meta.tsv.gz'

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
  summarise(run_accession          = paste0(unique(run_accession), collapse = ','),
            gtex_sra_run_accession = paste0(unique(gtex_sra_run_accession), collapse = ', '),
            gtex_accession         = paste0(unique(gtex_accession), collapse = ', '),
            sample_title           = paste0(unique(sample_title), collapse = ', ')) |> 
  ungroup()

fc_base  <- "featurecounts/"
fc_files <- paste0(fc_base, nmeta_filter$sample_accession, ".exon_counts.tsv")

# Filter metadata to existing count files
existing_mask <- file.exists(fc_files)
fc_files      <- fc_files[existing_mask]
nmeta_filter  <- nmeta_filter[existing_mask, ]

# ---------------------------------------------------------
# 2. Load exon count matrices, format coordinates & deduplicate
# ---------------------------------------------------------
# Read feature metadata from first file
fc_meta <- data.table::fread(fc_files[1]) %>% 
  dplyr::select(Geneid, Chr, Start, End, Strand, Length)

# Create pure coordinate row names without ENSE prefix (e.g. "chr1:100-200:+")
exon_coords <- paste0(fc_meta$Chr, ":", fc_meta$Start, "-", fc_meta$End, ":", fc_meta$Strand)

# Identify unique coordinate indices to eliminate duplicates
keep_idx <- !duplicated(exon_coords)

# Clean feature metadata and coordinates
fc_meta_clean     <- fc_meta[keep_idx, ]
exon_coords_clean <- exon_coords[keep_idx]

# Extract count column (col 7) and subset non-duplicate rows per sample
register(MulticoreParam(workers = 8))

counts_list <- bplapply(fc_files, function(f) {
  data.table::setDTthreads(1)
  data.table::fread(f, select = 7)[[1]][keep_idx]
})

counts_mat <- do.call(cbind, counts_list)
colnames(counts_mat) <- nmeta_filter$sample_accession
rownames(counts_mat) <- exon_coords_clean

# ---------------------------------------------------------
# 3. Load GTF Exon Annotations & align cleanly
# ---------------------------------------------------------
# Import GTF exon features, exclude NA exon IDs, and keep distinct entries
gtf_exons <- rtracklayer::import(gtf_path, feature.type = "exon") %>% 
  as_tibble() %>% 
  dplyr::filter(!is.na(exon_id)) %>% 
  dplyr::distinct(exon_id, .keep_all = TRUE)

# Left join preserves unique fc_meta_clean row order and count
exons_ordered <- fc_meta_clean %>% 
  dplyr::left_join(gtf_exons, by = c("Geneid" = "exon_id"))

# ---------------------------------------------------------
# 4. Construct SummarizedExperiment objects
# ---------------------------------------------------------
se_exon <- SummarizedExperiment(
  assays  = list(counts = counts_mat),
  rowData = as.data.frame(exons_ordered),
  colData = as.data.frame(nmeta_filter)
)

se_exon_cpm <- SummarizedExperiment(
  assays  = list(cpm = edgeR::cpm(counts_mat)),
  rowData = as.data.frame(exons_ordered),
  colData = as.data.frame(nmeta_filter)
)

# ---------------------------------------------------------
# 5. Save RDS and HDF5Array deliverables
# ---------------------------------------------------------
if (!dir.exists("counts")) dir.create("counts", recursive = TRUE)

saveRDS(se_exon,     file = 'counts/se_exon_human.20261013.rds')
saveRDS(se_exon_cpm, file = 'counts/se_exon_cpm_human.20261013.rds')

saveHDF5SummarizedExperiment(se_exon,     dir = 'counts/se_exon_human.20261013.hdf5',     replace = TRUE)
saveHDF5SummarizedExperiment(se_exon_cpm, dir = 'counts/se_exon_cpm_human.20261013.hdf5', replace = TRUE)
