# ============================================================================
# Exon-level SummarizedExperiments (counts, cpm) from featureCounts output
# (<input_dir>/<sample_accession>.exon_counts.tsv).
#
# usage:
#   Rscript build_exon_h5_se.R --species rat --meta META.tsv.gz \
#       [--config config.yaml] [--input_dir featurecounts] [--out_dir counts] [--version 20261120]
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(SummarizedExperiment)
  library(HDF5Array)
  library(edgeR)
  library(rtracklayer)
  library(BiocParallel)
})

.this <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(if (is.na(.this)) "scripts" else dirname(.this), "build_common.R"))
a <- parse_eiad_args(default_input_dir = "featurecounts")

# ---------------------------------------------------------
# 1. Metadata and count files
# ---------------------------------------------------------
nmeta_filter <- load_sample_meta(a)

fc_files <- file.path(a$input_dir, paste0(nmeta_filter$sample_accession, ".exon_counts.tsv"))
existing_mask <- file.exists(fc_files)
message(sum(existing_mask), " / ", length(existing_mask), " samples have exon count files")
if (!any(existing_mask)) stop("no *.exon_counts.tsv found under ", a$input_dir, call. = FALSE)
fc_files     <- fc_files[existing_mask]
nmeta_filter <- nmeta_filter[existing_mask, ]

# ---------------------------------------------------------
# 2. Load exon count matrices, format coordinates & deduplicate
# ---------------------------------------------------------
# feature metadata from the first file (all files share one featureCounts annotation)
fc_meta <- data.table::fread(fc_files[1]) %>%
  dplyr::select(Geneid, Chr, Start, End, Strand, Length)

# coordinate row names without the exon id prefix, e.g. "chr1:100-200:+"
exon_coords <- paste0(fc_meta$Chr, ":", fc_meta$Start, "-", fc_meta$End, ":", fc_meta$Strand)

keep_idx          <- !duplicated(exon_coords)
fc_meta_clean     <- fc_meta[keep_idx, ]
exon_coords_clean <- exon_coords[keep_idx]

register(MulticoreParam(workers = a$workers))

counts_list <- bplapply(fc_files, function(f) {
  data.table::setDTthreads(1)
  data.table::fread(f, select = 7)[[1]][keep_idx]
})

counts_mat <- do.call(cbind, counts_list)
colnames(counts_mat) <- nmeta_filter$sample_accession
rownames(counts_mat) <- exon_coords_clean

# ---------------------------------------------------------
# 3. GTF exon annotations
# ---------------------------------------------------------
gtf_exons <- rtracklayer::import(a$gtf, feature.type = "exon") %>% as_tibble()
if (!"exon_id" %in% names(gtf_exons)) stop("GTF has no exon_id attribute: ", a$gtf, call. = FALSE)
gtf_exons <- gtf_exons %>%
  dplyr::filter(!is.na(exon_id)) %>%
  dplyr::distinct(exon_id, .keep_all = TRUE)

# left join keeps fc_meta_clean row order and count
exons_ordered <- fc_meta_clean %>%
  dplyr::left_join(gtf_exons, by = c("Geneid" = "exon_id"))

frac_unmatched <- mean(is.na(exons_ordered$gene_id))
if (frac_unmatched > 0.05)
  warning(sprintf("%.1f%% of featureCounts exons not found in the GTF; was featureCounts run with this GTF?",
                  100 * frac_unmatched), call. = FALSE)
stopifnot(nrow(exons_ordered) == nrow(counts_mat))

# ---------------------------------------------------------
# 4. Construct SummarizedExperiment objects
# ---------------------------------------------------------
col_data <- as.data.frame(nmeta_filter)
row.names(col_data) <- col_data$sample_accession

se_exon <- SummarizedExperiment(
  assays  = list(counts = counts_mat),
  rowData = as.data.frame(exons_ordered),
  colData = col_data)

se_exon_cpm <- SummarizedExperiment(
  assays  = list(cpm = edgeR::cpm(counts_mat)),
  rowData = as.data.frame(exons_ordered),
  colData = col_data)

# ---------------------------------------------------------
# 5. Save RDS and HDF5Array deliverables
# ---------------------------------------------------------
saveRDS(se_exon,     file = eiad_out(a, "exon",     "rds"))
saveRDS(se_exon_cpm, file = eiad_out(a, "exon_cpm", "rds"))

saveHDF5SummarizedExperiment(se_exon,     dir = eiad_out(a, "exon",     "hdf5"), replace = TRUE)
saveHDF5SummarizedExperiment(se_exon_cpm, dir = eiad_out(a, "exon_cpm", "hdf5"), replace = TRUE)

message("done: ", a$species, " exon")
