# ============================================================================
# Gene- and transcript-level SummarizedExperiments (counts, abundance, length, cpm)
# from salmon quant output.
#
# usage:
#   Rscript build_tx_and_gene_h5_se.R --species rat --meta META.tsv.gz \
#       [--config config.yaml] [--input_dir salmon_quant] [--out_dir counts] [--version 20261120]
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(tximport)
  library(data.table)
  library(SummarizedExperiment)
  library(HDF5Array)
  library(edgeR)
})

.this <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(if (is.na(.this)) "scripts" else dirname(.this), "build_common.R"))
a <- parse_eiad_args(default_input_dir = "salmon_quant")

# ---------------------------------------------------------
# 1. Metadata and quant files
# ---------------------------------------------------------
nmeta_filter <- load_sample_meta(a)

salmon_files <- file.path(a$input_dir, nmeta_filter$sample_accession, "quant.sf")
has_quant    <- file.exists(salmon_files)
message(sum(has_quant), " / ", length(has_quant), " samples have quant.sf")
if (!any(has_quant)) stop("no quant.sf files found under ", a$input_dir, call. = FALSE)
salmon_files <- salmon_files[has_quant]
nmeta_filter <- nmeta_filter[has_quant, ]

# ---------------------------------------------------------
# 2. tx2gene (works for GENCODE '|' headers and plain Ensembl names)
# ---------------------------------------------------------
t2g <- make_tx2gene(a, salmon_files[1])
anno_tximport <- t2g %>% dplyr::select(target_id, gene_id)

# ---------------------------------------------------------
# 3. Read in counts
# ---------------------------------------------------------
fread_importer <- function(x) as.data.frame(data.table::fread(x))

txi <- tximport(salmon_files,
                type = "salmon", tx2gene = anno_tximport,
                dropInfReps = TRUE, importer = fread_importer)
txi_tx <- tximport(salmon_files,
                   type = "salmon", tx2gene = anno_tximport,
                   dropInfReps = TRUE, importer = fread_importer, txOut = TRUE)

colnames(txi$counts) <- colnames(txi$abundance) <- colnames(txi$length) <-
  colnames(txi_tx$counts) <- colnames(txi_tx$abundance) <- colnames(txi_tx$length) <-
  nmeta_filter$sample_accession

# ---------------------------------------------------------
# 4. Row annotation from the GTF
# ---------------------------------------------------------
gtf_genes <- as.data.frame(rtracklayer::import(a$gtf, feature.type = "gene"))
if (!"gene_name" %in% names(gtf_genes)) gtf_genes$gene_name <- NA_character_
# Ensembl rat leaves gene_name blank for many genes: fall back to the gene_id
gtf_genes$gene_name <- ifelse(is.na(gtf_genes$gene_name) | gtf_genes$gene_name == "",
                              gtf_genes$gene_id, gtf_genes$gene_name)
gtf_genes <- gtf_genes[!duplicated(gtf_genes$gene_id), ]

gtf_tx <- as.data.frame(rtracklayer::import(a$gtf, feature.type = "transcript"))
gtf_tx <- gtf_tx[!duplicated(gtf_tx$transcript_id), ]

genes_ordered <- tibble(gene_id = rownames(txi$counts)) %>%
  left_join(gtf_genes, by = "gene_id")

tx_ids <- t2g$transcript_id[match(rownames(txi_tx$counts), t2g$target_id)]
tx_ordered <- tibble(transcript_id = tx_ids) %>%
  left_join(gtf_tx, by = "transcript_id")

stopifnot(nrow(genes_ordered) == nrow(txi$counts),
          nrow(tx_ordered)    == nrow(txi_tx$counts))

# ---------------------------------------------------------
# 5. Construct SummarizedExperiment objects
# ---------------------------------------------------------
col_data <- as.data.frame(nmeta_filter)
row.names(col_data) <- col_data$sample_accession

se_gene <- SummarizedExperiment(
  assays  = list(counts = txi$counts, abundance = txi$abundance, length = txi$length),
  rowData = as.data.frame(genes_ordered),
  colData = col_data)

se_gene_cpm <- SummarizedExperiment(
  assays  = list(cpm = edgeR::cpm(txi$counts)),
  rowData = as.data.frame(genes_ordered),
  colData = col_data)

se_tx <- SummarizedExperiment(
  assays  = list(counts = txi_tx$counts, abundance = txi_tx$abundance, length = txi_tx$length),
  rowData = as.data.frame(tx_ordered),
  colData = col_data)

se_tx_cpm <- SummarizedExperiment(
  assays  = list(cpm = edgeR::cpm(txi_tx$counts)),
  rowData = as.data.frame(tx_ordered),
  colData = col_data)

# ---------------------------------------------------------
# 6. Save RDS and HDF5 deliverables
# ---------------------------------------------------------
saveRDS(se_gene,     file = eiad_out(a, "gene",     "rds"))
saveRDS(se_gene_cpm, file = eiad_out(a, "gene_cpm", "rds"))
saveRDS(se_tx,       file = eiad_out(a, "tx",       "rds"))
saveRDS(se_tx_cpm,   file = eiad_out(a, "tx_cpm",   "rds"))

saveHDF5SummarizedExperiment(se_gene,     dir = eiad_out(a, "gene",     "hdf5"), replace = TRUE)
saveHDF5SummarizedExperiment(se_gene_cpm, dir = eiad_out(a, "gene_cpm", "hdf5"), replace = TRUE)
saveHDF5SummarizedExperiment(se_tx,       dir = eiad_out(a, "tx",       "hdf5"), replace = TRUE)
saveHDF5SummarizedExperiment(se_tx_cpm,   dir = eiad_out(a, "tx_cpm",   "hdf5"), replace = TRUE)

message("done: ", a$species, " gene/tx")
