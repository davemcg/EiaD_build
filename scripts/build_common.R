# ============================================================================
# build_common.R -- shared helpers for the build_*_h5_se.R scripts
#
# Reads the reference block for a species from the existing pipeline config.yaml
#   references:
#     <species>:
#       annotation_path: ...
#       gtf: ...
# and takes everything else (metadata, input dir, output dir, version) as
# command line arguments.
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(dplyr)
  library(data.table)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# species key (config.yaml references:) -> value of the Organism column in the metadata
# add a line here for a new species, or pass --organism
ORGANISM <- c(human = "Homo sapiens",
              rat   = "Rattus norvegicus")

parse_eiad_args <- function(default_input_dir) {
  p <- OptionParser(
    usage = "usage: %prog --species human|rat --meta META.tsv.gz [options]",
    option_list = list(
      make_option("--species",   type = "character",
                  help = "key under references: in config.yaml (e.g. human, rat)"),
      make_option("--config",    type = "character", default = "config.yaml",
                  help = "pipeline config.yaml [default: %default]"),
      make_option("--meta",      type = "character",
                  help = "EiaD metadata tsv / tsv.gz (has Organism, sample_accession, ...)"),
      make_option("--input_dir", type = "character", default = default_input_dir,
                  help = "directory with this assay's per-sample inputs [default: %default]"),
      make_option("--out_dir",   type = "character", default = "counts",
                  help = "output directory [default: %default]"),
      make_option("--version",   type = "character", default = format(Sys.Date(), "%Y%m%d"),
                  help = "version stamp in output file names [default: today, %default]"),
      make_option("--organism",  type = "character", default = NULL,
                  help = "override the Organism value used to filter the metadata"),
      make_option("--workers",   type = "integer",
                  default = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")),
                  help = "parallel workers [default: SLURM_CPUS_PER_TASK or 8]")
    ))
  a <- parse_args(p)

  if (is.null(a$species) || is.null(a$meta)) {
    print_help(p); stop("--species and --meta are required", call. = FALSE)
  }
  if (!file.exists(a$config)) stop("config not found: ", a$config, call. = FALSE)
  if (!file.exists(a$meta))   stop("metadata not found: ", a$meta, call. = FALSE)

  ref <- yaml::read_yaml(a$config)$references[[a$species]]
  if (is.null(ref)) stop("species '", a$species, "' not found under references: in ", a$config, call. = FALSE)

  a$gtf      <- file.path(ref$annotation_path, ref$gtf)
  a$organism <- a$organism %||% unname(ORGANISM[a$species])
  if (is.null(a$organism) || is.na(a$organism))
    stop("no Organism value known for '", a$species, "'; add it to ORGANISM or pass --organism", call. = FALSE)
  if (!file.exists(a$gtf))      stop("GTF not found: ", a$gtf, call. = FALSE)
  if (!dir.exists(a$input_dir)) stop("input_dir not found: ", a$input_dir, call. = FALSE)

  dir.create(a$out_dir, showWarnings = FALSE, recursive = TRUE)
  message("species=", a$species, "  organism=", a$organism, "  version=", a$version)
  message("gtf=", a$gtf, "\ninput_dir=", a$input_dir, "  out_dir=", a$out_dir)
  a
}

# Collapse run-level metadata rows to one row per sample.
# (this block was pasted at the top of all three scripts)
load_sample_meta <- function(a) {
  m <- data.table::fread(a$meta) %>% dplyr::filter(Organism == a$organism)
  if (nrow(m) == 0) stop("no rows in metadata with Organism == '", a$organism, "'", call. = FALSE)
  m %>%
    dplyr::group_by(dplyr::across(!c(run_accession,
                                     gtex_sra_run_accession,
                                     gtex_accession,
                                     sample_title))) %>%
    dplyr::summarise(
      run_accession          = paste0(unique(run_accession), collapse = ","),
      gtex_sra_run_accession = paste0(unique(gtex_sra_run_accession), collapse = ", "),
      gtex_accession         = paste0(unique(gtex_accession), collapse = ", "),
      sample_title           = paste0(unique(sample_title), collapse = ", "),
      .groups = "drop")
}

# e.g. counts/se_gene_rat.20261120.hdf5
eiad_out <- function(a, assay, ext) {
  file.path(a$out_dir, sprintf("se_%s_%s.%s.%s", assay, a$species, a$version, ext))
}

# Transcript -> gene table, detecting the salmon target-name style so no config key is needed.
#   GENCODE (human): ENST...|ENSG...|OTTHUMG|OTTHUMT|tx_name|gene_name|len|type|   (has '|')
#   Ensembl (rat):   ENSRNOT00000012345.4                                         (no '|')
# Returns columns: target_id (salmon name), transcript_id (matches GTF transcript_id),
#                  gene_id, gene_name
make_tx2gene <- function(a, quant_file) {
  nm <- data.table::fread(quant_file, select = 1)[[1]]

  if (any(grepl("|", head(nm, 100), fixed = TRUE))) {
    parts <- data.table::tstrsplit(nm, "|", fixed = TRUE)
    t2g <- data.frame(target_id     = nm,
                      transcript_id = parts[[1]],
                      gene_id       = parts[[2]],
                      gene_name     = parts[[6]],
                      stringsAsFactors = FALSE)
  } else {
    g <- as.data.frame(rtracklayer::import(a$gtf, feature.type = "transcript"))
    if (!"gene_name" %in% names(g)) g$gene_name <- NA_character_
    g <- unique(g[, c("transcript_id", "gene_id", "gene_name")])
    g <- g[!duplicated(g$transcript_id), ]
    tid <- sub("\\.\\d+$", "", nm)                 # Ensembl cdna names carry a version, GTF ids do not
    i   <- match(tid, g$transcript_id)
    t2g <- data.frame(target_id     = nm,
                      transcript_id = tid,
                      gene_id       = g$gene_id[i],
                      gene_name     = g$gene_name[i],
                      stringsAsFactors = FALSE)
  }

  frac_missing <- mean(is.na(t2g$gene_id))
  if (frac_missing > 0) {
    warning(sprintf("%.2f%% of salmon targets have no gene_id (index and GTF may be different releases)",
                    100 * frac_missing), call. = FALSE)
    if (frac_missing > 0.05) stop("too many unmapped transcripts; check the salmon index vs the GTF", call. = FALSE)
    t2g <- t2g[!is.na(t2g$gene_id), ]
  }
  t2g
}
