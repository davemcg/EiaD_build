# ============================================================================
# Junction-level differential usage from STAR SJ.out.tab
# Includes BOTH annotated and novel junctions, with tiered filtering.
# Generates Fractional Usage (PSI5, PSI3), Raw Counts, and CPM matrices.
# ============================================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(data.table)
  library(dplyr)
  library(BiocParallel)
  library(GenomicFeatures)
  library(txdbmaker)
  library(GenomicRanges)
  library(rtracklayer)
  library(SummarizedExperiment)
  library(HDF5Array)
})

# Cross-platform parallel back-end handling (SnowParam for Windows / MulticoreParam for Unix)
BPPARAM <- if (.Platform$OS.type == "windows") {
  BiocParallel::SnowParam(workers = max(1L, parallel::detectCores() - 2L))
} else {
  #BiocParallel::MulticoreParam(workers = max(1L, parallel::detectCores() - 2L))
  BiocParallel::MulticoreParam(workers = 16)
}

cache_dir <- "counts/se_sj_human.20261007.hdf5"

# ---- 1. Parse all SJ.out.tab in parallel -------------------------------------
read_sj <- function(p) {
  if (!file.exists(p) || file.info(p)$size == 0) {
    return(data.table(jid = character(), annot = integer(), unique = integer()))
  }
  
  tryCatch({
    d <- fread(p, header = FALSE, col.names = c(
      "chr","start","end","strand","motif","annot","unique","multi","overhang"))
    d <- d[unique >= 4L &
             ((annot == 1L) |
                (annot == 0L & motif >= 1L & motif <= 6L &
                   overhang >= 12L & strand != 0L))]
    d[, strand_chr := c("*","+","-")[strand + 1L]]
    d[, jid := paste0(chr, ":", start, "-", end, ":", strand_chr)]
    d[, .(jid, annot, unique)]
  }, error = function(e) {
    warning(paste("Error processing file:", p, "-", e$message))
    return(data.table(jid = character(), annot = integer(), unique = integer()))
  })
}

meta <- data.table::fread('~/git/EiaD_build/data/2026_05_28_EiaD.meta.tsv.gz') %>%
  filter(Organism == 'Homo sapiens')

meta <- meta |> 
  group_by(across(!c(run_accession, 
                     gtex_sra_run_accession, 
                     gtex_accession,
                     sample_title))) |> 
  summarise(run_accession = paste0(unique(run_accession), collapse = ','),
            gtex_sra_run_accession = paste0(unique(gtex_sra_run_accession), collapse = ', '),
            gtex_accession =  paste0(unique(gtex_accession), collapse = ', '),
            sample_title  =  paste0(unique(sample_title), collapse = ', ')) |> 
  ungroup()

#star_dir <- "/Volumes/data/projects/nei/mcgaughey/EiaD_2026/STAR_align/"
star_dir <- "/data/OGVFB_BG/EiaD_2026/STAR_align/"
sj_exists <- list.files(
  star_dir,
  recursive = TRUE) %>%
  grep("SJ.out.tab", ., value = TRUE)
meta <- meta %>%
  filter(sample_accession %in% gsub("\\/.*","", sj_exists))
meta$sj_path <- paste0(
  star_dir,
  meta$sample_accession, "/SJ.out.tab")

sj_list <- bplapply(meta$sj_path, read_sj, BPPARAM = BPPARAM)
names(sj_list) <- meta$sample_accession

# ---- 2. Population-level prefilter on LONG layout (Prevents Integer Overflow) ---
large_dt <- rbindlist(sj_list, idcol = "sample_accession")
rm(sj_list); gc()

study_map <- setNames(meta$study_accession, meta$sample_accession)
large_dt[, study_accession := study_map[sample_accession]]

message("Calculating population filters across long table layout...")
junction_stats <- large_dt[unique >= 5L, .(
  n_samples = .N,
  n_studies = uniqueN(study_accession)
), by = jid]

keep_jids <- junction_stats[n_samples >= 20L | n_studies >= 2L, jid]
large_dt <- large_dt[jid %in% keep_jids]

annot_dt <- unique(large_dt[, .(jid, annot)], by = "jid")
setkey(annot_dt, jid)

# ---- 3. Build Raw Count Matrix via Integer Coordinates ----------------------
unique_jids    <- annot_dt$jid
unique_samples <- meta$sample_accession

mat <- matrix(0L, 
              nrow = length(unique_jids), 
              ncol = length(unique_samples),
              dimnames = list(unique_jids, unique_samples))

row_idx <- match(large_dt$jid, unique_jids)
col_idx <- match(large_dt$sample_accession, unique_samples)
mat[cbind(row_idx, col_idx)] <- large_dt$unique

rm(large_dt, junction_stats, row_idx, col_idx); gc()

annot_vec <- annot_dt[rownames(mat), annot]
message(sum(annot_vec == 1L), " annotated + ",
        sum(annot_vec == 0L), " novel junctions retained after filters")

# ---- 4. Assign junctions to genes via the GTF --------------------------------
coord_dt <- data.table(jid = rownames(mat))
coord_dt[, c("chr","range","strand_chr") := tstrsplit(jid, ":", fixed = TRUE)]
coord_dt[, c("start","end") :=
           tstrsplit(range, "-", fixed = TRUE, type.convert = TRUE)]
jgr <- GRanges(coord_dt$chr,
               IRanges(coord_dt$start, coord_dt$end),
               strand = coord_dt$strand_chr)
mcols(jgr)$jid   <- coord_dt$jid
mcols(jgr)$annot <- annot_vec

gtf <- "/data/mcgaugheyd/genomes/GRCh38/gencode.v47.primary_assembly.annotation.gtf.gz"
txdb <- makeTxDbFromGFF(gtf, format = "gtf")
gns  <- genes(txdb, single.strand.genes.only = TRUE)

gtf_genes <- import(gtf, feature.type = "gene")

if (!"gene_name" %in% names(mcols(gtf_genes))) {
  warning("GTF has no gene_name attribute; SYMBOL will equal GENEID")
  symbol_vec <- gtf_genes$gene_id
} else {
  symbol_vec <- ifelse(is.na(gtf_genes$gene_name) | gtf_genes$gene_name == "",
                       gtf_genes$gene_id, gtf_genes$gene_name)
}

# Lookup matching both versioned (ENSG...11) and unversioned (ENSG...) IDs
gene_symbol <- c(
  setNames(symbol_vec, gtf_genes$gene_id),
  setNames(symbol_vec, sub("\\..*", "", gtf_genes$gene_id))
)
rm(gtf_genes)

ov <- findOverlaps(jgr, gns, type = "any", ignore.strand = TRUE)
n_hits <- tabulate(queryHits(ov), nbins = length(jgr))
ov <- ov[n_hits[queryHits(ov)] == 1L]

j2g <- data.table(jid    = mcols(jgr)$jid[queryHits(ov)],
                  annot  = mcols(jgr)$annot[queryHits(ov)],
                  GENEID = names(gns)[subjectHits(ov)])
setkey(j2g, jid)

gene_tab   <- table(j2g$GENEID)
keep_genes <- names(gene_tab)[gene_tab >= 2L]
j2g        <- j2g[GENEID %in% keep_genes]

# Cleanly align raw count matrix rows to filtered gene-mapped junctions
mat        <- mat[j2g$jid, , drop = FALSE]

# ---- 4b. Localized PSI (PSI5 & PSI3) & CPM Calculation ----------------------
message("Computing CPM, PSI5 (Donor), and PSI3 (Acceptor) matrices...")

# Compute CPM (Counts Per Million) using edgeR
mat_cpm <- edgeR::cpm(mat)

# Parse raw splice coordinates and assign 5' and 3' anchors based on strand
j_info <- data.table(jid = rownames(mat))
j_info[, c("chr", "range", "strand_chr") := tstrsplit(jid, ":", fixed = TRUE)]
j_info[, c("start", "end") := tstrsplit(range, "-", fixed = TRUE, type.convert = TRUE)]
j_info[, range := NULL]

# + strand: 5' site = start, 3' site = end
# - strand: 5' site = end,   3' site = start
j_info[, ss5 := ifelse(strand_chr == "+", start, end)]
j_info[, ss3 := ifelse(strand_chr == "+", end, start)]

# --- 10bp Wiggle Clustering via Fast Vector Chaining ---
setorder(j_info, chr, strand_chr, ss5)
j_info[, cluster_5 := paste0(chr, "_", strand_chr, "_5s_", 
                             cumsum(c(1L, diff(ss5) > 10L))), by = .(chr, strand_chr)]

setorder(j_info, chr, strand_chr, ss3)
j_info[, cluster_3 := paste0(chr, "_", strand_chr, "_3s_", 
                             cumsum(c(1L, diff(ss3) > 10L))), by = .(chr, strand_chr)]

# Re-align metadata tracking rows exactly back to raw matrix row order
setkey(j_info, jid)
j_info <- j_info[rownames(mat)]

# --- Compute PSI5 Matrix (Donor Splicing Ratio) ---
totals_5  <- rowsum(mat, group = j_info$cluster_5)
denom_5   <- totals_5[j_info$cluster_5, , drop = FALSE]
mat_psi5  <- mat / denom_5
mat_psi5[is.na(mat_psi5)] <- 0L
mat_psi5  <- as.integer(round(mat_psi5 * 10000))
dim(mat_psi5)      <- dim(mat)
dimnames(mat_psi5) <- dimnames(mat)

# --- Compute PSI3 Matrix (Acceptor Splicing Ratio) ---
totals_3  <- rowsum(mat, group = j_info$cluster_3)
denom_3   <- totals_3[j_info$cluster_3, , drop = FALSE]
mat_psi3  <- mat / denom_3
mat_psi3[is.na(mat_psi3)] <- 0L
mat_psi3  <- as.integer(round(mat_psi3 * 10000))
dim(mat_psi3)      <- dim(mat)
dimnames(mat_psi3) <- dimnames(mat)

message("Assays generated successfully: counts, cpm, psi5, psi3")

# ---- 5. Build & save HDF5-backed SE for Shiny app + cache --------------------
row_dt <- data.table(
  jid    = rownames(mat),
  GENEID = j2g$GENEID,
  SYMBOL = unname(gene_symbol[j2g$GENEID]),
  annot  = j2g$annot
)
row_dt[, c("chr","range","strand") := tstrsplit(jid, ":", fixed = TRUE)]
row_dt[, c("start","end") :=
         tstrsplit(range, "-", fixed = TRUE, type.convert = TRUE)]
row_dt[, range := NULL]
row_dt[is.na(SYMBOL) | SYMBOL == "", SYMBOL := GENEID]

# Map localized splice site cluster IDs into row metadata
row_dt[, cluster_5 := j_info[jid, cluster_5]]
row_dt[, cluster_3 := j_info[jid, cluster_3]]

# Sort rows by Symbol & genomic coordinates for efficient HDF5 disk reads
setorder(row_dt, SYMBOL, chr, start, end)
mat      <- mat[row_dt$jid, , drop = FALSE]
mat_cpm  <- mat_cpm[row_dt$jid, , drop = FALSE]
mat_psi5 <- mat_psi5[row_dt$jid, , drop = FALSE]
mat_psi3 <- mat_psi3[row_dt$jid, , drop = FALSE]

se <- SummarizedExperiment(
  assays  = list(
    cpm    = mat_cpm,    # Counts Per Million assay
    counts = mat,        # Raw counts
    psi5   = mat_psi5,   # Localized Donor Assay
    psi3   = mat_psi3    # Localized Acceptor Assay
  ), 
  rowData = DataFrame(
    jid       = row_dt$jid,
    SYMBOL    = row_dt$SYMBOL,
    GENEID    = row_dt$GENEID,
    chr       = row_dt$chr,
    start     = row_dt$start,
    end       = row_dt$end,
    strand    = row_dt$strand,
    annot     = row_dt$annot,
    cluster_5 = row_dt$cluster_5, 
    cluster_3 = row_dt$cluster_3, 
    row.names = row_dt$jid
  ),
  colData = DataFrame(meta, row.names = meta$sample_accession)
)

dir.create(dirname(cache_dir), showWarnings = FALSE, recursive = TRUE)

# Save multi-assay container to HDF5
saveHDF5SummarizedExperiment(
  se, cache_dir, replace = TRUE,
  chunkdim = c(min(20L, nrow(mat_cpm)), ncol(mat_cpm)),
  level    = 5
)
