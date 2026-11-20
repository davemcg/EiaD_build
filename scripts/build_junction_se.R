# ============================================================================
# Junction-level differential usage from STAR SJ.out.tab
# Includes BOTH annotated and novel junctions, with tiered filtering.
# Cross-references with Salmon-based DTU calls from the edgeR pipeline.
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

BPPARAM   <- MulticoreParam(workers = max(1L, parallel::detectCores() - 2L))
cache_dir <- "data/2026_05_28_EiaD.sj_se_h5"

# ---- 1. Parse all SJ.out.tab in parallel -------------------------------------
# SJ.out.tab columns: chr, intron_start(1-based), intron_end, strand(0/1/2),
# motif(0=non-canonical,1=GT/AG,2=CT/AC,3=GC/AG,4=CT/GC,5=AT/AC,6=GT/AT),
# annot(0=novel,1=annotated), unique_reads, multi_reads, overhang
#
# Per-sample tiers:
#   Annotated:  unique >= 4                              (motif unrestricted -
#                                                         GTF-annotated junctions
#                                                         are kept even if STAR
#                                                         flags the motif as
#                                                         non-canonical)
#   Novel:      motif in 1:6, unique >= 4,
#               overhang >= 12, strand != 0
#   Overhang >= 12 strongly enriches for real novel splicing vs mismapping

read_sj <- function(p) {
  d <- fread(p, header = FALSE, col.names = c(
    "chr","start","end","strand","motif","annot","unique","multi","overhang"))
  d <- d[unique >= 4L &
           ((annot == 1L) |
              (annot == 0L & motif >= 1L & motif <= 6L &
                 overhang >= 12L & strand != 0L))]
  d[, strand_chr := c("*","+","-")[strand + 1L]]
  d[, jid := paste0(chr, ":", start, "-", end, ":", strand_chr)]
  d[, .(jid, annot, unique)]
}

meta <- data.table::fread('data/2026_05_28_EiaD.meta.tsv.gz') %>%
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

sj_exists <- list.files(
  "/Volumes/data/projects/nei/mcgaughey/EiaD_2026/STAR_align/",
  recursive = TRUE) %>%
  grep("SJ.out.tab", ., value = TRUE)
meta <- meta %>%
  filter(sample_accession %in% gsub("\\/.*","", sj_exists))
meta$sj_path <- paste0(
  "/Volumes/data/projects/nei/mcgaughey/EiaD_2026/STAR_align/",
  meta$sample_accession, "/SJ.out.tab")

sj_list <- bplapply(meta$sj_path, read_sj, BPPARAM = BPPARAM)
names(sj_list) <- meta$sample_accession

# ---- 2. Build junction x sample count matrix ---------------------------------
# Pull annotation status from the first sample that saw each junction (it's
# determined by the GTF used in alignment, so it's consistent across samples
# as long as every STAR run used the same --sjdbGTFfile).
all_jids <- unique(unlist(lapply(sj_list, `[[`, "jid"), use.names = FALSE))
annot_dt <- rbindlist(lapply(sj_list, function(d) d[, .(jid, annot)]))
annot_dt <- unique(annot_dt, by = "jid")
setkey(annot_dt, jid)

mat <- matrix(0L, nrow = length(all_jids), ncol = nrow(meta),
              dimnames = list(all_jids, meta$sample_accession))
for (s in meta$sample_accession) {
  d <- sj_list[[s]]
  mat[d$jid, s] <- d$unique
}


rm(sj_list); gc()

annot_vec <- annot_dt[rownames(mat), annot]
message("Before filtering, ",
        sum(annot_vec == 1L), " annotated + ",
        sum(annot_vec == 0L), " novel junctions")

# ---- 3. Population-level prefilter -------------------------------------------
# Keep a junction if EITHER:
#   (1) it has >=10 unique reads in >=20 samples, OR
#   (2) it has >=10 unique reads in samples from >=2 different studies.
# Applied uniformly to annotated and novel junctions.
samp_keep <- rowSums(mat >= 5L) >= 20L

# "Has >=10 reads in a study" = >=1 sample in that study has >=10 reads.
# vapply returns an nrow(mat) x n_studies logical matrix.
study_has <- vapply(
  split(seq_len(ncol(mat)), meta$study_accession),
  function(ix) rowSums(mat[, ix, drop = FALSE] >= 5L) > 0L,
  logical(nrow(mat))
)
n_studies_with_reads <- rowSums(study_has)
study_keep <- n_studies_with_reads >= 2L

keep      <- samp_keep | study_keep
mat       <- mat[keep, , drop = FALSE]
annot_vec <- annot_vec[keep]
message(sum(annot_vec == 1L), " annotated + ",
        sum(annot_vec == 0L), " novel junctions retained")

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

txdb <- makeTxDbFromGFF(
  "/Volumes/data/genomes/GRCh38/gencode.v47.annotation.gtf.gz", format = "gtf")
gns  <- genes(txdb, single.strand.genes.only = TRUE)

# TxDb doesn't carry gene_name, so pull it from the GTF directly.
# Falls back to gene_id when gene_name is missing/empty (some non-coding genes,
# or non-Ensembl GTFs that lack the attribute entirely).
gtf_genes <- import(
  "/Volumes/data/genomes/GRCh38/gencode.v47.annotation.gtf.gz",
  feature.type = "gene")
if (!"gene_name" %in% names(mcols(gtf_genes))) {
  warning("GTF has no gene_name attribute; SYMBOL will equal GENEID")
  symbol_vec <- gtf_genes$gene_id
} else {
  symbol_vec <- ifelse(is.na(gtf_genes$gene_name) | gtf_genes$gene_name == "",
                       gtf_genes$gene_id, gtf_genes$gene_name)
}
gene_symbol <- setNames(symbol_vec, gtf_genes$gene_id)
rm(gtf_genes)

# Junction-in-gene = intron lies within a gene body, on the same strand.
# Novel junctions outside any annotated gene (potential novel transcripts) and
# junctions spanning multiple genes (potential readthroughs) are excluded here.
ov <- findOverlaps(jgr, gns, type = "within", ignore.strand = FALSE)
n_hits <- tabulate(queryHits(ov), nbins = length(jgr))
ov <- ov[n_hits[queryHits(ov)] == 1L]

j2g <- data.table(jid    = mcols(jgr)$jid[queryHits(ov)],
                  annot  = mcols(jgr)$annot[queryHits(ov)],
                  GENEID = names(gns)[subjectHits(ov)])
setkey(j2g, jid)

gene_tab   <- table(j2g$GENEID)
keep_genes <- names(gene_tab)[gene_tab >= 2L]
j2g        <- j2g[GENEID %in% keep_genes]
mat        <- mat[j2g$jid, , drop = FALSE]

n_novel_per_gene <- j2g[, .(n_novel = sum(annot == 0L),
                            n_total = .N), by = GENEID]
message(nrow(mat), " junctions across ",
        nrow(n_novel_per_gene), " multi-junction genes (",
        sum(j2g$annot == 0L), " junctions novel)")

# ---- 5. Build & save HDF5-backed SE for Shiny app + cache --------------------
# Rich rowData: SYMBOL, GENEID, genomic coords, annotation flag.
# Sort rows by SYMBOL (then chr/start) so gene queries are contiguous reads.
# `annot` is pulled from j2g (which was filtered/reordered alongside `mat`)
# rather than the now-stale `annot_vec`.
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
row_dt[is.na(SYMBOL), SYMBOL := GENEID]

setorder(row_dt, SYMBOL, chr, start, end)
mat <- mat[row_dt$jid, , drop = FALSE]

se <- SummarizedExperiment(
  assays  = list(counts = mat),
  rowData = DataFrame(
    jid    = row_dt$jid,
    SYMBOL = row_dt$SYMBOL,
    GENEID = row_dt$GENEID,
    chr    = row_dt$chr,
    start  = row_dt$start,
    end    = row_dt$end,
    strand = row_dt$strand,
    annot  = row_dt$annot,    # 1=annotated, 0=novel
    row.names = row_dt$jid
  ),
  colData = DataFrame(meta, row.names = meta$sample_accession)
)

dir.create(dirname(cache_dir), showWarnings = FALSE, recursive = TRUE)
saveHDF5SummarizedExperiment(
  se, cache_dir, replace = TRUE,
  # # Small row chunks x all samples: a single-gene query (typ. 5-50 jcts)
  # # hits 1 chunk; full-matrix scans still stream cleanly.
  # chunkdim = c(min(100L, nrow(mat)), ncol(mat)),
  # level    = 6
)
