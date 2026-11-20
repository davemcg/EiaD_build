library(tidyverse)
library(SummarizedExperiment)

meta <- data.table::fread("data/eyeIntegration23_meta_2023_09_01.built.csv.gz")
counts <- data.table::fread("counts/gene_counts.csv.gz") |> data.frame()


table(meta$sample_accession %in% colnames(counts))

table( colnames(counts) %in% meta$sample_accession)

gmeta <- meta |> 
  group_by(across(!c(run_accession, 
                     gtex_sra_run_accession, 
                     gtex_accession,
                     sample_title))) |> 
  summarise(run_accession = paste0(unique(run_accession), collapse = ','),
            gtex_sra_run_accession = paste0(unique(gtex_sra_run_accession), collapse = ', '),
            gtex_accession =  paste0(unique(gtex_accession), collapse = ', '),
            sample_title  =  paste0(unique(sample_title), collapse = ', ')) |> 
  ungroup()

# align count and meta
ocounts <- counts[,gmeta$sample_accession] |> as.matrix()
row.names(ocounts) <- counts$Gene

# acquire gene names
gene_table <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, 
                                    keys=gsub('\\.\\d+','',row.names(ocounts)),
                                    columns=c("ENSEMBL","SYMBOL", "MAP","GENENAME", "ENTREZID"), keytype="ENSEMBL") |> 
  group_by(ENSEMBL) |> 
  summarise(
    SYMBOL = paste(unique(SYMBOL), collapse = ', '),
    GENENAME = paste(unique(GENENAME), collapse = ', '),
    ENTREZID = paste(unique(ENTREZID), collapse = ', '),
    MAP = paste(unique(MAP), collapse = ', ')
  )

table(gsub("\\.\\d+","", row.names(ocounts)) == gene_table$ENSEMBL)

# cpm 
## full
eiad <- SummarizedExperiment(assays = list(cpm = scuttle::calculateCPM(ocounts)),
                             colData = DataFrame(gmeta),
                             rowData = gene_table)

saveRDS(eiad, file = 'data/2023.full.eiad.cpm.rse.rds')

## eye
eiad <- SummarizedExperiment(assays = 
                               list(cpm = scuttle::calculateCPM(ocounts[, gmeta |> filter(Cohort == 'Eye') |> pull(sample_accession)])),
                             colData = DataFrame(gmeta |> filter(Cohort == 'Eye')),
                             rowData = gene_table)

saveRDS(eiad, file = 'data/2023.eye.eiad.cpm.rse.rds')

## gtex
eiad <- SummarizedExperiment(assays = 
                               list(cpm = scuttle::calculateCPM(ocounts[, gmeta |> filter(Cohort == 'Body') |> pull(sample_accession)])),
                             colData = DataFrame(gmeta |> filter(Cohort == 'Body')),
                             rowData = gene_table)

saveRDS(eiad, file = 'data/2023.body.eiad.cpm.rse.rds')

# counts 
## full
eiad <- SummarizedExperiment(assays = list(counts = (ocounts)),
                             colData = DataFrame(gmeta),
                             rowData = gene_table)

saveRDS(eiad, file = 'data/2023.full.eiad.counts.rse.rds')

## eye
eiad <- SummarizedExperiment(assays = 
                               list(counts = (ocounts[, gmeta |> filter(Cohort == 'Eye') |> pull(sample_accession)])),
                             colData = DataFrame(gmeta |> filter(Cohort == 'Eye')),
                             rowData = gene_table)

saveRDS(eiad, file = 'data/2023.eye.eiad.counts.rse.rds')

## gtex
eiad <- SummarizedExperiment(assays = 
                               list(counts = (ocounts[, gmeta |> filter(Cohort == 'Body') |> pull(sample_accession)])),
                             colData = DataFrame(gmeta |> filter(Cohort == 'Body')),
                             rowData = gene_table)

saveRDS(eiad, file = 'data/2023.body.eiad.counts.rse.rds')

