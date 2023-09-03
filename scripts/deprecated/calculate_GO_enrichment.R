library(tidyverse)
library(limma)
library(doParallel)

# set dir
# load arguments
args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
threads <- args[2]
limma_object <- args[3]
output <- args[4]

registerDoParallel(as.numeric(threads))
load(limma_object)


# Function to return gene list that are abs(logFC) > 1 and adj.P.Val < 0.05 for the different contrasts
## just up vs 'control' (right side of contrast)
contrast_UPgene_lists <- function(limma_data, logFC){
  vector_of_comparisons <- colnames(limma_data)
  out_list <- list()
  for (i in 1:length(vector_of_comparisons)){
    stats <-  topTable(limma_data, coef=i, number = 30000, adjust.method = 'fdr', p.value = 0.05)
    if(nrow(stats)==0){
      out_list[vector_of_comparisons[i]] <- list('')
      next}
    stats_cut_down <- stats[stats[,'logFC'] > logFC,]
    out_list[vector_of_comparisons[i]] <- list(row.names(stats_cut_down))
    }
  out_list
}
## just down vs 'control' 
contrast_DOWNgene_lists <- function(limma_data, logFC){
  vector_of_comparisons <- colnames(limma_data)
  out_list <- list()
  for (i in 1:length(vector_of_comparisons)){
    stats <-  topTable(limma_data, coef=i, number = 30000, adjust.method = 'fdr', p.value = 0.05)
    if(nrow(stats)==0){
      out_list[vector_of_comparisons[i]] <- list('')
      next}
    stats_cut_down <- stats[stats[,'logFC'] < -logFC,]
    out_list[vector_of_comparisons[i]] <- list(row.names(stats_cut_down))
    }
  out_list
}

up_gene_lists <- contrast_UPgene_lists(efit_all, 2) 
down_gene_lists <- contrast_DOWNgene_lists(efit_all, 2)

background_genes_all_by_all <- topTable(efit_all, number=300000) %>% rownames_to_column('Gene') %>% dplyr::select(Gene)

go_maker <- function(comparison){
	cat(comparison)
    suppressPackageStartupMessages({
        library(clusterProfiler)
        library(org.Hs.eg.db)
    })
ids <- bitr(background_genes_all_by_all$Gene, fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")

set_up <- bitr(up_gene_lists[[comparison]],
               fromType="SYMBOL",
               toType="ENTREZID",
               OrgDb = "org.Hs.eg.db")
set_down <- bitr(down_gene_lists[[comparison]],
               fromType="SYMBOL",
               toType="ENTREZID",
               OrgDb = "org.Hs.eg.db")

ego_up <- enrichGO(gene          = set_up$ENTREZID,
                   universe      = ids$ENTREZID %>% unique(),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   pool = TRUE,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.1,
                   keyType = 'ENTREZID',
                   readable      = TRUE)

ego_down <- enrichGO(gene          = set_down$ENTREZID,
                     universe      = ids$ENTREZID %>% unique(),
                     OrgDb         = org.Hs.eg.db,
                     ont           = "ALL",
                     pool = TRUE,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.1,
                     qvalueCutoff = 0.1,
                     keyType = 'ENTREZID',
                     readable      = TRUE)

up_and_down <- bind_rows(ego_up %>%
                             data.frame() %>%
                             mutate(geneID = gsub('/', '<br>', geneID)) %>%
                             mutate_at(.vars=vars(pvalue,p.adjust,qvalue), .funs=funs(formatC(., format = "e", digits = 2))) %>%
                             mutate(Test = 'Up', Set = comparison),
                         ego_down %>%
                             data.frame() %>%
                             mutate(geneID = gsub('/', '<br>', geneID)) %>%
                             mutate_at(.vars=vars(pvalue,p.adjust,qvalue), .funs=funs(formatC(., format = "e", digits = 2))) %>%
                             mutate(Test = 'Down', Set = comparison))

up_and_down
}

comparisons <- efit_all$contrasts %>% colnames()
all_vs_all_go <- foreach(i=1:length(comparisons), .combine = rbind) %dopar% {
	go_run <- tryCatch(go_maker(comparisons[i]), error = function(e) data.frame("ONTOLOGY" = NA, "ID" = NA, "Description" = NA, "GeneRatio" = NA, "BgRatio" = NA, "pvalue" = NA, "p.adjust" = NA, "qvalue" = NA, "geneID" = NA, "Count" = NA, "Test" = NA, Set = comparisons[i]))
}

save(all_vs_all_go, file = output)
