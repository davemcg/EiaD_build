
setwd('/data/swamyvs/autoRNAseq')
library(broom)
library(UpSetR)
library(tidyverse)
library(limma)
library(DT)
#library(ReporteRs)
#load('data/limma_voom_DE_all_by_all_synthesticSet2.Rdata')
#load('differential_expression_data.Rdata')
source('scripts/GO_enrichment.R')
source('scripts/GO_term_finder.R')
#//change efit_vinny to efit for davids data
# for GO enrichement, big six
#background_genes <- topTable(big_six[[1]], number=30000) %>% rownames_to_column('Gene') %>% dplyr::select(Gene)#// need big 6
# for GO enrichment, all by all
load('results/diffexp_efit.Rdata')
efit_all=efit_all_vinny
background_genes_all_by_all <- topTable(efit_all, number=30000)%>% rownames_to_column('Gene') %>% dplyr::select(Gene)
#options("ReporteRs-fontsize"=8, "ReporteRs-default-font"="Monaco")
set_maker2 <- function(list.of.all){
    out = list()
    # returns setdiff of list(s) against all lists
    # example: length(setdiff(intersect(test1_comps[[1]],test1_comps[[6]]),unlist(test1_comps[c(2,3,4,5)])))
    combinations <- combn(c(1,2,3,4,5,6),2)
    for (i in 1:6){
        print(names(list.of.all)[i])
        set2 <- setdiff(1:6, i)
        out[[names(list.of.all)[i]]] <- setdiff(list.of.all[[i]], unlist(list.of.all[set2])) 
    }
    for (i in 1:ncol(combinations)){
        set1 <- combinations[, i]
        print(paste(names(list.of.all)[set1], collapse = ', '))
        set2 <- setdiff(1:6, set1)
        out[[paste(names(list.of.all)[set1], collapse = ', ')]] <- setdiff(intersect(list.of.all[[set1[1]]], list.of.all[[set1[2]]]), unlist(list.of.all[set2]))
    }  
    out
}

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
# just down vs 'control' 
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


go_printer <- function(go_frame, num_to_print, title){
    out <- go_frame %>% arrange(as.numeric(`P value`)) %>% 
        head(n=num_to_print)  %>%  FlexTable(body.cell.props = cellProperties( padding = 2 ), 
                                             header.par.props = parProperties(text.align = "left" ), 
                                             body.par.props = parProperties(text.align = 'left'), header.columns = FALSE) %>%
        addHeaderRow(value = c(title), colspan = c(8), par.properties = parProperties(text.align='left')) %>% 
        addHeaderRow(value = colnames(go_frame)) %>% 
        setFlexTableWidths(widths = c(1, 1.4, 1.2, 1, 0.5, 0.7, 0.7, 6) )  %>% 
        setFlexTableBorders(inner.vertical = no_border, inner.horizontal = no_border,
                            outer.vertical = no_border, outer.horizontal = big_border )
    print(out)
}

# alt go printer to try and make a bit more sense out of the 80 bajillion results we get for retina vs rpe
go_printer2 <- function(go_frame, num_to_print, title){
    out <- go_frame %>% filter(as.numeric(`P value (FDR)`) < 0.05, as.numeric(Count) > 9, `Odds Ratio` > 4) %>% 
        head(n=num_to_print) %>% 
        FlexTable(body.cell.props = cellProperties( padding = 2 ), 
                  header.par.props = parProperties(text.align = "left" ), 
                  body.par.props = parProperties(text.align = 'left'), header.columns = FALSE) %>%
        addHeaderRow(value = c(title), colspan = c(8), par.properties = parProperties(text.align='left')) %>% 
        addHeaderRow(value = colnames(go_frame)) %>% 
        setFlexTableWidths(widths = c(1, 1.4, 1.2, 1, 0.5, 0.7, 0.7, 6) )  %>% 
        setFlexTableBorders(inner.vertical = no_border, inner.horizontal = no_border,
                            outer.vertical = no_border, outer.horizontal = big_border )
    print(out)
}

all_go_tester <- function(test_set, test_name, background_genes){
    column_classes <- rep(x='character',11)
    column_names <- c("GO ID","P value","P value (FDR)","Odds Ratio","Expected Count","Count","Size","Term","Test", "Set", "Ontology")
    empty_frame <- read.table(text='', colClasses = column_classes, col.names = column_names)
    all_go_terms = empty_frame
    
    for (i in names(test_set)){
        if (length(test_set[[i]])<6){
            print(paste(i, 'has less than 6 genes in it, skipping GO enrichment'))
        } 
        else {
            print(i)
            for (ontology in c('BP','MF')){
                go_run <- GO_enrichment(test_set[[i]], background_genes, ontology)
                go_run <- go_run %>% data.frame() %>% mutate(Test=test_name, Set=i, Ontology=ontology)
                colnames(go_run)[1] <- 'GO.ID'
                print(colnames(go_run)==colnames(all_go_terms))
                all_go_terms <- rbind(all_go_terms, go_run)
            }
        }
    }
    
    colnames(all_go_terms) <- c("GO ID","P value","P value (FDR)","Odds Ratio","Expected Count","Count","Size","Term","Test", "Set", "Ontology")
    all_go_terms
}




up_gene_lists <- contrast_UPgene_lists(efit_all, 2)
# only keep body comparisons
# up_gene_lists <- up_gene_lists[grep('Body',names(contrast_UPgene_lists(efit_all, 2)))]
down_gene_lists <- contrast_DOWNgene_lists(efit_all, 2)
# only keep body comparisons
# down_gene_lists <- down_gene_lists[grep('Body',names(contrast_DOWNgene_lists(efit_all, 2)))]
go_up <- all_go_tester(up_gene_lists,'Up', background_genes_all_by_all)
save(go_up,file = 'results/go_up.Rdata')
go_down <- all_go_tester(down_gene_lists,'Down', background_genes_all_by_all)
save(go_down,file = 'results/go_down_david.Rdata')
all_vs_all_go <- rbind(go_up, go_down)

all_vs_all_go %>% mutate(Set = ifelse(Test=='Down', gsub('_vs_',' < ', Set),gsub('_vs_',' > ', Set) )) %>% group_by(Set, Test, Ontology) %>% filter(as.numeric(`P value (FDR)`)<0.01) %>%  summarise(Count=n())

