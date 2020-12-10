library(tximport)
library(dplyr)
library(data.table)
library(limma)
library(edgeR)
library(qsmooth)
library(readr)
library(tidyverse)
args=commandArgs(trailingOnly = T)
save(args, file = 'testing/diffexp_args.Rdata')
setwd(args[1])
sample_table <- read_tsv(args[2])
colnames(sample_table)[1] <- 'sample_accession'
lstpms_smoothed <- read.csv(args[3], row.names = 1)  #'results/smoothed_filtered_tpms.csv',row.names = 1)
load(args[4]) # core_tight
de_comparison_contrast_names_file <- args[5]
output_limma_object_file <- args[6]
output_list_of_df_file <- args[7]
colnames(lstpms_smoothed) <- str_replace_all(colnames(lstpms_smoothed), '^E\\.MTAB\\.', 'E-MTAB-')
#colnames(sample_table) <- c('sample_acession','run','paired','tissue','subtissue','origin')
sample_table <- filter(sample_table, sample_accession %in% colnames(lstpms_smoothed))
warning('Following Samples not in sample table but in expression matrix')
print(colnames(lstpms_smoothed)[!colnames(lstpms_smoothed) %in% sample_table$sample_accession])
keep_samples <- intersect(sample_table$sample_accession, colnames(lstpms_smoothed))
sample_table <- sample_table %>% filter(sample_accession %in% keep_samples)
lstpms_smoothed <- lstpms_smoothed[,sample_table$sample_accession]
eye_samples <- filter(sample_table, Tissue %in% c('Retina','RPE','Cornea','ESC','Lens','Choroid.Plexus'))
body_samples <- filter(sample_table, !Tissue %in% c('Retina','RPE','Cornea','ESC', 'Lens', 'RetinalEpithelium', 'Choroid.Plexus'))
set.seed(123421)
gtex_sample <- sample_n(body_samples,176)
# getting different samples than david, going to use his list
# DM: nah, let's use yours
# load('ref/david_gtex_subsamples.RData')# l> gtex_sub_samples
# gtex_sample <- filter(sample_table,sample%in%gtex_sub_samples)
gtex_sample$Tissue <- gtex_sample$Sub_Tissue <-  'Body'
deg_sample_table <- rbind(eye_samples,gtex_sample,body_samples %>% mutate(Sub_Tissue = Tissue))
deg_sample_table <- left_join(deg_sample_table, core_tight %>% select(sample_accession, mapping_rate), by = 'sample_accession')
deg_counts <- lstpms_smoothed[,deg_sample_table$sample_accession]
cs= colSums(deg_counts) == 0
if (any(cs)){
  bs = colnames(deg_counts)[cs]
  warning('following samples have 0 counts  and will be removed')
  print(bs)
  deg_counts=deg_counts[,!cs]
  deg_sample_table=deg_sample_table %>% filter(!sample_accession %in% bs)

}

Sub_Tissue <- deg_sample_table$Sub_Tissue %>% as.factor()
deg_sample_table$mapping_rate[is.na( deg_sample_table$mapping_rate)] <- 0
mapping_rate <- deg_sample_table$mapping_rate %>% as.numeric()
study <- deg_sample_table$study_accession %>% as.factor()
design_eye_and_gtex <- model.matrix(~0 + mapping_rate + Sub_Tissue )
colnames(design_eye_and_gtex) <- c('mapping_rate', levels(Sub_Tissue))

y_eye_gtex <- DGEList(deg_counts)
y_eye_gtex <- calcNormFactors(y_eye_gtex)# calcNormFactors automatically performs TMM normalization when a DGELisgt is the input
v_eye_gtex <- voom(y_eye_gtex, design_eye_and_gtex)
vfit_eye_gtex <- lmFit(v_eye_gtex, design_eye_and_gtex)



conts <- combn(unique(Sub_Tissue),2) %>%
  t() %>%
  data.table() %>%
  mutate(name=paste(V1,V2,sep='_vs_'),
         contrast=paste(V1,V2,sep='-'),
         tissue_2=strsplit(V2,'_|-')%>%lapply(function(x)x[[1]]),
         tissue_1=strsplit(V1,'_|-')%>%lapply(function(x)x[[1]]),
         type_2= strsplit(V2,'_|-')%>%lapply(function(x) ifelse(x[[1]] %in% grep("_|-",Sub_Tissue, value=T, invert = T) %>% unique(),'Adult',x[[2]])),
         type_1= strsplit(V1,'_|-')%>%lapply(function(x) ifelse(x[[1]] %in% grep("_|-",Sub_Tissue, value=T, invert = T) %>% unique(),'Adult',x[[2]])),
         cont_name= paste0(tissue_1,' (',type_1,') vs ',tissue_2,' (',type_2,')')) %>%
  filter(grepl('RPE|Cornea|Retin|Lens|ESC|RetinalEpithelium|Choroi', ignore.case=T, contrast))  # remove all notEye-notEye comparisons

de_comparison_contrast_names <- conts$contrast
names(de_comparison_contrast_names) <- gsub('\\.', ' ', conts$cont_name)
save(de_comparison_contrast_names,file=de_comparison_contrast_names_file)
# https://support.bioconductor.org/p/9228/
design.pairs <-
  function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n-1))
      for (j in (i+1):n) {
        k <- k+1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i],"-",levels[j],sep="")
      }
    design
  }
cont.matrix_all <- design.pairs(c('mapping_rate', conts$V1,conts$V2) %>% unique())
cont.matrix_all <- cont.matrix_all[,grep('RPE|Cornea|Retin|Lens|ESC|RetinalEpithelium|Choroi', ignore.case=T, colnames(cont.matrix_all))] 


vfit_all <- lmFit(v_eye_gtex, design_eye_and_gtex)
cont.matrix_all <- cont.matrix_all[colnames(vfit_all),] # reorder rows to match vfit_all
vfit_all <- contrasts.fit(vfit_all, contrasts=cont.matrix_all)
efit_all <- eBayes(vfit_all)

# turn limma data into a list of dataframes
limma_de_data = list()
for (i in colnames(efit_all)){
  limma_de_data[[i]] <- topTable(efit_all, coef=i, adjust.method = 'fdr', number=300000)
}
#names(limma_de_data) <- gsub('-','_vs_', names(limma_de_data))

save(efit_all, file = output_limma_object_file)
save(limma_de_data, file = output_list_of_df_file)
