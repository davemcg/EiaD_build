library(tximport)
library(dplyr)
library(data.table)
library(limma)
library(edgeR)
library(qsmooth)
library(readr)
args=commandArgs(trailingOnly = T)

setwd(args[1])
sample_table <- read.table(args[2],stringsAsFactors = F, header = F,sep = '\t')
lstpms_smoothed <- read.csv(args[3], row.names = 1)  #'results/smoothed_filtered_tpms.csv',row.names = 1)
de_comparison_contrast_names_file <- args[4]
output_limma_object_file <- args[5]
output_list_of_df_file <- args[6]

colnames(sample_table) <- c('sample','run','paired','tissue','subtissue','origin')
sample_table <- filter(sample_table,sample%in%colnames(lstpms_smoothed))
eye_samples <- filter(sample_table,tissue%in%c('Retina','RPE','Cornea','ESC','Lens'))
body_samples <- filter(sample_table,!tissue%in%c('Retina','RPE','Cornea','ESC', 'Lens', 'RetinalEpithelium'))
set.seed(123421)
gtex_sample <- sample_n(body_samples,176)
# getting different samples than david, going to use his list
# DM: nah, let's use yours
# load('ref/david_gtex_subsamples.RData')# l> gtex_sub_samples
# gtex_sample <- filter(sample_table,sample%in%gtex_sub_samples)
gtex_sample$tissue <- gtex_sample$subtissue <-  'Body'
deg_sample_table <- rbind(eye_samples,gtex_sample)
deg_counts <- lstpms_smoothed[,deg_sample_table$sample]
subtissue <- deg_sample_table$subtissue%>%as.factor()
design_eye_and_gtex <- model.matrix(~0 + subtissue )
colnames(design_eye_and_gtex) <- levels(subtissue)

y_eye_gtex <- DGEList(deg_counts)
y_eye_gtex <- calcNormFactors(y_eye_gtex)
v_eye_gtex <- voom(y_eye_gtex, design_eye_and_gtex)
vfit_eye_gtex <- lmFit(v_eye_gtex, design_eye_and_gtex)



conts <- combn(unique(subtissue),2) %>%
  t() %>%
  data.table() %>%
  mutate(name=paste(V2,V1,sep='_vs_'),
         contrast=paste(V2,V1,sep='-'),
         tissue_2=strsplit(V2,'_')%>%lapply(function(x)x[[1]]),
         tissue_1=strsplit(V1,'_')%>%lapply(function(x)x[[1]]),
         type_2= strsplit(V2,'_')%>%lapply(function(x) ifelse(x[[1]]=='Body','Adult',x[[2]])),
         type_1= strsplit(V1,'_')%>%lapply(function(x) ifelse(x[[1]]=='Body','Adult',x[[2]])),
         cont_name= paste0(tissue_2,' (',type_2,') vs ',tissue_1,' (',type_1,')'))
de_comparison_contrast_names <- conts$name
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
cont.matrix_all <- design.pairs(c(conts$V1,conts$V2) %>% unique())

vfit_all <- lmFit(v_eye_gtex, design_eye_and_gtex)
cont.matrix_all <- cont.matrix_all[colnames(vfit_all),] # reorder to match vfit_all
vfit_all <- contrasts.fit(vfit_all, contrasts=cont.matrix_all)
efit_all <- eBayes(vfit_all)

# turn limma data into a list of dataframes
limma_de_data = list()
for (i in colnames(efit_all)){
  limma_de_data[[i]] <- topTable(efit_all, coef=i, adjust.method = 'fdr', number=300000)
}
names(limma_de_data) <- gsub('-','_vs_', names(limma_de_data))

save(efit_all, file = output_limma_object_file)
save(limma_de_data, file = output_list_of_df_file)
