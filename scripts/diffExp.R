#setwd('~/NIH/autoRNAseq/')
setwd('/data/swamyvs/autoRNAseq')
library(tximport)
library(dplyr)
library(limma)
library(edgeR)
library(qsmooth)
sample_table <- read.table('samplesrun_0928.tissue',stringsAsFactors = F, header = F,sep = '\t')
lstpms_smoothed <- read.csv('smoothed_filtered_tpms.csv',row.names = 1)
colnames(sample_table) <- c('sample','run','paired','tissue','subtissue','origin')
sample_table <- filter(sample_table,sample%in%colnames(lstpms_smoothed))
eye_samples <- filter(sample_table,tissue%in%c('Retina','RPE','Cornea','ESC'))
body_samples <- filter(sample_table,!tissue%in%c('Retina','RPE','Cornea','ESC'))
set.seed(123421)
gtex_sample <- sample_n(body_samples,176)
#getting different samples than david, going to use his list
load('david_gtex_subsamples.RData')# l> gtex_sub_samples
gtex_sample <- filter(sample_table,sample%in%gtex_sub_samples)
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



# combn(unique(subtissue),2) %>% 
#     t() %>% 
#     data.table() %>% 
#     mutate(name=paste(V2,V1,sep='_vs_'),
#            contrast=paste(V2,V1,sep='-'),
#            all_makeContrasts=paste(name,'=\"',contrast,'\",\n',sep='')) %>% 
#     .[['all_makeContrasts']] %>% cat()
cont.matrix_all <- makeContrasts(RPE_Cell.Line_vs_RPE_Stem.Cell.Line="RPE_Cell.Line-RPE_Stem.Cell.Line",
    ESC_Stem.Cell.Line_vs_RPE_Stem.Cell.Line="ESC_Stem.Cell.Line-RPE_Stem.Cell.Line",
    Retina_Adult.Tissue_vs_RPE_Stem.Cell.Line="Retina_Adult.Tissue-RPE_Stem.Cell.Line",
    RPE_Adult.Tissue_vs_RPE_Stem.Cell.Line="RPE_Adult.Tissue-RPE_Stem.Cell.Line",
    RPE_Fetal.Tissue_vs_RPE_Stem.Cell.Line="RPE_Fetal.Tissue-RPE_Stem.Cell.Line",
    Cornea_Adult.Tissue_vs_RPE_Stem.Cell.Line="Cornea_Adult.Tissue-RPE_Stem.Cell.Line",
    Cornea_Fetal.Tissue_vs_RPE_Stem.Cell.Line="Cornea_Fetal.Tissue-RPE_Stem.Cell.Line",
    Cornea_Cell.Line_vs_RPE_Stem.Cell.Line="Cornea_Cell.Line-RPE_Stem.Cell.Line",
    Retina_Stem.Cell.Line_vs_RPE_Stem.Cell.Line="Retina_Stem.Cell.Line-RPE_Stem.Cell.Line",
    Body_vs_RPE_Stem.Cell.Line="Body-RPE_Stem.Cell.Line",
    ESC_Stem.Cell.Line_vs_RPE_Cell.Line="ESC_Stem.Cell.Line-RPE_Cell.Line",
    Retina_Adult.Tissue_vs_RPE_Cell.Line="Retina_Adult.Tissue-RPE_Cell.Line",
    RPE_Adult.Tissue_vs_RPE_Cell.Line="RPE_Adult.Tissue-RPE_Cell.Line",
    RPE_Fetal.Tissue_vs_RPE_Cell.Line="RPE_Fetal.Tissue-RPE_Cell.Line",
    Cornea_Adult.Tissue_vs_RPE_Cell.Line="Cornea_Adult.Tissue-RPE_Cell.Line",
    Cornea_Fetal.Tissue_vs_RPE_Cell.Line="Cornea_Fetal.Tissue-RPE_Cell.Line",
    Cornea_Cell.Line_vs_RPE_Cell.Line="Cornea_Cell.Line-RPE_Cell.Line",
    Retina_Stem.Cell.Line_vs_RPE_Cell.Line="Retina_Stem.Cell.Line-RPE_Cell.Line",
    Body_vs_RPE_Cell.Line="Body-RPE_Cell.Line",
    Retina_Adult.Tissue_vs_ESC_Stem.Cell.Line="Retina_Adult.Tissue-ESC_Stem.Cell.Line",
    RPE_Adult.Tissue_vs_ESC_Stem.Cell.Line="RPE_Adult.Tissue-ESC_Stem.Cell.Line",
    RPE_Fetal.Tissue_vs_ESC_Stem.Cell.Line="RPE_Fetal.Tissue-ESC_Stem.Cell.Line",
    Cornea_Adult.Tissue_vs_ESC_Stem.Cell.Line="Cornea_Adult.Tissue-ESC_Stem.Cell.Line",
    Cornea_Fetal.Tissue_vs_ESC_Stem.Cell.Line="Cornea_Fetal.Tissue-ESC_Stem.Cell.Line",
    Cornea_Cell.Line_vs_ESC_Stem.Cell.Line="Cornea_Cell.Line-ESC_Stem.Cell.Line",
    Retina_Stem.Cell.Line_vs_ESC_Stem.Cell.Line="Retina_Stem.Cell.Line-ESC_Stem.Cell.Line",
    Body_vs_ESC_Stem.Cell.Line="Body-ESC_Stem.Cell.Line",
    RPE_Adult.Tissue_vs_Retina_Adult.Tissue="RPE_Adult.Tissue-Retina_Adult.Tissue",
    RPE_Fetal.Tissue_vs_Retina_Adult.Tissue="RPE_Fetal.Tissue-Retina_Adult.Tissue",
    Cornea_Adult.Tissue_vs_Retina_Adult.Tissue="Cornea_Adult.Tissue-Retina_Adult.Tissue",
    Cornea_Fetal.Tissue_vs_Retina_Adult.Tissue="Cornea_Fetal.Tissue-Retina_Adult.Tissue",
    Cornea_Cell.Line_vs_Retina_Adult.Tissue="Cornea_Cell.Line-Retina_Adult.Tissue",
    Retina_Stem.Cell.Line_vs_Retina_Adult.Tissue="Retina_Stem.Cell.Line-Retina_Adult.Tissue",
    Body_vs_Retina_Adult.Tissue="Body-Retina_Adult.Tissue",
    RPE_Fetal.Tissue_vs_RPE_Adult.Tissue="RPE_Fetal.Tissue-RPE_Adult.Tissue",
    Cornea_Adult.Tissue_vs_RPE_Adult.Tissue="Cornea_Adult.Tissue-RPE_Adult.Tissue",
    Cornea_Fetal.Tissue_vs_RPE_Adult.Tissue="Cornea_Fetal.Tissue-RPE_Adult.Tissue",
    Cornea_Cell.Line_vs_RPE_Adult.Tissue="Cornea_Cell.Line-RPE_Adult.Tissue",
    Retina_Stem.Cell.Line_vs_RPE_Adult.Tissue="Retina_Stem.Cell.Line-RPE_Adult.Tissue",
    Body_vs_RPE_Adult.Tissue="Body-RPE_Adult.Tissue",
    Cornea_Adult.Tissue_vs_RPE_Fetal.Tissue="Cornea_Adult.Tissue-RPE_Fetal.Tissue",
    Cornea_Fetal.Tissue_vs_RPE_Fetal.Tissue="Cornea_Fetal.Tissue-RPE_Fetal.Tissue",
    Cornea_Cell.Line_vs_RPE_Fetal.Tissue="Cornea_Cell.Line-RPE_Fetal.Tissue",
    Retina_Stem.Cell.Line_vs_RPE_Fetal.Tissue="Retina_Stem.Cell.Line-RPE_Fetal.Tissue",
    Body_vs_RPE_Fetal.Tissue="Body-RPE_Fetal.Tissue",
    Cornea_Fetal.Tissue_vs_Cornea_Adult.Tissue="Cornea_Fetal.Tissue-Cornea_Adult.Tissue",
    Cornea_Cell.Line_vs_Cornea_Adult.Tissue="Cornea_Cell.Line-Cornea_Adult.Tissue",
    Retina_Stem.Cell.Line_vs_Cornea_Adult.Tissue="Retina_Stem.Cell.Line-Cornea_Adult.Tissue",
    Body_vs_Cornea_Adult.Tissue="Body-Cornea_Adult.Tissue",
    Cornea_Cell.Line_vs_Cornea_Fetal.Tissue="Cornea_Cell.Line-Cornea_Fetal.Tissue",
    Retina_Stem.Cell.Line_vs_Cornea_Fetal.Tissue="Retina_Stem.Cell.Line-Cornea_Fetal.Tissue",
    Body_vs_Cornea_Fetal.Tissue="Body-Cornea_Fetal.Tissue",
    Retina_Stem.Cell.Line_vs_Cornea_Cell.Line="Retina_Stem.Cell.Line-Cornea_Cell.Line",
    Body_vs_Cornea_Cell.Line="Body-Cornea_Cell.Line",
    Body_vs_Retina_Stem.Cell.Line="Body-Retina_Stem.Cell.Line",
    levels=design_eye_and_gtex)






vfit_all <- lmFit(v_eye_gtex, design_eye_and_gtex)
vfit_all <- contrasts.fit(vfit_all, contrasts=cont.matrix_all)
efit_all_vinny <- eBayes(vfit_all)
DE_55_vinny <- topTableF(efit_all,number=300000, adjust.method = 'fdr') %>% data.frame()
rownames(DE_55_vinny) <- rownames(deg_counts)
load('david_diffexp.RData')
davetop55 <- topTableF(david_efit_all,number=300000, adjust.method = 'fdr') %>% data.frame()

vin_de55 <- DE_55[rownames(DE_55)%in%rownames(davetop55),]

davetop55 <- davetop55[rownames(davetop55)%in%rownames(vin_de55),]
davetop55 <- topTableF(david_efit_all,number=300000, adjust.method = 'fdr') %>% data.frame()
vin_de55 <- DE_55[rownames(DE_55)%in%rownames(davetop55),]
davetop55 <- davetop55[rownames(davetop55)%in%rownames(vin_de55),]
for(i in 1:59)print(cor(davetop55[,i],vin_de55[,i]))

