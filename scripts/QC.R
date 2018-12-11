library(tximport)
library(tidyverse)
library(qsmooth)
library(Rtsne)
library(ggplot2)
library(dbscan)
library(edgeR)
library(RSQLite)
library(SRAdb)
args=commandArgs(trailingOnly = T)


sample_metadata = args[1]
gtf_file = args[2]
working_dir = args[3]
level = args[4] # transcript or gene level quantification
output_file = args[5]
sqlfile <- args[6]
setwd(working_dir)

# Qsmooth > remove median counts > remove lowly expressed genes > tSNE > DBSCAN

# load gtf to construct gene <-> transcript_mapping
# load sample tissue classification for qsmooth

sample_design <- read.table(sample_metadata, stringsAsFactors = F, header=F, sep = '\t')
gtf <- rtracklayer::readGFF(gtf_file) %>% dplyr::filter(type=='transcript')
anno <- gtf[,c("gene_id", "gene_name", "transcript_id")]
colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','tissue','sub-tissue','origin')
##generate gene, tx and core tight files

sra_con <- dbConnect(RSQLite::SQLite(),sqlfile)
que <- "SELECT * FROM sra WHERE sample_accession='BLANK'"
eye_tissue <- c('Retina','RPE','Cornea','ESC')
core_tight <- lapply(sample_design$sample_accession,function(x) dbGetQuery(sra_con, gsub('BLANK',x,que)))%>%do.call(rbind,.)%>%
  distinct%>%select(study_accession ,study_title,study_abstract,sample_accession,sample_attribute)%>%
  left_join(sample_design[,c(1,4,5)],by='sample_accession')
core_tight$origin <-  strsplit(core_tight[,'sub-tissue'],'_')%>%sapply(function(x)ifelse(x[[1]]%in%eye_tissue,x[[2]],'Tissue'))
save(core_tight,file = 'results/core_tight.Rdata')

genes= anno$gene_name %>% unique
save(genes, file='results/gene_names.Rdata')

tx=anno[,c("gene_name", "transcript_id")]%>%distinct%>%split(.[,2])%>%sapply(function(x) paste0(x[,1],' (',x[,2],')'))  
save(tx, file='results/tx_gene_names.Rdata')

##read in quant files
files0 <-paste0('RE_quant_files/',sample_design$sample_accession, '/quant.sf')
samplenames <- strsplit(files0,'/' )%>% sapply(function(x) x[2])
sample_design <-  filter(sample_design, sample_accession%in%samplenames)
#load('tpms.Rdata')
#txi.counts <- tximport(files=files0,tx2gene =  anno[,3:2],type = "salmon")
# load data at the transcript level or merge to the gene level
if (level == 'transcript') {
<<<<<<< HEAD
	txi.lsTPMs <- tximport(files=files0, txOut=TRUE, type = "salmon", countsFromAbundance = "lengthScaledTPM")
=======
	txi.lsTPMs_tx <- tximport(files=files0,txOut = T, type = "salmon", countsFromAbundance = "lengthScaledTPM")
	txi.lsTPMs <- tximport(files=files0, tx2gene =  anno[,3:2], type = "salmon", countsFromAbundance = "lengthScaledTPM")
	tpms_tx <- as.data.frame(txi.lsTPMs_tx$counts)
	colnames(tpms_tx) <- samplenames
>>>>>>> 6ef6baff4c91e660be0954da60dd8d1287e5aa9b
} else {
	txi.lsTPMs <- tximport(files=files0, tx2gene =  anno[,3:2], type = "salmon", countsFromAbundance = "lengthScaledTPM")
}
#save(txi,file = 'tpms.Rdata')
tpms <- as.data.frame(txi.lsTPMs$counts)
#counts <- as.data.frame(txi.counts$counts)
colnames(tpms)  <-  samplenames
# ## remove samples that failed to meet mapping %
# bad_mapping <- read.table('logs/bad_mapping',stringsAsFactors = F)
# sample_design <- sample_design[!sample_design$sample_accession%in%bad_mapping$V2,]
# samplenames <- sample_design$sample_accession
# tpms <- tpms[,samplenames]
tpms[is.na(tpms)] <- 0
keep_genes <- which(rowSums(tpms)>=ncol(tpms))# revove gene w less than an average count of .5
tpms <- tpms[keep_genes,]
sample_medians <- apply(tpms,2,function(x) median(x))

# more genes => lower medians, looks like 3 givees results most similar to david, see comparingQC.r
keep_median <- which(sample_medians>2)# criteria for removing samples with low counts
# removed <- samplenames[which(sample_medians<=50)]
# intersect(d.removed,removed)%>%length
# View(sample_design[removed,])
# View(sample_design[setdiff(removed,d.removed),])
tpms <- tpms[,keep_median]
if (level == 'transcript'){
  #to keep the samples the same between the gene and tx level, filter samples only at the gene level
  #going from cut of 1 > .5 keeps about 40k tx's
  tpms <- tpms_tx
  tpms[is.na(tpms)] <- 0
  keep_genes <- which(rowSums(tpms)>=ncol(tpms)*1)# revove gene w less than an average count of 1
  tpms <- tpms[keep_genes,keep_median]
  
}
#normalize for library size
norm <- DGEList(tpms)
norm <- calcNormFactors(norm)
norm_counts <- norm$counts
#extract scaling factor for each sample and  multiply
correction <- norm$samples %>% data.frame() %>% .[['norm.factors']]
lsTPM_librarySize <- norm_counts %*% diag(correction)
colnames(lsTPM_librarySize) <- colnames(tpms)

#quantile normalize samples
sample_design <- sample_design %>% filter(sample_accession  %in% colnames(lsTPM_librarySize))
qs <- qsmooth(object = lsTPM_librarySize,groupFactor = as.factor(sample_design$tissue))
lstpms_smoothed <- as.data.frame(qsmoothData(qs))

colnames(lstpms_smoothed) <- colnames(lsTPM_librarySize)

#cluster with tSNE, then run clustered data throught
tpms_smoothed_filtered <- lstpms_smoothed
print(dim(tpms_smoothed_filtered))
set.seed(23235)
tsne_out <- Rtsne(as.matrix(log2(t(tpms_smoothed_filtered)+1)),perplexity = 40, check_duplicates = FALSE, theta=0.0 )
tsne_plot <- data.frame(tsne_out$Y,sample_design[sample_design$sample_accession%in%colnames(tpms_smoothed_filtered),])
# ggplot(tsne_plot,aes(x=X1,y=X2, col=tissue))+
#   geom_point(size=2)+
#   theme_minimal()

#find outliers based on tsne grouping. calculate the center of each group, and look for points n standard deviations away
e_Dist <- function(p1,p2) return(sqrt(sum((p1-p2)^2)))
tsne_plot$outlier <- NA
for(t in tsne_plot$tissue){
  tsne.tissue <- filter(tsne_plot,tissue==t)
  center <- c(mean(tsne.tissue$X1),mean(tsne.tissue$X2))
  dist <- apply(tsne.tissue[,1:2],1,function(x) e_Dist(x,center))
  n=4 #number of standard deviations a point is allowed to be from the center
  allowed <- c(mean(dist)-n*sd(dist),mean(dist)+n*sd(dist))
  outliers <- dist<allowed[1] | dist> allowed[2]
  tsne_plot[tsne_plot$tissue==t,]$outlier <- outliers
}
tsne_plot$outlier[is.na(tsne_plot$outlier)] <- F
# ggplot(tsne_plot,aes(x=X1,y=X2,col=tissue, shape=outlier))+
#   geom_point(size=3)+
#   ggtitle('outlier from tSNE data')+
#   theme_minimal()
trimmed_counts_smoothed <- tpms_smoothed_filtered[,!tsne_plot$outlier]
readr::write_csv(trimmed_counts_smoothed, path = output_file)

k <- filter(sample_design, sample_accession%in%colnames(trimmed_counts_smoothed), tissue%in%c('Retina','RPE','Cornea'))#%>%table(.[,'tissue'])
table(k$tissue)
