setwd('/data/swamyvs/autoRNAseq')
#setwd('~/NIH/autoRNAseq')
library(tximport)
library(dplyr)
library(qsmooth)
library(Rtsne)
library(ggplot2)
library(dbscan)
library(edgeR)
library(rtracklayer)
args=commandArgs(trailingOnly = T)#for gtf
e_Dist <- function(p1,p2) return(sqrt(sum((p1-p2)^2)))

###write.
# Qsmooth > remove median counts > remove lowly expressed genes > tSNE > DBSCAN
sample_design <- read.table(args[1],stringsAsFactors = F,header=F, sep = '\t')
gtf <- readGFF(args[2])%>%dplyr::filter(type=='transcript')
anno <- gtf[,c("gene_id", "gene_name", "transcript_id")]
colnames(sample_design) <- c('sample_accession', 'run_accession', 'paired','tissue','sub-tissue','origin')
files0 <-paste0('RE_quant_files/',sample_design$sample_accession, '/quant.sf')

# remove samples which failed to align
files0 <- files0[file.exists(files0)]
samplenames <- strsplit(files0,'/' )%>% sapply(function(x) x[2])
sample_design <-  filter(sample_design, sample_accession%in%samplenames)
#load('tpms.Rdata')
#txi.counts <- tximport(files=files0,tx2gene =  anno[,3:2],type = "salmon")
txi.lsTPMs <- tximport(files=files0,tx2gene =  anno[,3:2],type = "salmon", countsFromAbundance = "lengthScaledTPM")
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
keep_genes <- which(rowSums(tpms)>=ncol(tpms)/2)# revove gene w less than an average count of .5
tpms <- tpms[keep_genes,]
sample_medians <- apply(tpms,2,function(x) median(x))

# more genes => lower medians, looks like 3 givees results most similar to david, see comparingQC.r
keep_median <- which(sample_medians>50)# criteria for removing samples with low counts
# removed <- samplenames[which(sample_medians<=50)]
# intersect(d.removed,removed)%>%length
# View(sample_design[removed,])
# View(sample_design[setdiff(removed,d.removed),])


tpms <- tpms[,keep_median]

#normalize for liberry size
norm <- DGEList(tpms)
norm <- calcNormFactors(norm)
norm_counts <- norm$counts
#extract scaling factor for each sample and  multiply
correction <- norm$samples %>% data.frame() %>% .[['norm.factors']]
lsTPM_librarySize <- norm_counts %*% diag(correction)
colnames(lsTPM_librarySize) <- colnames(tpms)

#quantile normalize samples
sample_design <- sample_design[sample_design$sample_accession%in%colnames(lsTPM_librarySize),]
qs <- qsmooth(object = lsTPM_librarySize,groupFactor = as.factor(sample_design$tissue))
lstpms_smoothed <- as.data.frame(qsmoothData(qs))

colnames(lstpms_smoothed) <- colnames(lsTPM_librarySize)
################
# load('/volumes/McGaughey_S/Human_eyeIntegration_paper/data/lengthScaledTPM_processed_2017_02.Rdata')
# david_lsTPMSqsft <- as.data.frame(lengthScaledTPM_processed)
# g <- intersect(rownames(david_lsTPMSqsft),rownames(lstpms_smoothed))
# s <- intersect(colnames(david_lsTPMSqsft),colnames(lstpms_smoothed))
# # 18701 genes in common, 849 samples
# cor(cbind(lstpms_smoothed[g,s[1]],david_lsTPMSqsft[g,s[1]]),method = 'spearman')
# k <- vector(mode = 'numeric',length = length(s))
# names(k) <- s
# for(i in s)k[i] <- cor(cbind(lstpms_smoothed[g,i],david_lsTPMSqsft[g,i]),method = 'spearman')[2,1]








#cluster with tSNE, then run clustered data throught
tpms_smoothed_filtered <- lstpms_smoothed
set.seed(23235)
tsne_out <- Rtsne(as.matrix(log2(t(tpms_smoothed_filtered)+1)),perplexity = 40, check_duplicates = FALSE, theta=0.0 )
tsne_plot <- data.frame(tsne_out$Y,sample_design[sample_design$sample_accession%in%colnames(tpms_smoothed_filtered),])
# ggplot(tsne_plot,aes(x=X1,y=X2, col=tissue))+
#   geom_point(size=2)+
#   theme_minimal()



#find outliers based on tsne grouping. calculate the center of each group, and look for points n standard deviations away
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
write.csv(trimmed_counts_smoothed,'results/smoothed_filtered_tpms.csv')


