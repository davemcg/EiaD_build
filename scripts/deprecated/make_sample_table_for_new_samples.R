setwd('~/NIH/autoRNAseq')
library(RSQLite)
library(SRAdb)
library(tidyverse)
library(stringr)
library(dplyr)
#getSRAdbFile(destdir='ref/',destfile='SRAmetadb.11-06-18.sqlite.gz') # do periodically. latest from 2017-01-19, 2017-04-17
sqlfile <- '~/NIH/autoRNAseq/ref/SRAmetadb.11-06-18.sqlite' 
sra_con <- dbConnect(RSQLite::SQLite(),sqlfile)
new_samples <- read.table('samples_geo',stringsAsFactors = F,sep='\n')[,1]%>%as.character()
res <- list()
for(i in 1:nrow(new_samples)){
    targ <- new_samples[i,]
    samps <- strsplit(targ$Relevant.Sample.Accessions,',')%>%unlist
    n <- length(samps)
    df <- data.frame(sample=samps, run=NA, paired=NA,tissue=rep(targ$tissue,n),subtissue=rep(targ$Sub.Tissue,n),origin=NA, stringsAsFactors = F)
    res[[length(res)+1]] <- df
}
tab <- do.call(rbind,res)
tab <- tab[!is.na(tab$sample),]
tab$sample <- gsub(' ','',tab$sample)
tab <- filter(tab, !sample%in%c('SRS1747723','SRS1588794','SRS1588796','SRS1747748','SRS1747747','SRS1747749'))# aleady in sample
#tab <- filter(tab, !grepl('SRX',sample)) # not in db
que <- "SELECT study_accession FROM sra WHERE sample_accession='BLANK'"
info <- lapply(tab$V1,function(x) dbGetQuery(sra_con, gsub('BLANK',x,que)))
to_dl <- do.call(rbind, info)
write(to_dl, 'new_runs_tdl',sep = 'n')
info_runs <- lapply(info,function(x) paste(x[,1], collapse = ','))
info_paired <- lapply(info, function(x) ifelse(any(grepl('PAIRED',x[,2])),'y','n'))
tab$run <- info_runs%>%unlist()
tab$paired <- info_paired%>%unlist
tab$origin <- 'new'
old_tab <- read.table('sampleTable1015_tissues.tab',stringsAsFactors = F, header = F, sep = '\t')
old_tab$V6 <- 'old'
colnames(old_tab) <- colnames(tab)
full_tab <- rbind(old_tab,tab)
full_tab <- full_tab[!duplicated(full_tab[,1]),]
write.table(full_tab,'sampleTable1106_tissues.tab',row.names = F,col.names = F,sep = '\t',quote = F)
b <- !full_tab$tissue%in%c('Retina','Cornea','RPE','ESC')
full_tab$tissue[b] <- 'body'
full_tab$subtissue[b] <- 'body'
write.table(full_tab,'sampleTable1106.tab',row.names = F,col.names = F,sep = '\t',quote = F)
#######

new <- read.csv('samples_geo.csv', header = F, stringsAsFactors = F)
old <- read.table('sampleTable1106_tissues.tab', stringsAsFactors = F)
comb <- rbind(old,new, stringsAsFactors=F)
write.table(comb, 'sampleTable1127_tissues.tab', sep = '\t', row.names = F, col.names = F, quote = F)
b <- !comb$V4%in%c('Retina','Cornea','RPE','ESC')
comb$V4[b] <- 'body'
comb$V5[b] <- 'body'
write.table(comb, 'sampleTable1127.tab', sep = '\t', row.names = F, col.names = F, quote = F)
