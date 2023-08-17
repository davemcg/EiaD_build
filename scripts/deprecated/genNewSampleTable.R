setwd('NIH/autoRNAseq/')
library(dplyr)
load('comparisions/davidstart.Rdata')
load('/Volumes/McGaughey_S/Human_eyeIntegration_paper/data/core_tight.Rdata')
st <- filter(core_tight,sample_accession%in%tissues$sample_accession)
curs <- read.table('samplesrun_0928.tissue',stringsAsFactors = F,header=F,sep='\t')
to_Add <- filter(st,!sample_accession%in%curs$V1)
writeLines(to_Add$sample_accession,'sampleidtoadd',sep = '\n')
#run cmd on term

# while read p; do sqlite3 ref/SRAmetadb_072418.sqlite "SELECT library_layout FROM sra WHERE sample_accession='${p}'"; done<sampleidtoadd > paired
# while read p; do sqlite3 ref/SRAmetadb_072418.sqlite "SELECT run_accession FROM sra WHERE sample_accession='${p}'"; done<sampleidtoadd > run_ac


of <- read.table('outfile',sep = '\t',stringsAsFactors = F)
of <- filter(of,V1%in%to_Add$sample_accession)
p1 <- data.frame(of, tissue=to_Add$Tissue,subtissue=to_Add$Sub_Tissue,origin=to_Add$Origin,stringsAsFactors = F)
p1[grep('PAIRED',p1$V3),'V3'] <-'y'
p1[grep('SINGLE',p1$V3),'V3'] <-'n'
p1[grep('E-MTAB',p1$V1),'V2'] <- '.'   
p1[grep('E-MTAB',p1$V1),'V3'] <- 'y'   
p1$tissue <-  gsub(' ','.',p1$tissue)
p1$subtissue <- gsub(' ','.',p1$subtissue)
p1$origin <- gsub(' ','.',p1$origin)
colnames(p1) <- colnames(curs)
full <- rbind(curs,p1)
i <- 'SRR4425263
SRR4425267
SRR4425263
SRR4425262
SRR4425267
SRR4425266
SRR4425262
SRR4425264
SRR4425269
SRR4425266
SRR4425269
SRR4425268
SRR4425268
SRR4425265
SRR4425265
SRR4425264'
single_ended <- strsplit(i,'\n')%>%unlist
for( i in single_ended){
    full[grep(i,full$V2),2:3] <- c(i,'n')
}

write.table(full,'sampleTable1015_tissues.tab',row.names = F,col.names = F,quote = F,sep = '\t')
full[!full$V4%in%c('Retina','RPE','ESC','Cornea'),c('V4','V5')] <- 'body'
write.table(full,'sampleTable1015.tab',row.names = F,col.names = F,quote = F,sep = '\t')





