notes <- '
-for a single subtisue, there are different sites present in diffferent comparisons, might want to look into that


'
library(dplyr)



i_files <- dir('rmats_out')

subtissues_PE <- c("Retina_Adult.Tissue", "RPE_Cell.Line", "ESC_Stem.Cell.Line" , "RPE_Adult.Tissue" )# add body back in  at some point
k <- combn(subtissues_PE,2,simplify = F)
i_event <- 'SE.MATS.JC.txt'


 
combine_PE_SE <- function(combination,event,files){
    target_files <- files[grepl(combination[1],files)]%>%.[grepl(combination[2],.)]
    if(length(target_files)==1){
    tmp <- paste('rmats_out',target_files,event, sep = '/')%>%read.table(header = T,sep = '\t',stringsAsFactors = F)
    path <- paste0('rmats_comb/',combination[1],'_VS_',combination[2])
    dir.create(path = path)
    write.table(tmp,paste(path,event,sep='/'),row.names = F,col.names = T, quote = F,sep = '\t')
    return(0)
    }else if(length(target_files)==0){
      print('REEEEEEEEEEEEEE')
      return(1)
    }
    names(target_files) <- grepl('_PE',target_files)%>%ifelse('PE','SE')
    samp_PE <- paste('rmats_out',target_files[grep('_PE',target_files)],event,sep = '/')%>%read.table(header = T,sep = '\t',stringsAsFactors = F)
    samp_SE <- paste('rmats_out',target_files[grep('_SE',target_files)],event,sep = '/')%>%read.table(,header = T,sep = '\t',stringsAsFactors = F)
    samp_PE[,13:16] <- apply(samp_PE[,13:14],2, function(x) sapply(x,function(y) strsplit(y,',')%>%unlist%>%as.numeric%>%sum) )
    samp_SE[,13:16] <- apply(samp_SE[,13:14],2, function(x) sapply(x,function(y) strsplit(y,',')%>%unlist%>%as.numeric%>%sum) )
    
    #the first tissue is the first one in combination, the second tissue is the second 
    st1_se <- paste0(c('IJC_SAMPLE_','SJC_SAMPLE_'),  grep(combination[1],strsplit(target_files['SE'],'VS')%>%unlist))
    st2_se <- paste0(c('IJC_SAMPLE_','SJC_SAMPLE_'),  grep(combination[2],strsplit(target_files['SE'],'VS')%>%unlist))
    st1_pe <- paste0(c('IJC_SAMPLE_','SJC_SAMPLE_'),  grep(combination[1],strsplit(target_files['PE'],'VS')%>%unlist))
    st2_pe <- paste0(c('IJC_SAMPLE_','SJC_SAMPLE_'),  grep(combination[2],strsplit(target_files['PE'],'VS')%>%unlist))
    #test <- full_join(samp_SE,samp_PE, by=c("chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"))
    good_cols <- c("GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE",'IJC_SAMPLE_1','SJC_SAMPLE_1','IJC_SAMPLE_2','SJC_SAMPLE_2',"PValue","FDR")
    samp_PE <- samp_PE[,c("GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE",st1_pe,st2_pe,"PValue","FDR")]
    colnames(samp_PE) <- good_cols
    samp_SE <- samp_SE[,c("GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE",st1_se,st2_se,"PValue","FDR")]
    colnames(samp_SE) <- good_cols
      
    z_merge <- full_join(samp_PE,samp_SE,by=c("chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"))
    # fill in na values for info
    z_merge[is.na(z_merge$IJC_SAMPLE_1.x),c(1:2,11:14)]  <-z_merge[is.na(z_merge$IJC_SAMPLE_1.x),17:22] 
    # fill in na p-values by just replicatng p-value from sample with valid values, so when we average, it will stay the same
    z_merge$PValue.x[is.na(z_merge$PValue.x)] <- z_merge$PValue.y[is.na(z_merge$PValue.x)]
    z_merge$PValue.y[is.na(z_merge$PValue.y)] <- z_merge$PValue.x[is.na(z_merge$PValue.y)]
    new_pvalue <- rowMeans(z_merge[c('PValue.x','PValue.y')])
    new_fdr <- p.adjust(new_pvalue,method = "BH")
    final <- data.frame(z_merge[,1:14],new_pvalue,new_fdr,stringsAsFactors = F)
    colnames(final) <- good_cols
    path <- paste0('rmats_comb/',combination[1],'_VS_',combination[2])
    dir.create(path = path)
    write.table(final,paste(path,event,sep='/'),row.names = F,col.names = T, quote = F,sep = '\t')
}

for(i in 1:length(k)) {
  i_combination <- k[[i]]
  combine_PE_SE(combination = i_combination,event = i_event,files = i_files )
}

# files <- dir('rmats_out')
# subtissue <- 'Retina_Adult.Tissue'
# files.st <- files[grep(subtissue,files)]
# event <- 'SE.MATS.JC.txt'
# first <- TRUE


##Combine all different comparisons for a specific tissue

combine_rmats_output <- function(files,subtissue,event,first=TRUE){
    files.st <- files[grep(subtissue,files)]
    for(comparison in files.st){
        #generate the first comparison
        if(first==TRUE){
              first <- FALSE
              test1 <- paste('rmats_out',comparison,event, sep = '/')%>% read.table(,header = T,sep = '\t',stringsAsFactors = F)
              st_counts <- c('IJC_SAMPLE_1','SJC_SAMPLE_1')
              comp <- paste0(c("PValue","FDR"),'.',comparison )
              if(strsplit(comparison,'VS')%>%unlist%>%grepl(subtissue,.)%>%.[2]) st_counts <- c('IJC_SAMPLE_2','SJC_SAMPLE_2')
              cols <- c( "GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE",st_counts,"PValue","FDR" )
              test1 <- test1[,cols]
              test1[,11:12] <- apply(test1[,11:12],2, function(x) sapply(x,function(y) strsplit(y,',')%>%unlist%>%as.numeric%>%sum) )
              colnames(test1)<- c( "GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE", 'IJC_SAMPLE_1','SJC_SAMPLE_1',comp)
          } else{# now add rest of comparisons to first
              test2 <- paste('rmats_out',comparison,event, sep = '/')%>% read.table(,header = T,sep = '\t',stringsAsFactors = F)
              st_counts <- c('IJC_SAMPLE_1','SJC_SAMPLE_1')
              # sample be first or second sample , so account for that
              if(strsplit(comparison,'VS')%>%unlist%>%grepl(subtissue,.)%>%.[2]) st_counts <- c('IJC_SAMPLE_2','SJC_SAMPLE_2')
              cols <- c( "GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE",st_counts,"PValue","FDR" )
              test2 <- test2[,cols]
              # counts are presented as a comma sep list, so split and sum for total count for a tissue
              test2[,11:12] <- apply(test2[,11:12],2, function(x) sapply(x,function(y) strsplit(y,',')%>%unlist%>%as.numeric%>%sum) )
              comp <- paste0(c("PValue","FDR"),'.',comparison)
              colnames(test2)<- c( "GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE", 'IJC_SAMPLE_1','SJC_SAMPLE_1',comp)
              #join old and new dfs together, then fill in any events only foun in new,  and the format
              test_join <- full_join(test1,test2, by=c("chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE"))
              end=ncol(test_join)
              test_join[is.na(test_join$IJC_SAMPLE_1.x),c(1:2,11:12)] <-test_join[is.na(test_join$IJC_SAMPLE_1.x),c((end-5):(end-4),(end-1):end)] 
              test_join <-test_join[,-c(end-2,end-3,end-4,end-5)]
              colnames(test_join) <- c(colnames(test1),comp)
              test1 <- test_join
        }
    
    }
}





