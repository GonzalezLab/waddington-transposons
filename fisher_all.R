suppressPackageStartupMessages({
library(data.table)
library(dplyr)
library(rstatix)
})

args <- commandArgs(TRUE)
program <- args[2]
class<-args[3]
type<-args[4]

if (program == "PopTE2") {
  
  fisher <- read.table(args[1],sep="\t",header = F,fill = NA)
  if ( ncol(fisher) > 9) {
  colnames(fisher) <- c("Line", "TE", "chr", "start","end", "Cond", "Analysis", "Frequency","ReadsPresence","ReadsTotal")
  fisher$ReadsAbsence <- fisher$ReadsTotal - fisher$ReadsPresence
  rownames(fisher) <- fisher$Cond
  fisher_subset <- fisher[fisher$Cond %in% c("A", "N", "C", "P0"), c("ReadsPresence", "ReadsAbsence")]
  
  if ( nrow(na.omit(fisher_subset)) > 1 ) {
  fisher_test(as.matrix(na.omit(fisher_subset)))
  pairwise_fisher_test(as.matrix(na.omit(fisher_subset)), detailed = F,p.adjust.method = "fdr")
  
  resultFisher<-pairwise_fisher_test(as.matrix(na.omit(fisher_subset)), detailed = F,p.adjust.method = "fdr")
  resultFisher<-left_join(resultFisher,fisher[,c("Cond","Frequency")], by = c("group1"="Cond"))
  resultFisher<-left_join(resultFisher,fisher[,c("Cond","Frequency")], by = c("group2"="Cond"))
  resultFisher<-dplyr::rename(resultFisher,"Frequency1" = "Frequency.x")
  resultFisher<-dplyr::rename(resultFisher,"Frequency2" = "Frequency.y")
  
  #resultFisher <- merge(resultFisher,fisher[,c("Cond","Frequency")])
  finalTab<-cbind(unique(fisher[,c("Line","TE","chr","start","end")]),fisher_test(as.matrix(na.omit(fisher_subset))),resultFisher)
  
  
  write.table(finalTab,file=paste0("fisher_PopTE2_all_",class,"_",type,".tab"),append = T,quote=F,row.names = F,col.names = F,sep="\t")
  
  if (!is.na(fisher[fisher$Cond=="C",]$Frequency) && !is.na(fisher[fisher$Cond=="P0",]$Frequency)) {
  if (fisher[fisher$Cond=="C",]$Frequency > fisher[fisher$Cond=="P0",]$Frequency ) {
  FC <- (fisher[fisher$Cond==type,]$Frequency - fisher[fisher$Cond=="C",]$Frequency)/(fisher[fisher$Cond=="C",]$Frequency - fisher[fisher$Cond=="P0",]$Frequency)
  
  resultFC <- cbind(reshape(fisher[,c("Line","TE","chr","start","end","Cond","Frequency")], idvar = c("Line", "TE", "chr", "start", "end"), timevar = "Cond", direction = "wide"),FC)

  write.table(resultFC,file=paste0("FC_PopTE2_all_",class,"_",type,".tab"),append = T,quote=F,row.names = F,col.names = F,sep="\t")
  } else {
    FC <- (fisher[fisher$Cond==type,]$Frequency - fisher[fisher$Cond=="C",]$Frequency)/(fisher[fisher$Cond=="P0",]$Frequency-fisher[fisher$Cond=="C",]$Frequency)
    
    resultFC <- cbind(reshape(fisher[,c("Line","TE","chr","start","end","Cond","Frequency")], idvar = c("Line", "TE", "chr", "start", "end"), timevar = "Cond", direction = "wide"),FC)
    
    write.table(resultFC,file=paste0("FC_PopTE2_all_",class,"_",type,".tab"),append = T,quote=F,row.names = F,col.names = F,sep="\t")
  }
  } else {
    FC <- "NA"
    
    resultFC <- cbind(reshape(fisher[,c("Line","TE","chr","start","end","Cond","Frequency")], idvar = c("Line", "TE", "chr", "start", "end"), timevar = "Cond", direction = "wide"),FC)
    
    write.table(resultFC,file=paste0("FC_PopTE2_all_",class,"_",type,".tab"),append = T,quote=F,row.names = F,col.names = F,sep="\t")
  }
  
  
  
  } else {
    df <- data.frame(
      n = "NA",
      p = "NA",
      p.signif = "NA",
      group1 = "NA",
      group2 = "NA",
      n = "NA",
      p = "NA",
      p.adj = "NA",
      p.adj.signif = "NA",check.names = FALSE
    )
    colnames(fisher) <- c("Line", "TE","chr","start","end", "Cond", "Analysis", "Frequency","ReadsPresence","ReadsTotal")
    
    finalTab<-cbind(unique(fisher[,c("Line","TE","chr","start","end")]),df,fisher[,c("Frequency","Cond")])
    write.table(finalTab,file=paste0("fisher_PopTE2_all_",class,"_",type,".tab"),append = T,quote=F,row.names = F,col.names = F,sep="\t")
  } 
  } else {
      df <- data.frame(
        n = "NA",
        p = "NA",
        p.signif = "NA",
        group1 = "NA",
        group2 = "NA",
        n = "NA",
        p = "NA",
        p.adj = "NA",
        p.adj.signif = "NA",check.names = FALSE
      )
      colnames(fisher) <- c("Line", "TE","chr","start","end", "Cond", "Analysis", "Frequency")
      
      finalTab<-cbind(unique(fisher[,c("Line","TE","chr","start","end")]),df,fisher[,c("Frequency","Cond")])
      
      write.table(finalTab,file=paste0("fisher_PopTE2_all_",class,"_",type,".tab"),append = T,quote=F,row.names = F,col.names = F,sep="\t")
  }
  
} else if  (program =="TEMP2") {
  
  
  fisher <- read.table(args[1],sep="\t",header = F,fill = NA)
  
  if ( ncol(fisher) > 9) {
  colnames(fisher) <- c("Line", "TE","chr","start","end", "Cond", "Analysis", "Frequency","ReadsPresence","ReadsAbsence")
  rownames(fisher) <- fisher$Cond
  fisher_subset <- fisher[fisher$Cond %in% c("A", "N", "C", "P0"), c("ReadsPresence", "ReadsAbsence")]
  if ( nrow(na.omit(fisher_subset)) > 1 ) {
  fisher_test(as.matrix(na.omit(fisher_subset)))
  resultFisher<-pairwise_fisher_test(as.matrix(na.omit(fisher_subset)), detailed = F,p.adjust.method = "fdr")
  resultFisher<-left_join(resultFisher,fisher[,c("Cond","Frequency")], by = c("group1"="Cond"))
  resultFisher<-left_join(resultFisher,fisher[,c("Cond","Frequency")], by = c("group2"="Cond"))
  resultFisher<-dplyr::rename(resultFisher,"Frequency1" = "Frequency.x")
  resultFisher<-dplyr::rename(resultFisher,"Frequency2" = "Frequency.y")
  
  #resultFisher <- merge(resultFisher,fisher[,c("Cond","Frequency")])
  finalTab<-cbind(unique(fisher[,c("Line","TE","chr","start","end")]),fisher_test(as.matrix(na.omit(fisher_subset))),resultFisher)
  write.table(finalTab,file=paste0("fisher_TEMP2_all_",class,"_",type,".tab"),append = T,quote=F,row.names = F,col.names = F,sep="\t")
  
  
  if (!is.na(fisher[fisher$Cond=="C",]$Frequency) && !is.na(fisher[fisher$Cond=="P0",]$Frequency)) {
    if (fisher[fisher$Cond=="C",]$Frequency > fisher[fisher$Cond=="P0",]$Frequency ) {
      FC <- (fisher[fisher$Cond==type,]$Frequency - fisher[fisher$Cond=="C",]$Frequency)/(fisher[fisher$Cond=="C",]$Frequency - fisher[fisher$Cond=="P0",]$Frequency)
      
      resultFC <- cbind(reshape(fisher[,c("Line","TE","chr","start","end","Cond","Frequency")], idvar = c("Line", "TE", "chr", "start", "end"), timevar = "Cond", direction = "wide"),FC)
      
      write.table(resultFC,file=paste0("FC_TEMP2_all_",class,"_",type,".tab"),append = T,quote=F,row.names = F,col.names = F,sep="\t")
    } else {
      FC <- (fisher[fisher$Cond==type,]$Frequency - fisher[fisher$Cond=="C",]$Frequency)/(fisher[fisher$Cond=="P0",]$Frequency-fisher[fisher$Cond=="C",]$Frequency)
      
      resultFC <- cbind(reshape(fisher[,c("Line","TE","chr","start","end","Cond","Frequency")], idvar = c("Line", "TE", "chr", "start", "end"), timevar = "Cond", direction = "wide"),FC)
      
      write.table(resultFC,file=paste0("FC_TEMP2_all_",class,"_",type,".tab"),append = T,quote=F,row.names = F,col.names = F,sep="\t")
    }
  } else {
    FC <- "NA"
    
    resultFC <- cbind(reshape(fisher[,c("Line","TE","chr","start","end","Cond","Frequency")], idvar = c("Line", "TE", "chr", "start", "end"), timevar = "Cond", direction = "wide"),FC)
    
    write.table(resultFC,file=paste0("FC_TEMP2_all_",class,"_",type,".tab"),append = T,quote=F,row.names = F,col.names = F,sep="\t")
  }
  
  
  
  }
  else {
      df <- data.frame(
        n = "NA",
        p = "NA",
        p.signif = "NA",
        group1 = "NA",
        group2 = "NA",
        n = "NA",
        p = "NA",
        p.adj = "NA",
        p.adj.signif = "NA",check.names = FALSE
      )
      colnames(fisher) <- c("Line", "TE", "chr","start","end", "Cond", "Analysis","Frequency","ReadsPresence","ReadsAbsence")
      
      finalTab<-cbind(unique(fisher[,c("Line","TE","chr","start","end")]),df,fisher[,c("Frequency","Cond")])
      
      write.table(finalTab,file=paste0("fisher_TEMP2_all_",class,"_",type,".tab"),append = T,quote=F,row.names = F,col.names = F,sep="\t")
      
  }
  
  } else {
    df <- data.frame(
      n = "NA",
      p = "NA",
      p.signif = "NA",
      group1 = "NA",
      group2 = "NA",
      n = "NA",
      p = "NA",
      p.adj = "NA",
      p.adj.signif = "NA",check.names = FALSE
    )
    colnames(fisher) <- c("Line", "TE", "chr","start","end", "Cond", "Analysis","Frequency")
    
    finalTab<-cbind(unique(fisher[,c("Line","TE","chr","start","end")]),df,fisher[,c("Frequency","Cond")])
    
    write.table(finalTab,file=paste0("fisher_TEMP2_all_",class,"_",type,".tab"),append = T,quote=F,row.names = F,col.names = F,sep="\t")
  }
  
}
