library(data.table)
library(dplyr)

args <- commandArgs(TRUE)
file <- args[1]

data<-fread(file)

TEfamily <- fread("../../TE_annotation/TE_library_family.csv")

dataFam<-left_join(data,TEfamily[,c("Consensus_ID","Family")], by = c("V5" = "Consensus_ID"))

write.table(dataFam, paste0(file,".family"),quote = F, col.names = F, sep = "\t",row.names = F)