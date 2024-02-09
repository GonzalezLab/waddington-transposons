
library(data.table)
library(dplyr)

args <- commandArgs(TRUE)
file <- args[1]

data<-fread(file)

TEfamily <- fread("../../TE_annotation/TE_library_family.csv")

data<-tidyr::separate(data, col = `Transposon:Start:End:Strand`, sep = ":", into="Consensus_ID")


dataFamtmp<-left_join(data,TEfamily[,c("Consensus_ID","Family")], by = "Consensus_ID")
data<-fread(file)
dataFam <- cbind(data,dataFamtmp$Family)
write.table(dataFam, paste0(file,".family"),quote = F, col.names = T, sep = "\t",row.names = F)
