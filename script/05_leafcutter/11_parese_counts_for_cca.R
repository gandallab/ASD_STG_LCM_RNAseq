# neuron
rm(list=ls());gc()
options(stringsAsFactors = F)

library(SummarizedExperiment)
library(tidyverse)

load("working_data/summarizedExperiment/fc_2pass/se_neuron_CPM_outlierRemoved.RData")
datMeta=as.data.frame(colData(se_neuron))

datMeta1=datMeta[,c("Diagnosis" , "Sex" , "Age" , "Type_RNAseqRunNumber" , "seqPC1" , "seqPC2" , "seqPC3")]
datMeta2=datMeta1[order(datMeta1$Diagnosis, decreasing = T),]

datSp=read.delim("working_data/leafcutter/v2/cca/block.cell.combined_perind.counts.gz",sep = " ",row.names = 1)
#tmp=datSp1[1:10,]
table(make.names(rownames(datMeta2)) %in% colnames(datSp))
rownames(datMeta2)=make.names(rownames(datMeta2))
datSp1=datSp[,rownames(datMeta2)]
ind=sapply(rownames(datSp1),function(x){
  str_split_fixed(x,":",Inf)[1,1]
})
ind1=(ind %in% paste("chr",seq(1,22),sep = ""))
datSp2=datSp1[ind1,]

datSp3=rownames_to_column(datSp2,"chrom")
#write.table(datMeta2,file = "working_data/leafcutter/v2/neuron.group.txt",quote = F,sep = "\t",row.names = T)
write.table(datSp3,file = "working_data/leafcutter/v2/cca/neuron_unfilter.perind.counts",quote = F,sep = " ",row.names = F)
