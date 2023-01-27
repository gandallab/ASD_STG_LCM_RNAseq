rm(list = ls());gc()
options(stringsAsFactors = F)

library(tidyverse)

inputs=list.files("./working_data/star_rsem/cell/multiqc_data/",glob2rx("multiqc_*.txt"),full.names = T)


for (input in inputs){
  dat=read.delim(input)
  print(c(input,nrow(dat)))
  dat$Sample=sapply(dat$Sample,function(x){str_split_fixed(x,"\\.",Inf)[1,1]})
  if (exists("output")){
    output=merge(output,dat,by="Sample")
  }else{
    output=dat
  }
}
rownames(output)=output$Sample;output$Sample=NULL
datMeta=read.delim("raw_data/metaData/metadata.lcm.txt",row.names = 1)
datMeta=datMeta[datMeta$Type == "Oligo",]
output=output[rownames(datMeta),]
output=as.matrix(output);mode(output)="numeric"
output1=Filter(function(x){(length(unique(x)) > 1) & (sum(is.na(x)) == 0)},as.data.frame(output))
output2=as.matrix(output1)
library(caret)
coline=cor(output2,method = "pe")
to_remove=findCorrelation(coline,cutoff = 0.98,names = T)
output3=output2[,!colnames(output2) %in% to_remove]

pca=prcomp(scale(output3,center = T,scale = T))
pcadata=pca$x
colnames(pcadata)=paste("seq",colnames(pcadata),sep = "")

library(SummarizedExperiment)


table(rownames(datMeta) %in% rownames(pcadata))
datMeta=merge(datMeta,pcadata,by="row.names")
datMeta$Sample=datMeta$Row.names;datMeta$Row.names=NULL

counts=read.delim("working_data/star_rsem/cell/genes.results.txt",row.names = 1,check.names = F)

counts=counts[,datMeta$Sample]
counts1=counts[rowSums(counts>0) > 0.8*ncol(counts),]
table(datMeta$Sample == colnames(counts1))
counts2=round(as.matrix(counts1))

se_oligo=SummarizedExperiment(assays = list(counts=counts2),colData=datMeta)
save(file = "working_data/summarizedExperiment/newrsem/se_oligo.RData",se_oligo)
