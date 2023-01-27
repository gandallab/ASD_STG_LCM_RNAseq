rm(list = ls());gc()
options(stringsAsFactors = F)

library(tidyverse)

inputs=list.files("./working_data/multiqc_2pass/blocks/",glob2rx("multiqc_*.txt"),full.names = T)


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
load("working_data/summarizedExperiment/se_blocks.RData")
se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
datMeta=as.data.frame(colData(se_blocks))
datMeta=datMeta[,-c(16:82)]
datMeta=merge(datMeta,pcadata,by="row.names")
datMeta=datMeta[!is.na(datMeta$Diagnosis),]
rownames(datMeta)=datMeta$Row.names;datMeta$Row.names=NULL

rm(output)
counts=list.files("C:/Users/Pan Zhang/Downloads/fc/pass2/countsfile_block/","*counts.txt")
for (count in counts){
  dat=read.delim(paste("C:/Users/Pan Zhang/Downloads/fc/pass2/countsfile_block/",count,sep = ""),skip = 1)
  dat=dat[,c(1,7,8)]
  colnames(dat)[3]=gsub(".counts.txt","",count)
  if (exists("output")){
    if (sum(output$Geneid != dat$Geneid) > 0){
      print(count)
      next()
    }
    output=cbind(output,dat[,3,drop=F])
  }
  else{
    output=dat
  }
}

rownames(output)=output$Geneid;output$Geneid=NULL;output$gene_name=NULL
output1=as.matrix(output)

table(datMeta$Sample %in% colnames(output1))
output2=output1[,datMeta$Sample]
output3=output2[rowSums(output2>0)>0.8*ncol(output2),]
se_blocks=SummarizedExperiment(assays = list(counts=output3),colData=datMeta)
save(file = "working_data/summarizedExperiment/fc_2pass//se_blocks.RData",se_blocks)
