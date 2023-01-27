library(SummarizedExperiment)
library(tidyverse)
library(psych)

rm(list = ls());gc()
options(stringsAsFactors = F)

load("working_data/wgcna/voom.forWGCNA.input.neuron.fc.manualRemove.RData")
load("working_data/dge/dge.neuron.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")
snor=annot[annot$genetype == "snoRNA",]
dat1=dat[rownames(dat) %in% snor$geneid,]
dat2=dat1[dat1$log2FoldChange < -0.2 & dat1$pvalue < 0.5,]

rownames(datExpr)=sapply(rownames(datExpr),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
datExpr1=datExpr[rownames(dat2),]

# datExpr2=datExpr1
# rownames(datExpr2)=annot$genename[match(rownames(datExpr2),annot$geneid)] 
write.table(datExpr1,file = "working_data/snoRNA/snorna.individual.txt",quote = F,sep = "\t",row.names = T)

pca=prcomp(scale(t(datExpr1),center = T,scale = F))
pcadata=pca$x

pcadata1=merge(datMeta,pcadata,by="row.names")

cor(pcadata[,1],colMeans(datExpr1))

write.table(pcadata,file = "working_data/snorna.pca.txt",quote = F,sep = "\t",row.names = T)

tmp=corr.test(pcadata,t(datExpr1),adjust = "none")
tmp1=t(rbind(tmp$r,tmp$p))
tmp1[,3]=row.names(tmp1)
tmp1=cor(pcadata,t(datExpr1))
cor.test(pcadata[,1],datExpr1[1,])

