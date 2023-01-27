rm(list = ls());gc()
options(stringsAsFactors = F)
library(tidyverse)

load("working_data/dge/dge.neuron.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")
dat$genenames=annot$genename[match(rownames(dat),annot$geneid)]
dat=dat[dat$padj < 0.05,]

gsea=read.delim("working_data/final/tables/dge.neuron.gsea.txt")

for (i in 1:nrow(gsea)){
  print(c(gsea[i,2],intersect(str_split_fixed(gsea[i,13],",",Inf)[1,],dat$genenames)))
}
