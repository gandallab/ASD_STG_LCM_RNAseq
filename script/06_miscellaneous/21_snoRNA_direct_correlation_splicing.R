library(tidyverse)

rm(list = ls());gc()
options(stringsAsFactors = F)

dat=read.delim("working_data/snorna_directCor.txt",header = F)
dat$fdr=p.adjust(dat$V3)
dat$cluster=sapply(dat$V1,function(x){
  a=str_split_fixed(x,":",Inf)[1,]
  paste("chr",a[1],":",a[length(a)],sep = "")
})


sg=read.delim("working_data/leafcutter/v2/neuron_cluster_significance.txt")

table(dat$cluster %in% sg$cluster)

sg=sg[sg$cluster %in% dat$cluster,]

table(sg1$status)

dat$gene=sg$genes[match(dat$cluster,sg$cluster)]
write.table(dat,file = "working_data/snoRNA_directCorSplicing_annotated.txt",quote = F,sep = "\t",row.names = F)
dat1=dat[dat$fdr < 0.1,]
sg1=sg[sg$cluster %in% dat1$cluster,]
