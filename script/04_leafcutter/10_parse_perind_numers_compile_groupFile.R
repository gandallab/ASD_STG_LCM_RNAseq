rm(list=ls());gc()
options(stringsAsFactors = F)

library(SummarizedExperiment)
library(tidyverse)

# blocks

load("se_blocks.RData")
datMeta=as.data.frame(colData(se_blocks))

datMeta1=datMeta[,c("Diagnosis" , "Sex" , "Age" , "RIN" , "X260.280" , "seqPC1" , "seqPC2" , "seqPC3")]
datMeta2=datMeta1[order(datMeta1$Diagnosis, decreasing = F),]

datSp=read.delim("block.cell.combined_perind_numers.counts",sep = " ",row.names = 1)
#rs=rowSums(datSp)
#table(rs == 0)
#tmp=datSp[1:10,]
table(rownames(datMeta2) %in% colnames(datSp))
datSp1=datSp[,rownames(datMeta2)]
datSp2=data.frame(clusterid=sapply(rownames(datSp1),function(x){str_split_fixed(x,":",Inf)[1,4]}),
                  counts=rowSums(datSp1)) 
datSp3 = datSp2 %>%
  group_by(clusterid) %>%
  summarise(totalcounts=sum(counts))

ind = (datSp2$clusterid %in% datSp3$clusterid[datSp3$totalcounts > 0])

datSp4=datSp1[ind,]
write.table(datMeta2,file = "blocks.group.txt",quote = F,sep = "\t",row.names = T)
write.table(datSp4,file = "blocks.perind_numers.counts",quote = F,sep = " ",row.names = T)

# neuron
rm(list=ls());gc()

load("se_neuron.RData")
datMeta=as.data.frame(colData(se_neuron))

datMeta1=datMeta[,c("Diagnosis" , "Sex" , "Age" , "Type_RNAseqRunNumber" , "seqPC1" , "seqPC2" , "seqPC3")]
datMeta2=datMeta1[order(datMeta1$Diagnosis, decreasing = T),]

datSp=read.delim("block.cell.combined_perind_numers.counts",sep = " ",row.names = 1)
#tmp=datSp1[1:10,]
table(make.names(rownames(datMeta2)) %in% colnames(datSp))
rownames(datMeta2)=make.names(rownames(datMeta2))
datSp1=datSp[,rownames(datMeta2)]
datSp2=data.frame(clusterid=sapply(rownames(datSp1),function(x){str_split_fixed(x,":",Inf)[1,4]}),
                  counts=rowSums(datSp1)) 
datSp3 = datSp2 %>%
  group_by(clusterid) %>%
  summarise(totalcounts=sum(counts))

ind = (datSp2$clusterid %in% datSp3$clusterid[datSp3$totalcounts > 0])

datSp4=datSp1[ind,]

write.table(datMeta2,file = "neuron.group.txt",quote = F,sep = "\t",row.names = T)
write.table(datSp1,file = "neuron_unfilter.perind_numers.counts",quote = F,sep = " ",row.names = T)


