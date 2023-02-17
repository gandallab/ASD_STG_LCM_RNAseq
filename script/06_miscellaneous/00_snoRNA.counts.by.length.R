rm(list=ls());gc()
options(stringsAsFactors = F)

library(DESeq2)
library(ggplot2)
library(SummarizedExperiment)


load("working_data/summarizedExperiment/se_blocks_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/se_neuron_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/se_oligo_CPM_outlierRemoved.RData")

se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
se_oligo$Diagnosis = factor(se_oligo$Diagnosis, levels= c("Control", "Autism"))
se_neuron$Diagnosis = factor(se_neuron$Diagnosis, levels= c("Control", "Autism"))


# blocks
datMeta=colData(se_blocks)
countData = assays(se_blocks)$log2CPM
annot=as.data.frame(rowData(se_blocks))
table(rownames(annot) == rownames(countData))
snor=as.data.frame(countData[annot$gene_type == "snoRNA",])
snor$geneid=rownames(snor)
snor1 = snor %>% 
  gather("sample","log2CPM",1:55)
annot.snor=annot[annot$gene_type == "snoRNA",]
snor1$geneid=factor(snor1$geneid,levels = rownames(annot.snor)[order(annot.snor$width,decreasing = F)])
snor1$Diagnosis=datMeta$Diagnosis[match(snor1$sample,datMeta$Sample)]
snor2=snor1 %>%
  group_by(geneid) %>%
  summarise(medianlog2CPM=median(log2CPM)) %>%
  ungroup()
ggplot(snor1)+
  geom_point(aes(geneid,log2CPM,col=Diagnosis),alpha=0.5)+
  geom_line(data = snor2,aes(geneid,medianlog2CPM,group=1))+
  xlab("order by length")+
  theme(axis.text.x = element_text(angle = 90))

# neuron
datMeta=colData(se_neuron)
countData = assays(se_neuron)$log2CPM
annot=as.data.frame(rowData(se_neuron))
table(rownames(annot) == rownames(countData))
snor=as.data.frame(countData[annot$gene_type == "snoRNA",])
snor$geneid=rownames(snor)
snor1 = snor %>% 
  gather("sample","log2CPM",1:39)
annot.snor=annot[annot$gene_type == "snoRNA",]
snor1$geneid=factor(snor1$geneid,levels = rownames(annot.snor)[order(annot.snor$width,decreasing = F)])
snor1$Diagnosis=datMeta$Diagnosis[match(snor1$sample,datMeta$Sample)]
snor2=snor1 %>%
  group_by(geneid) %>%
  summarise(medianlog2CPM=median(log2CPM)) %>%
  ungroup()
ggplot(snor1)+
  geom_point(aes(geneid,log2CPM,col=Diagnosis),alpha=0.5)+
  geom_line(data = snor2,aes(geneid,medianlog2CPM,group=1))+
  xlab("order by length")+
  theme(axis.text.x = element_text(angle = 90))

# oligo
datMeta=colData(se_oligo)
countData = assays(se_oligo)$log2CPM
annot=as.data.frame(rowData(se_oligo))
table(rownames(annot) == rownames(countData))
snor=as.data.frame(countData[annot$gene_type == "snoRNA",])
snor$geneid=rownames(snor)
snor1 = snor %>% 
  gather("sample","log2CPM",1:18)
annot.snor=annot[annot$gene_type == "snoRNA",]
snor1$geneid=factor(snor1$geneid,levels = rownames(annot.snor)[order(annot.snor$width,decreasing = F)])
snor1$Diagnosis=datMeta$Diagnosis[match(snor1$sample,datMeta$Sample)]
snor2=snor1 %>%
  group_by(geneid) %>%
  summarise(medianlog2CPM=median(log2CPM)) %>%
  ungroup()
ggplot(snor1)+
  geom_point(aes(geneid,log2CPM,col=Diagnosis),alpha=0.5)+
  geom_line(data = snor2,aes(geneid,medianlog2CPM,group=1))+
  xlab("order by length")+
  theme(axis.text.x = element_text(angle = 90))

