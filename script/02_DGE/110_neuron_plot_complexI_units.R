rm(list=ls());gc()
options(stringsAsFactors = F)

library(ggplot2)
library(ggrepel)
library(scales)
library(tidyverse)


annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")


# blocks
load("working_data/dge/dge.blocks.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
dat$genename=annot$genename[match(rownames(dat),annot$geneid)]
dat$genetype=annot$genetype[match(rownames(dat),annot$geneid)]
block=dat

# neuron
load("working_data/dge/dge.neuron.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
dat$genename=annot$genename[match(rownames(dat),annot$geneid)]
dat$genetype=annot$genetype[match(rownames(dat),annot$geneid)]
neuron=dat

block_nduf=block[grep("^NDUF",block$genename),c("genename","log2FoldChange")]
neuron_nduf=neuron[grep("^NDUF",neuron$genename),c("genename","log2FoldChange")]
colnames(block_nduf)=c("genename","block_fc")
colnames(neuron_nduf)=c("genename","neuron_fc")
dat=merge(block_nduf,neuron_nduf,by="genename")
ggplot(dat)+
  geom_point(aes(block_fc,neuron_fc))+
  xlim(-1.1,0.6)+
  ylim(-1.1,0.6)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept=0)+
  geom_abline(slope = 1,intercept = 0)+
  theme_bw()

dat=dat[order(dat$neuron_fc),]

dat1=dat %>%
  gather("type","fc",2:3) %>%
  arrange(fc)

dat1$genename=factor(dat1$genename,levels = dat$genename)
ggplot(dat1)+
  geom_bar(aes(genename,fc,fill=type),stat = "identity",position = "dodge")+
  xlab("")+
  ylab("log2 Fold Change")+
  theme_bw()+
  theme(text=element_text(family = "Arial",size = 16,face = "bold"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
