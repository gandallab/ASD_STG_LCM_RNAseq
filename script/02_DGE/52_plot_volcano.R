rm(list=ls());gc()
options(stringsAsFactors = F)

library(ggrepel)
library(scales)
library(tidyverse)

asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")


# blocks
load("working_data/dge/dge.blocks.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
dat$genename=annot$genename[match(rownames(dat),annot$geneid)]

fdrdat=dat[dat$padj < 0.05,]
pdat=dat[dat$padj > 0.05 & dat$pvalue < 0.05,]
restdat=dat[dat$pvalue > 0.05,]
labeldat=dat[dat$genename %in% c("HSPA1A","HSPA1B","DNAJB1","DNAJB4","HSPB1","HSPB8","BAG3","PTGES3","KCNH3","KCNIP1"),]


ggplot()+
  geom_point(data=fdrdat[fdrdat$log2FoldChange > 0,],aes(log2FoldChange,-log10(pvalue)),color="gold",size=2)+
  geom_point(data=fdrdat[fdrdat$log2FoldChange < 0,],aes(log2FoldChange,-log10(pvalue)),color="royalblue",size=2)+
  geom_point(data=pdat,aes(log2FoldChange,-log10(pvalue)),color="gray",size=1.5)+
  geom_point(data = restdat,aes(log2FoldChange,-log10(pvalue)),color="gray",size=1.5)+
  geom_text_repel(data = fdrdat[(!fdrdat$genename %in% labeldat$genename) & (fdrdat$padj < 0.02),],aes(log2FoldChange,-log10(pvalue),label=genename),segment.size = 0.1,col="black",fontface="bold",size=3)+
  geom_text_repel(data = labeldat,aes(log2FoldChange,-log10(pvalue),label=genename),segment.size = 0.1,col="red",fontface="bold",size=3.5)+
  scale_x_continuous(trans = "asinh")+
  theme_classic()+
  theme(text = element_text(family = "Arial",face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=16))


# neuron
rm(list=ls());gc()
options(stringsAsFactors = F)

library(ggrepel)
library(scales)
library(tidyverse)

asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")

load("working_data/dge/dge.neuron.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
dat$genename=annot$genename[match(rownames(dat),annot$geneid)]
dat$genetype=annot$genetype[match(rownames(dat),annot$geneid)]

dat.sno=dat[dat$genetype == "snoRNA",]
dat.sno=dat.sno[!is.na(dat.sno$genetype),]
dat.sno$genename[dat.sno$pvalue > 0.05]=""
dat.nonsno=dat[dat$genetype != "snoRNA",]

fdrdat=dat[dat$padj < 0.05,]
pdat=dat[dat$padj > 0.05 & dat$pvalue < 0.05,]
restdat=dat[dat$pvalue > 0.05,]
labeldat=dat[dat$genename %in% c("GAD1","GAD2","JUN","JUNB","FOS","SOX9","S1PR1","PPP1R16B","KCNJ2","NFKBID"),]

ggplot()+
  geom_point(data=dat.nonsno,aes(log2FoldChange,-log10(pvalue)),color="gray",size=1.5)+
  geom_point(data=dat.sno,aes(log2FoldChange,-log10(pvalue)),color="red",size=2)+
  geom_hline(linetype=2,size = 1,col="black",yintercept = c(1.3,-log10(max(dat$pvalue[dat$padj < 0.05]))))+
  geom_text_repel(data = dat.sno,aes(log2FoldChange,-log10(pvalue),label=genename),min.segment.length = unit(0, "lines"),col="red",fontface="bold")+
  scale_x_continuous(trans = "asinh")+
  theme_classic()+
  theme(text = element_text(family = "Arial",face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=16))

ggplot()+
  geom_point(data=fdrdat[fdrdat$log2FoldChange > 0,],aes(log2FoldChange,-log10(pvalue)),color="gold",size=2.5)+
  geom_point(data=fdrdat[fdrdat$log2FoldChange < 0,],aes(log2FoldChange,-log10(pvalue)),color="royalblue",size=2.5)+
  geom_point(data=pdat,aes(log2FoldChange,-log10(pvalue)),color="gray",size=2)+
  geom_point(data = restdat,aes(log2FoldChange,-log10(pvalue)),color="gray",size=1.5)+
  geom_text_repel(data = fdrdat[(!fdrdat$genename %in% labeldat$genename) & (fdrdat$padj < 0.03),],aes(log2FoldChange,-log10(pvalue),label=genename),segment.size = 0.1,col="black",size=3,fontface="bold")+
  geom_text_repel(data = labeldat,aes(log2FoldChange,-log10(pvalue),label=genename),segment.size = 0.1,col="red",fontface="bold",size=3.5)+
  scale_x_continuous(trans = "asinh")+
  theme_classic()+
  theme(text = element_text(family = "Arial",face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=16))

