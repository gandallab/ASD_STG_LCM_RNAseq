options(stringsAsFactors = F)
library(tidyverse)

# block
rm(list=ls());gc()

load("working_data/dge/dge.blocks.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]

dat=dat[dat$padj <= 0.05,]

dat$direct = "up"
dat$direct[dat$log2FoldChange < 0]="down"
dat$FC= 2 ^ dat$log2FoldChange

ggplot()+
  geom_histogram(data = dat[dat$direct == "up",],aes(x=FC,y=..count..),binwidth = 0.1,fill="gold",color="black")+
  geom_histogram(data = dat[dat$direct == "down",],aes(x= 1/FC,y=-..count..),binwidth = 0.1,fill="royalblue",color="black")+
  xlab("RNA-seq Fold Change")+
  ylab("Number of DE genes ( FDR < 0.05)")+
  scale_x_continuous(breaks = seq(1,4,0.5))+
  scale_y_continuous(breaks = seq(-10,25,5))+
  theme_classic()+
  theme(text = element_text(family = "Arial",face = "bold"),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18))


# neuron
rm(list=ls());gc()

load("working_data/dge/dge.neuron.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]

dat=dat[dat$padj <= 0.05,]

dat$direct = "up"
dat$direct[dat$log2FoldChange < 0]="down"
dat$FC= 2 ^ dat$log2FoldChange

ggplot()+
  geom_histogram(data = dat[dat$direct == "up",],aes(x=FC,y=..count..),binwidth = 0.3,fill="gold",color="black")+
  geom_histogram(data = dat[dat$direct == "down",],aes(x= 1/FC,y=-..count..),binwidth = 0.3,fill="royalblue",color="black")+
  xlab("RNA-seq Fold Change")+
  ylab("Number of DE genes ( FDR < 0.05)")+
  scale_x_continuous(breaks = seq(1,10,1))+
  scale_y_continuous(breaks = seq(-10,10,5))+
  theme_classic()+
  theme(text = element_text(family = "Arial",face = "bold"),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18))