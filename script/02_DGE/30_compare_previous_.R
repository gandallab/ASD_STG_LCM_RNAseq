library(tidyverse)


# compare to Parikshak
rm(list = ls());gc()
options(stringsAsFactors = F)

load("dge.blocks.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
old=read.delim("parikshak_natrue_cortex.txt")
old1=old[,c(1,9)]
dat=dat[dat$pvalue < 0.05,]
dat1=dat[,c("log2FoldChange"),drop=F]
rownames(dat1)=sapply(rownames(dat1),function(x){str_split_fixed(x,"\\.",Inf)[1,1]})

plotdata=merge(dat1,old1,by.x="row.names",by.y="ENSEMBL.ID")
library(ggplot2)
ggplot(plotdata,aes(log2.FC..ASD.vs.CTL,log2FoldChange))+
  geom_point(shape=1,alpha=0.5)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  xlab("ParikshakNature")+
  ylab("LCM")+
  geom_text(aes(-0.25,-2,label="Blocks vs Parikshak: Spearman r = 0.22"))
cor(plotdata$log2FoldChange,plotdata$log2.FC..ASD.vs.CTL,method = "sp")



# compare among Block, Neuron and Oligo
rm(list = ls());gc()
options(stringsAsFactors = F)

load("dge/dge.blocks.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
#dat=dat[dat$pvalue < 0.05,]
dat1=dat[,c("log2FoldChange","stat","pvalue"),drop=F]
colnames(dat1)=paste(colnames(dat1),"blocks",sep = "_")

load("dge.neuron.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
#dat=dat[dat$pvalue < 0.05,]
dat2=dat[,c("log2FoldChange","stat","pvalue"),drop=F]
colnames(dat2)=paste(colnames(dat2),"neuron",sep = "_")

plotdata=merge(dat1,dat2,by="row.names")

plotdata=plotdata[plotdata$pvalue_blocks < 0.1 | plotdata$pvalue_neuron < 0.1,]
plotdata=plotdata[plotdata$pvalue_blocks < 0.1,]
plotdata=plotdata[plotdata$pvalue_neuron < 0.1,]
plotdata=plotdata[plotdata$pvalue_blocks < 0.1 & plotdata$pvalue_neuron < 0.1,]
plotdata=plotdata[!is.na(plotdata$Row.names),]

library(ggplot2)

ggplot(plotdata,aes(log2FoldChange_blocks,log2FoldChange_neuron))+
  geom_point(shape=1,alpha=0.5)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  ylim(-1.5,2)
  #xlab("Psychencode")+
  #ylab("LCM")+
  geom_text(aes(-1,1.8,label="Neuron vs Block: Spearman r = 0.13"))
#cor(plotdata$log2FoldChange_blocks,plotdata$log2FoldChange_neuron,method = "sp")
cor.test(plotdata$stat_blocks,plotdata$stat_neuron,method = "sp")

ggplot(plotdata,aes(stat_blocks,stat_neuron,color=..count..))+
  geom_hex(bins=60)+
  geom_abline(slope = 1,intercept = 0,linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  xlab("Blocks t statistics")+
  ylab("Neuron t statistics")+
  scale_fill_gradient(low = "gray",high = "magenta")+
  scale_color_gradient(low = "gray",high = "magenta")+
  scale_x_continuous(limits = c(-8,8))+
  scale_y_continuous(limits = c(-8,8))+
  theme_classic(base_size = 16,base_family = "Arial")

ggplot()+
  geom_hex(data=plotdata,aes(stat_blocks,stat_neuron,color=..count..),bins=60)+
  geom_smooth(data = plotdata,aes(stat_blocks,stat_neuron),method = "lm")+
  geom_abline(slope = 1,intercept = 0,linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  xlab("Blocks t statistics")+
  ylab("Neuron t statistics")+
  scale_fill_gradient(low = "gray",high = "magenta")+
  scale_color_gradient(low = "gray",high = "magenta")+
  scale_x_continuous(limits = c(-8,8))+
  scale_y_continuous(limits = c(-8,8))+
  theme_classic(base_size = 16,base_family = "Arial")


# compare to pan-cortical data (Jill)

rm(list = ls());gc()
options(stringsAsFactors = F)
library(tidyverse)
load("dge.blocks.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
old=read.delim("ttable_ASD_BA41_42_22.csv",sep = ",")
#old=old[old$P.Value < 0.1,]
old1=old[,c("X","logFC","t","external_gene_name","P.Value")]
old1$gene_id=sapply(old1$X,function(x){str_split_fixed(x,"\\.",Inf)[1,1]})
#dat=dat[dat$pvalue < 0.1,]
dat1=dat[,c("log2FoldChange","stat","pvalue"),drop=F]
rownames(dat1)=sapply(rownames(dat1),function(x){str_split_fixed(x,"\\.",Inf)[1,1]})
plotdata=merge(dat1,old1,by.x="row.names",by.y="gene_id")

plotdata=plotdata[plotdata$pvalue < 0.1 | plotdata$P.Value < 0.1,]
plotdata=plotdata[plotdata$pvalue < 0.1,]
plotdata=plotdata[plotdata$P.Value < 0.1,]
plotdata=plotdata[plotdata$pvalue < 0.1 & plotdata$P.Value < 0.1,]

ggplot(plotdata,aes(logFC,log2FoldChange))+
  geom_point(shape=1,alpha=0.5)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  xlab("logFC BA_41_42_22")+
  ylab("logFC LCM")+
  geom_text(aes(0.5,-2,label="Blocks vs BA41_42_22: Spearman r = 0.45"),size=4.5)+
  theme(text = element_text(family = "Arial"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size=10))

ggplot(plotdata,aes(t,stat))+
  geom_point(shape=1,alpha=0.5)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  xlab("logFC BA_41_42_22")+
  ylab("logFC LCM")+
  geom_text(aes(0.5,-2,label="Blocks vs BA41_42_22: Spearman r = 0.45"),size=4.5)+
  theme(text = element_text(family = "Arial"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size=10))
cor(plotdata$t,plotdata$stat,method = "sp")
#cor(plotdata$logFC,plotdata$log2FoldChange,method = "sp")

ggplot(plotdata,aes(t,stat,color=..count..))+
  geom_hex(bins=60)+
  geom_abline(slope = 1,intercept = 0,linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  xlab("BA41_42_22 ASD vs CTL t statistics")+
  ylab("STG ASD vs CTL t statistics")+
  scale_fill_gradient(low = "gray",high = "magenta")+
  scale_color_gradient(low = "gray",high = "magenta")+
  scale_x_continuous(limits = c(-8,8))+
  scale_y_continuous(limits = c(-8,8))+
  theme_classic(base_size = 16,base_family = "Arial")


ggplot()+
  geom_hex(data=plotdata,aes(t,stat,color=..count..),bins=60)+
  geom_smooth(data = plotdata,aes(t,stat),method = "lm")+
  geom_abline(slope = 1,intercept = 0,linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  xlab("BA41_42_22 ASD vs CTL t statistics")+
  ylab("STG ASD vs CTL t statistics")+
  scale_fill_gradient(low = "gray",high = "magenta")+
  scale_color_gradient(low = "gray",high = "magenta")+
  scale_x_continuous(limits = c(-8,8))+
  scale_y_continuous(limits = c(-8,8))+
  theme_classic(base_size = 16,base_family = "Arial")



