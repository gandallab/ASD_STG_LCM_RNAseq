library(tidyverse)
rm(list = ls());gc()
options(stringsAsFactors = F)


# compare to MSBB

rm(list = ls());gc()
options(stringsAsFactors = F)
library(tidyverse)
#load("dge.blocks.fc2pass.RData")
load("dge.neuron.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
old=read.delim("MSBB_Differential_Expression_(diagnosis).tsv")
old1=old[old$Comparison == "AD_STG - CT_STG", c("ensembl_gene_id","logFC","t","P.Value","adj.P.Val","Direction","hgnc_symbol")]
#dat=dat[dat$pvalue < 0.1,]
dat1=dat[,c("log2FoldChange","stat","pvalue","padj"),drop=F]
rownames(dat1)=sapply(rownames(dat1),function(x){str_split_fixed(x,"\\.",Inf)[1,1]})
plotdata=merge(dat1,old1,by.x="row.names",by.y="ensembl_gene_id")
plotdata=plotdata[!is.na(plotdata$Row.names),]

plotdata=plotdata[plotdata$pvalue < 0.1 | plotdata$P.Value < 0.1,]
plotdata=plotdata[plotdata$padj < 0.1,]
plotdata=plotdata[plotdata$adj.P.Val < 0.1,]
plotdata=plotdata[plotdata$pvalue < 0.1 & plotdata$P.Value < 0.1,]
plotdata=plotdata[plotdata$padj < 0.1 & plotdata$adj.P.Val < 0.1,]
plotdata=plotdata[!is.na(plotdata$Row.names),]

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
cor(plotdata$logFC,plotdata$log2FoldChange,method = "sp")

ggplot(plotdata,aes(stat,t,color=..count..))+
  geom_hex(bins=60)+
  geom_abline(slope = 1,intercept = 0,linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  xlab("ASD DGE t statistics")+
  ylab("AD DGE t statistics")+
  scale_fill_gradient(low = "gray",high = "magenta")+
  scale_color_gradient(low = "gray",high = "magenta")+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(-10,10))+
  theme_classic(base_size = 16,base_family = "Arial")


ggplot()+
  geom_hex(data=plotdata,aes(t,stat,color=..count..),bins=60)+
  geom_smooth(data = plotdata,aes(t,stat),method = "lm")+
  geom_vline(xintercept = 0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  xlab("ASD DGE t statistics")+
  ylab("AD DGE t statistics")+
  scale_fill_gradient(low = "gray",high = "magenta")+
  scale_color_gradient(low = "gray",high = "magenta")+
  scale_x_continuous(limits=c(-10,10))+
  scale_y_continuous(limits = c(-10,10))+
  theme_classic(base_size = 16,base_family = "Arial")

ggvenn(data=list(ASD=plotdata$Row.names[plotdata$padj < 0.1],AD=plotdata$Row.names[plotdata$adj.P.Val < 0.1]),show_percentage=FALSE)

