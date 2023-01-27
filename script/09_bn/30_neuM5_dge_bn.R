rm(list = ls());gc()
options(stringsAsFactors = F)
library(tidyverse)
library(bnlearn)

# WGCNA inputs

load("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/neuron.fc2pass.recut.RData")

modules = networks$merged$colors
m5genes=rownames(networks$datExpr)[modules == 5]
# m5genes=sapply(m5genes,function(x) str_split_fixed(x,"\\.",Inf)[1,1])
# annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt",header = T,sep = "\t")
# m5genename=annot$genename[match(m5genes,annot$geneid)]

load("working_data/dge/dge.neuron.fc2pass.RData")
dat=as.data.frame(dge.deseq2)


sigGenes=dat[dat$padj < 0.05,]
sigGenes=sigGenes[!is.na(sigGenes$padj),]
# rownames(sigGenes)=sapply(rownames(sigGenes),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
# sigGenes$genename=annot$genename[match(rownames(sigGenes),annot$geneid)]

finalgenes=intersect(m5genes,rownames(sigGenes))


datExpr=networks$datExpr
inputdata=datExpr[finalgenes,]
rownames(inputdata)=sapply(rownames(inputdata),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
inputdata=as.data.frame(inputdata)
#sachs = read.table("C:/Users/Pan Zhang/Downloads/sachs.interventional.txt/sachs.data.txt", header = TRUE)

annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt",header = T,sep = "\t")
inputdata$genename=annot$genename[match(rownames(inputdata),annot$geneid)]
rownames(inputdata)=inputdata$genename
inputdata$genename=NULL
inputdata=as.data.frame(t(inputdata))
datMeta=networks$datMeta
inputdata$ASD=datMeta$Diagnosis[match(rownames(inputdata),rownames(datMeta))]
inputdata$ASD=factor(inputdata$ASD,levels = c("Control","Autism"))
inputdata_disc = discretize(inputdata, method = "hartemink", breaks = 3, ibreaks = 60, idisc = "quantile")
save(inputdata_disc,file = "working_data/bn/neuM5_dge_bn_input.RData")

load("working_data/bn/neuM5_dge_bn_input.RData")
# bl=data.frame(from="ASD",to=colnames(inputdata_disc))
# bl=bl[bl$to != "ASD",]
inputdata_disc$ASD=NULL
boot = boot.strength(data = inputdata_disc, R = 500, algorithm = "hc",
                     algorithm.args = list(score = "bde", iss = 10))
save(boot,file = "working_data/bn/neuM5_dge_bn_bootstrap500_noASD.RData")

load("working_data/bn/neuM5_dge_bn_bootstrap500.RData")
boot[(boot$strength > 0.7) & (boot$direction >= 0.5), ]
avg.boot = averaged.network(boot,threshold = 0.6) # 0.6/0.57 gives a connected graph
graphviz.plot(avg.boot,layout = "dot",shape = "ellipse")
score(avg.boot, inputdata_disc, type = "bde", iss = 10)
arcs=as.data.frame(avg.boot$arcs)
arcs=merge(arcs,boot,by=c("from","to"))
write.table(arcs,file = "working_data/bn/neuM5_dge_bn_arcTable_withASD.txt",quote = F,sep = "\t",row.names = F)
