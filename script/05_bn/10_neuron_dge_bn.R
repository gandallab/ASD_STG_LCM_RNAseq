rm(list = ls());gc()
options(stringsAsFactors = F)
#library(tidyverse)
library(bnlearn)

# DGE result

load("working_data/dge/dge.neuron.fc2pass.RData")
dat=as.data.frame(dge.deseq2)


sigGenes=dat[dat$padj < 0.02,]
sigGenes=sigGenes[!is.na(sigGenes$padj),]

# WGCNA inputs

load("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/neuron.fc2pass.recut.RData")

datExpr=networks$datExpr
inputdata=datExpr[rownames(sigGenes),]
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
save(inputdata_disc,file = "working_data/bn/neuron_dge_bn_input_fdr002.RData")

load("working_data/bn/neuron_dge_bn_input_fdr002.RData")
bl=data.frame(from="ASD",to=colnames(inputdata_disc))
boot = boot.strength(data = inputdata_disc, R = 50, algorithm = "hc",
                     algorithm.args = list(score = "bde", iss = 10,blacklist=bl))
save(boot,file = "working_data/bn/neuron_dge_bn_fdr002_bootstrap50_blacklist.RData")

load("working_data/bn/neuron_dge_bn_fdr002_bootstrap50.RData")

boot[(boot$strength > 0.7) & (boot$direction >= 0.5), ]
avg.boot = averaged.network(boot,threshold = 0.5)
graphviz.plot(avg.boot)

# random start

nodes = names(inputdata_disc)
start = random.graph(nodes = nodes, method = "melancon", num = 50)
netlist = lapply(start, function(net) {
  hc(inputdata_disc, score = "bde", iss = 10, start = net) })
rnd = custom.strength(netlist, nodes = nodes)
rnd[(rnd$strength > 0.85) & (rnd$direction >= 0.5), ]
avg.start = averaged.network(rnd, threshold = 0.85)
