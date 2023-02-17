rm(list = ls());gc()
options(stringsAsFactors = F)
library(tidyverse)
library(bnlearn)

# WGCNA inputs

load("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/neuron.fc2pass.recut.RData")

inputdata=networks$MEs$eigengenes

datMeta=networks$datMeta

inputdata$ASD=datMeta$Diagnosis[match(rownames(inputdata),rownames(datMeta))]
inputdata$ASD=factor(inputdata$ASD,levels = c("Control","Autism"))
inputdata$ME0=NULL

inputdata_disc = discretize(inputdata, method = "hartemink", breaks = 3, ibreaks = 60, idisc = "quantile")
save(inputdata_disc,file = "working_data/bn/neuron_eigengene_bn_input.RData")

boot = boot.strength(data = inputdata_disc, R = 500, algorithm = "hc",
                     algorithm.args = list(score = "bde", iss = 10))
save(boot,file = "working_data/bn/neuron_eigengene_bn_bootstrap500.RData")



load("working_data/bn/neuron_eigengene_bn_input.RData")
load("working_data/bn/neuron_eigengene_bn_bootstrap500.RData")

boot[(boot$strength > 0.7) & (boot$direction >= 0.5), ]
avg.boot = averaged.network(boot,threshold = 0.67) # 0.67 gives a connected graph
graphviz.plot(avg.boot,layout = "dot",shape = "ellipse")
score(avg.boot, inputdata_disc, type = "bde", iss = 10)
arcs=as.data.frame(avg.boot$arcs)
arcs=merge(arcs,boot,by=c("from","to"))

write.table(arcs,file = "working_data/bn/neuron_bn_eigengene_arcTable.txt",quote = F,sep = "\t",row.names = F)

library(WGCNA)
maptable=data.frame(node=paste("Neu-M",1:18,sep = ""),color=labels2colors(1:18))
write.table(maptable,file = "working_data/bn/neuron_bn_eigengene_colorTable.txt",quote = F,sep = "\t",row.names = F)
cv5=bn.cv(inputdata_disc,k=5,bn = avg.boot,loss = "pred-lw",loss.args = list(target = "ASD"))
dat=data.frame(OBS = unlist(lapply(cv5, `[[`, "observed")),
PRED = unlist(lapply(cv5, `[[`, "predicted")))
with(dat,table(OBS, PRED))
