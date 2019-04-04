rm(list=ls());gc()
options(stringsAsFactors = F)

library(ggplot2)
library(SummarizedExperiment)


load("working_data/summarizedExperiment/se_blocks_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/se_neuron_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/se_oligo_CPM_outlierRemoved.RData")

se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
se_oligo$Diagnosis = factor(se_oligo$Diagnosis, levels= c("Control", "Autism"))
se_neuron$Diagnosis = factor(se_neuron$Diagnosis, levels= c("Control", "Autism"))

# blocks
datMeta=as.data.frame(colData(se_blocks))
X = model.matrix(~Diagnosis + Sex + Age + seqPC2 + seqPC3,colData(se_blocks))
Y = assays(se_blocks)$log2CPM
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
#datExpr = Y - t(as.matrix(X[,4]) %*% t(as.matrix(beta[4,])))
datExpr = Y - t(X[,c(5:ncol(X))] %*% beta[c(5:nrow(beta)),])
table(rownames(datMeta) == colnames(datExpr))
save(file = "./working_data/wgcna/voom.forWGCNA.input.blocks.RData",datExpr,datMeta)

# neuron
datMeta=as.data.frame(colData(se_neuron))
datExpr = assays(se_neuron)$log2CPM
table(rownames(datMeta) == colnames(datExpr))
save(file = "./working_data/wgcna/voom.forWGCNA.input.neuron.RData",datExpr,datMeta)

# oligo
datMeta=as.data.frame(colData(se_oligo))
datExpr = assays(se_oligo)$log2CPM
table(rownames(datMeta) == colnames(datExpr))
save(file = "./working_data/wgcna/voom.forWGCNA.input.oligo.RData",datExpr,datMeta)
