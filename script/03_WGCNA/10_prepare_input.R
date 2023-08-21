rm(list=ls());gc()
options(stringsAsFactors = F)

library(ggplot2)
library(SummarizedExperiment)


load("se_blocks.RData")
load("se_neuron.RData")


se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
se_neuron$Diagnosis = factor(se_neuron$Diagnosis, levels= c("Control", "Autism"))



# blocks
datMeta=as.data.frame(colData(se_blocks))
X = model.matrix(~Diagnosis + Sex + Age + RIN + X260.280 + seqPC1 + seqPC2 + seqPC3,colData(se_blocks))
Y = assays(se_blocks)$log2CPM
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
#datExpr = Y - t(as.matrix(X[,4]) %*% t(as.matrix(beta[4,])))
datExpr = Y - t(X[,c(5:ncol(X))] %*% beta[c(5:nrow(beta)),])
table(rownames(datMeta) == colnames(datExpr))
save(file = "voom.forWGCNA.input.blocks.fc.RData",datExpr,datMeta)

# neuron
datMeta=as.data.frame(colData(se_neuron))
X = model.matrix(~Diagnosis + Sex + Age + Type_RNAseqRunNumber + seqPC1 + seqPC2 + seqPC3,colData(se_neuron))
Y = assays(se_neuron)$log2CPM
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
#datExpr = Y - t(as.matrix(X[,5]) %*% t(as.matrix(beta[5,])))
datExpr = Y - t(X[,c(5:ncol(X))] %*% beta[c(5:nrow(beta)),])
table(rownames(datMeta) == colnames(datExpr))
save(file = "voom.forWGCNA.input.neuron.fc2pass.RData",datExpr,datMeta)
