rm(list=ls());gc()
options(stringsAsFactors = F)

library(ggplot2)
library(SummarizedExperiment)


load("working_data/summarizedExperiment/fc/se_blocks_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/fc/se_neuron_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/fc/se_oligo_CPM_outlierRemoved.RData")

se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
se_oligo$Diagnosis = factor(se_oligo$Diagnosis, levels= c("Control", "Autism"))
se_neuron$Diagnosis = factor(se_neuron$Diagnosis, levels= c("Control", "Autism"))

se_neuron=se_neuron[,!se_neuron$Sample %in% c("AO_B02","6221_N")]

# blocks
datMeta=as.data.frame(colData(se_blocks))
X = model.matrix(~Diagnosis + Sex + Age + RIN + X260.280 + seqPC1 + seqPC2 + seqPC3,colData(se_blocks))
Y = assays(se_blocks)$log2CPM
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
#datExpr = Y - t(as.matrix(X[,4]) %*% t(as.matrix(beta[4,])))
datExpr = Y - t(X[,c(5:ncol(X))] %*% beta[c(5:nrow(beta)),])
table(rownames(datMeta) == colnames(datExpr))
save(file = "./working_data/wgcna/voom.forWGCNA.input.blocks.fc.RData",datExpr,datMeta)

# neuron
datMeta=as.data.frame(colData(se_neuron))
X = model.matrix(~Diagnosis + Sex + Age + Type_RNAseqRunNumber + seqPC1 + seqPC2 + seqPC3,colData(se_neuron))
Y = assays(se_neuron)$log2CPM
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
#datExpr = Y - t(as.matrix(X[,5]) %*% t(as.matrix(beta[5,])))
datExpr = Y - t(X[,c(5:ncol(X))] %*% beta[c(5:nrow(beta)),])
table(rownames(datMeta) == colnames(datExpr))
save(file = "./working_data/wgcna/voom.forWGCNA.input.neuron.fc2pass.RData",datExpr,datMeta)
mds=cmdscale(dist(t(datExpr)))
colnames(mds)=c("PC1","PC2")
dat=cbind(mds,datMeta)

ggplot(dat)+
  geom_point(aes(PC1,PC2,col=Diagnosis,shape=Sex),size=4)+
  #  geom_text_repel(aes(PC1,PC2,label=rownames(dat)),size=4)+
  theme_classic()+
  theme(text = element_text(family = "Arial"))
pca=prcomp(scale(t(datExpr),center = T,scale = F))
pcadata=pca$x
pcarot=pca$rotation
# oligo
datMeta=as.data.frame(colData(se_oligo))
X = model.matrix(~Diagnosis + Sex + Age + seqPC1 + seqPC2 + seqPC3,colData(se_oligo))
Y = assays(se_oligo)$log2CPM
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
#datExpr = Y - t(as.matrix(X[,4]) %*% t(as.matrix(beta[4,])))
datExpr = Y - t(X[,c(5:ncol(X))] %*% beta[c(5:nrow(beta)),])
table(rownames(datMeta) == colnames(datExpr))
save(file = "./working_data/wgcna/voom.forWGCNA.input.oligo.fc.RData",datExpr,datMeta)



# sft
rm(list=ls());gc()
load("working_data/wgcna/voom.forWGCNA.input.neuron.fc.manualRemove.RData")
powers = c(seq(1,9,by=1),seq(10,30,by=2))
sft = pickSoftThreshold(data= t(datExpr), networkType = "signed", corFnc="bicor",verbose=2,powerVector=powers,blockSize = 30000)

par(mfrow=c(2,1))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2")
abline(h=0.8, col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")
