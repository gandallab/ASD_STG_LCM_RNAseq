rm(list=ls());gc()
options(stringsAsFactors = F)

library(tidyverse)

dat=read.delim("working_data/leafcutter/v2/cca/neuron_unfilter.perind.counts.gz.qqnorm_chr12")
rownames(dat)=dat$ID
dat1=dat[,-c(1:4)]

pca=prcomp(scale(t(dat1),center = T,scale = F))
pcadata=pca$x

datMeta=read.delim("working_data/leafcutter/v2/neuron.group.txt",header = F)

table(rownames(pcadata) == datMeta$V1)

plotdata=merge(datMeta,pcadata,by.x="V1",by.y="row.names")

ggplot(plotdata)+
  geom_point(aes(PC1,PC2,col=V2,shape=V3),size=4)+
  #  geom_text_repel(aes(PC1,PC2,label=rownames(dat)),size=4)+
  theme_classic()+
  theme(text = element_text(family = "Arial"))

X = model.matrix(~ V2 + V3 + V4 + V5 + V6 + V7 + V8,datMeta)
Y = dat1
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
#datExpr = Y - t(as.matrix(X[,5]) %*% t(as.matrix(beta[5,])))
datSp = Y - t(X[,c(5:ncol(X))] %*% beta[c(5:nrow(beta)),])
table(datMeta$V1 == colnames(datSp))

load("working_data/wgcna/voom.forWGCNA.input.neuron.fc.manualRemove.RData")
dim(datExpr)
table(make.names(colnames(datExpr)) %in% colnames(datSp))
colnames(datExpr)=make.names(colnames(datExpr))
datExpr1=datExpr[,colnames(datSp)]
datExpr2=t(datExpr1)

datSp1=t(datSp)

write.table(datExpr2,file = "working_data/cca/chr12_geneexpr.txt",quote = F,sep = "\t",row.names = T)
write.table(datSp1,file = "working_data/cca/chr12_sp.txt",quote = F,sep = "\t",row.names = T)

library(CCA)

res.regul <- estim.regul(datExpr2, datSp1, plt = TRUE,grid1 = seq(0.05,0.5,l=10),grid2 = seq(0.05,0.5,l=10))
