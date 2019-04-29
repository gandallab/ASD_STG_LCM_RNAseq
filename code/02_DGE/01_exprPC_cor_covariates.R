rm(list=ls());gc()
options(stringsAsFactors = F)

library(edgeR)
library(tidyverse)
library(SummarizedExperiment)
library(caret)
library(export)


# blocks
load("working_data/summarizedExperiment/se_blocks_CPM_outlierRemoved.RData")

datMeta=Filter(function(x)(length(unique(x))>1), as.data.frame(colData(se_blocks)))
missing = colSums(is.na(datMeta))
datMeta=datMeta[,missing == 0]
ind=sapply(datMeta,function(x){is.character(x) & length(unique(x)) == nrow(datMeta)})
datMeta=datMeta[,!ind]
candidates="Diagnosis , Age , Sex , Type_RNAseqRunNumber , RNAseqPool , RIN , X260.280 , X260.230 , PositionOnSequencingPlate , seqPC1 , seqPC2 , seqPC3 , seqPC4 , seqPC5 , seqPC6 , seqPC7 , seqPC8 ,seqPC9 , seqPC10"
#candidates=sapply(strsplit(candidates, '[, ]+'), function(x) toString(dQuote(x)))
datMeta1=datMeta[,colnames(datMeta) %in% strsplit(candidates, '[, ]+')[[1]]]

exPCA=prcomp(scale(t(assays(se_blocks)$log2CPM),center = TRUE, scale = FALSE))
exPCAdata=exPCA$x
exPCAdata=exPCAdata[,1:10]

# exPCA1=cmdscale(dist(t(assays(se_blocks)$log2CPM)),k = 10,eig = TRUE)
# exPCAdata1=exPCA1$points
#comb=expand.grid(colnames(exPCAdata),colnames(datMeta1))
cross=data.frame()
for (i in 1:ncol(exPCAdata)){
  for (j in 1:ncol(datMeta1)){
    vec1=exPCAdata[,i]
    vec2=datMeta1[,j]
    mylm=lm(vec1 ~ vec2)
    cross=rbind(cross,data.frame(exPC=colnames(exPCAdata)[i],covariate=colnames(datMeta1)[j],r.squared=summary(mylm)$adj.r.squared))
  }
}

cross$correlation=sqrt(cross$r.squared)
cross$correlation[is.na(cross$correlation)]=0
cross$exPC=factor(cross$exPC,levels = colnames(exPCAdata))
cross$covariate=factor(cross$covariate,levels = colnames(datMeta1))
ggplot(cross,aes(exPC,covariate,fill=correlation))+
  geom_tile()+
#  scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = 0.3)+
  scale_fill_gradient(low = "white",high = "red")

finalCandidates=c("RIN","X260.280","seqPC2","seqPC3")
datMeta2=datMeta1[,finalCandidates]
cor(datMeta2) # seqPC2 is correlated to both RIN and X260.280
graph2ppt(file="working_data/figures/blocks.exPC.covariates.correlation.pptx",width=10.5,height=7.1,append=FALSE)


# neuron
rm(list=ls());gc()
load("working_data/summarizedExperiment/se_neuron_CPM_outlierRemoved.RData")

datMeta=Filter(function(x)(length(unique(x))>1), as.data.frame(colData(se_neuron)))
missing = colSums(is.na(datMeta))
datMeta=datMeta[,missing == 0]
ind=sapply(datMeta,function(x){is.character(x) & length(unique(x)) == nrow(datMeta)})
datMeta=datMeta[,!ind]
candidates="Diagnosis , Age , Sex , Type_RNAseqRunNumber , RNAseqPool , PositionOnSequencingPlate , seqPC1 , seqPC2 , seqPC3 , seqPC4 , seqPC5 , seqPC6 , seqPC7 , seqPC8 ,seqPC9 , seqPC10"
#candidates=sapply(strsplit(candidates, '[, ]+'), function(x) toString(dQuote(x)))
datMeta1=datMeta[,colnames(datMeta) %in% strsplit(candidates, '[, ]+')[[1]]]
exPCA=prcomp(scale(t(assays(se_neuron)$log2CPM),center = TRUE, scale = FALSE))
exPCAdata=exPCA$x
exPCAdata=exPCAdata[,1:10]
#comb=expand.grid(colnames(exPCAdata),colnames(datMeta1))
cross=data.frame()
for (i in 1:ncol(exPCAdata)){
  for (j in 1:ncol(datMeta1)){
    vec1=exPCAdata[,i]
    vec2=datMeta1[,j]
    mylm=lm(vec1 ~ vec2)
    cross=rbind(cross,data.frame(exPC=colnames(exPCAdata)[i],covariate=colnames(datMeta1)[j],r.squared=summary(mylm)$adj.r.squared))
  }
}

cross$correlation=sqrt(cross$r.squared)
cross$correlation[is.na(cross$correlation)]=0
cross$exPC=factor(cross$exPC,levels = colnames(exPCAdata))
cross$covariate=factor(cross$covariate,levels = colnames(datMeta1))
ggplot(cross,aes(exPC,covariate,fill=correlation))+
  geom_tile()+
  #  scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = 0.3)+
  scale_fill_gradient(low = "white",high = "red")
   ## seqPC2,seqPC3
graph2ppt(file="working_data/figures/neuron.exPC.covariates.correlation.pptx",width=10.5,height=7.1,append=FALSE)

# oligo
rm(list=ls());gc()
load("working_data/summarizedExperiment/se_oligo_CPM_outlierRemoved.RData")

datMeta=Filter(function(x)(length(unique(x))>1), as.data.frame(colData(se_oligo)))
missing = colSums(is.na(datMeta))
datMeta=datMeta[,missing == 0]
ind=sapply(datMeta,function(x){is.character(x) & length(unique(x)) == nrow(datMeta)})
datMeta=datMeta[,!ind]
candidates="Diagnosis , Age , Sex , Type_RNAseqRunNumber , RNAseqPool , PositionOnSequencingPlate , seqPC1 , seqPC2 , seqPC3 , seqPC4 , seqPC5 , seqPC6 , seqPC7 , seqPC8 ,seqPC9 , seqPC10"
#candidates=sapply(strsplit(candidates, '[, ]+'), function(x) toString(dQuote(x)))
datMeta1=datMeta[,colnames(datMeta) %in% strsplit(candidates, '[, ]+')[[1]]]
exPCA=prcomp(scale(t(assays(se_oligo)$log2CPM),center = TRUE, scale = FALSE))
exPCAdata=exPCA$x
exPCAdata=exPCAdata[,1:10]
#comb=expand.grid(colnames(exPCAdata),colnames(datMeta1))
cross=data.frame()
for (i in 1:ncol(exPCAdata)){
  for (j in 1:ncol(datMeta1)){
    vec1=exPCAdata[,i]
    vec2=datMeta1[,j]
    mylm=lm(vec1 ~ vec2)
    cross=rbind(cross,data.frame(exPC=colnames(exPCAdata)[i],covariate=colnames(datMeta1)[j],r.squared=summary(mylm)$adj.r.squared,pvalue=summary(mylm)$coefficients[2,4]))
  }
}

cross$correlation=sqrt(cross$r.squared)
cross$correlation[is.na(cross$correlation)]=0
cross$exPC=factor(cross$exPC,levels = colnames(exPCAdata))
cross$covariate=factor(cross$covariate,levels = colnames(datMeta1))
ggplot(cross,aes(exPC,covariate,fill=correlation))+
  geom_tile()+
  #  scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = 0.3)+
  scale_fill_gradient(low = "white",high = "red")
  ## seqPC2
graph2ppt(file="working_data/figures/oligo.exPC.covariates.correlation.pptx",width=10.5,height=7.1,append=FALSE)
