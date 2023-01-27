rm(list=ls());gc()
options(stringsAsFactors = F)

library(edgeR)
library(ggplot2)
library(SummarizedExperiment)
library(caret)


# blocks
load("working_data/summarizedExperiment/se_blocks_CPM_outlierRemoved.RData")

datMeta=Filter(function(x)(length(unique(x))>1), as.data.frame(colData(se_blocks)))
missing = colSums(is.na(datMeta))
datMeta=datMeta[,missing == 0]
ind=sapply(datMeta,function(x){is.character(x) & length(unique(x)) == nrow(datMeta)})
datMeta=datMeta[,!ind]
X=model.matrix(~0 + Diagnosis + Age + Sex + Type_RNAseqRunNumber + RNAseqPool + 
                 RIN + X260.280 + X260.230 + PositionOnSequencingPlate + seqPC1 + 
                 seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + 
                 seqPC9 + seqPC10,data = datMeta)
paircorr=findCorrelation(cor(X),cutoff = 0.9,names = T)
X=X[,!colnames(X) %in% paircorr]


source("D:/scripts/local/scripts//MARS (2).R")
e1 = runMARS(gene.expr = assays(se_blocks)$log2CPM, covars = X,n.cores = 1, n.replicates = 100, allow.interaction = F)
dat = rbind(e1)
idx = which(dat$var %in% names(table(dat$var))[table(dat$var) >= mean(table(dat$var))])
ggplot(dat[idx,], aes(x=reorder(var, abs.beta), y=abs.beta)) + geom_boxplot() + coord_flip() + xlab("")+theme_classic() # seqPC2,seqPC3

## add SVs
rm(list=ls());gc()
options(stringsAsFactors = F)
library(sva)
load("working_data/summarizedExperiment/se_blocks_CPM_outlierRemoved.RData")
se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
mod=model.matrix(~Diagnosis + Age + Sex + seqPC2 + seqPC3,data=colData(se_blocks))
mod0=model.matrix(~Age + Sex + seqPC2 + seqPC3,data=colData(se_blocks))
svaResult=svaseq(cpm(calcNormFactors(DGEList(assays(se_blocks)$counts), method = 'TMM'), log = F),mod = mod,mod0 = mod0)
sv=svaResult$sv
colnames(sv)=paste("sv",seq(1,ncol(sv)),sep = "")
colData(se_blocks)=cbind(colData(se_blocks),sv)

datMeta=Filter(function(x)(length(unique(x))>1), as.data.frame(colData(se_blocks)))
missing = colSums(is.na(datMeta))
datMeta=datMeta[,missing == 0]
ind=sapply(datMeta,function(x){is.character(x) & length(unique(x)) == nrow(datMeta)})
datMeta=datMeta[,!ind]
X=model.matrix(~0 + Diagnosis + Age + Sex + Type_RNAseqRunNumber + RNAseqPool + 
                 RIN + X260.280 + X260.230 + PositionOnSequencingPlate + seqPC1 + 
                 seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + 
                 seqPC9 + seqPC10 + sv1+sv2+sv3+sv4+sv5+sv6+sv7+sv8+sv9+sv10+sv11,data = datMeta)
paircorr=findCorrelation(cor(X),cutoff = 0.9,names = T)
linearcomb=findLinearCombos(X) 
colnames(X)[linearcomb$remove] # all SVs are linear combindations of all existing variates

# neurons
rm(list=ls());gc()
options(stringsAsFactors = F)

load("working_data/summarizedExperiment/se_neuron_CPM_outlierRemoved.RData")

datMeta=Filter(function(x)(length(unique(x))>1), as.data.frame(colData(se_neuron)))
missing = colSums(is.na(datMeta))
datMeta=datMeta[,missing == 0]
ind=sapply(datMeta,function(x){is.character(x) & length(unique(x)) == nrow(datMeta)})
datMeta=datMeta[,!ind]
X=model.matrix(~0 + Diagnosis + Age + Sex + Type_RNAseqRunNumber + RNAseqPool + 
                 PositionOnSequencingPlate + seqPC1 + 
                 seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + 
                 seqPC9 + seqPC10,data = datMeta)
paircorr=findCorrelation(cor(X),cutoff = 0.9,names = T)
X=X[,!colnames(X) %in% paircorr]


source("D:/scripts/local/scripts//MARS (2).R")
e1 = runMARS(gene.expr = assays(se_neuron)$log2CPM, covars = X,n.cores = 1, n.replicates = 100, allow.interaction = F)
dat = rbind(e1)
idx = which(dat$var %in% names(table(dat$var))[table(dat$var) >= mean(table(dat$var))])
ggplot(dat[idx,], aes(x=reorder(var, abs.beta), y=abs.beta)) + geom_boxplot() + coord_flip() + xlab("")+theme_classic() # just intercept


## add SVs
rm(list=ls());gc()
options(stringsAsFactors = F)
library(sva)
load("working_data/summarizedExperiment/se_neuron_CPM_outlierRemoved.RData")
se_neuron$Diagnosis = factor(se_neuron$Diagnosis, levels= c("Control", "Autism"))
mod=model.matrix(~Diagnosis + Age + Sex,data=colData(se_neuron))
mod0=model.matrix(~Age + Sex,data=colData(se_neuron))
svaResult=svaseq(cpm(calcNormFactors(DGEList(assays(se_neuron)$counts), method = 'TMM'), log = F),mod = mod,mod0 = mod0)
sv=svaResult$sv
colnames(sv)=paste("sv",seq(1,ncol(sv)),sep = "")
colData(se_neuron)=cbind(colData(se_neuron),sv)

datMeta=Filter(function(x)(length(unique(x))>1), as.data.frame(colData(se_neuron)))
missing = colSums(is.na(datMeta))
datMeta=datMeta[,missing == 0]
ind=sapply(datMeta,function(x){is.character(x) & length(unique(x)) == nrow(datMeta)})
datMeta=datMeta[,!ind]
X=model.matrix(~0 + Diagnosis + Age + Sex + Type_RNAseqRunNumber + RNAseqPool + 
                 PositionOnSequencingPlate + seqPC1 + 
                 seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + 
                 seqPC9 + seqPC10 + sv1+sv2+sv3+sv4+sv5+sv6+sv7+sv8+sv9+sv10+sv11+sv12+sv13,data = datMeta)
paircorr=findCorrelation(cor(X),cutoff = 0.9,names = T)
linearcomb=findLinearCombos(X) 
colnames(X)[linearcomb$remove] # all SVs are linear combindations of all existing variates

# oligo
rm(list=ls());gc()
options(stringsAsFactors = F)

load("working_data/summarizedExperiment/se_oligo_CPM_outlierRemoved.RData")

datMeta=Filter(function(x)(length(unique(x))>1), as.data.frame(colData(se_oligo)))
missing = colSums(is.na(datMeta))
datMeta=datMeta[,missing == 0]
ind=sapply(datMeta,function(x){is.character(x) & length(unique(x)) == nrow(datMeta)})
datMeta=datMeta[,!ind]
X=model.matrix(~0 + Diagnosis + Age + Sex + RNAseqPool + 
                 PositionOnSequencingPlate + seqPC1 + 
                 seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + 
                 seqPC9 + seqPC10,data = datMeta)
paircorr=findCorrelation(cor(X),cutoff = 0.9,names = T)
X=X[,!colnames(X) %in% paircorr]


source("D:/scripts/local/scripts//MARS (2).R")
e1 = runMARS(gene.expr = assays(se_oligo)$log2CPM, covars = X,n.cores = 1, n.replicates = 100, allow.interaction = F)
dat = rbind(e1)
idx = which(dat$var %in% names(table(dat$var))[table(dat$var) >= mean(table(dat$var))])
ggplot(dat[idx,], aes(x=reorder(var, abs.beta), y=abs.beta)) + geom_boxplot() + coord_flip() + xlab("")+theme_classic() # only intercept


## add SVs
rm(list=ls());gc()
options(stringsAsFactors = F)
library(sva)
load("working_data/summarizedExperiment/se_oligo_CPM_outlierRemoved.RData")
se_oligo$Diagnosis = factor(se_oligo$Diagnosis, levels= c("Control", "Autism"))
mod=model.matrix(~Diagnosis + Age + Sex,data=colData(se_oligo))
mod0=model.matrix(~Age + Sex,data=colData(se_oligo))
svaResult=svaseq(cpm(calcNormFactors(DGEList(assays(se_oligo)$counts), method = 'TMM'), log = F),mod = mod,mod0 = mod0)
sv=svaResult$sv
colnames(sv)=paste("sv",seq(1,ncol(sv)),sep = "")
colData(se_oligo)=cbind(colData(se_oligo),sv)

datMeta=Filter(function(x)(length(unique(x))>1), as.data.frame(colData(se_oligo)))
missing = colSums(is.na(datMeta))
datMeta=datMeta[,missing == 0]
ind=sapply(datMeta,function(x){is.character(x) & length(unique(x)) == nrow(datMeta)})
datMeta=datMeta[,!ind]
X=model.matrix(~0 + Diagnosis + Age + Sex + RNAseqPool + 
                 PositionOnSequencingPlate + seqPC1 + 
                 seqPC2 + seqPC3 + seqPC4 + seqPC5 + seqPC6 + seqPC7 + seqPC8 + 
                 seqPC9 + seqPC10 + sv1+sv2+sv3+sv4+sv5,data = datMeta)
paircorr=findCorrelation(cor(X),cutoff = 0.9,names = T)
linearcomb=findLinearCombos(X) 
colnames(X)[linearcomb$remove]# all SVs are linear combindations of all existing variates

### Final model
blocks: ~Diagnosis + Age + Sex + seqPC2 + seqPC3
neuron: ~Diagnosis + Age + Sex
oligo: ~Diagnosis + Age + Sex
