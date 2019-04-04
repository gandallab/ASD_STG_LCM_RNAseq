rm(list=ls());gc()
options(stringsAsFactors = F)


library(edgeR)
library(ggplot2)
library(SummarizedExperiment)


load("working_data/summarizedExperiment/se_blocks_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/se_neuron_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/se_oligo_CPM_outlierRemoved.RData")

se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
se_oligo$Diagnosis = factor(se_oligo$Diagnosis, levels= c("Control", "Autism"))
se_neuron$Diagnosis = factor(se_neuron$Diagnosis, levels= c("Control", "Autism"))


# blocks
y<-DGEList(counts=assays(se_blocks)$counts)
y <- calcNormFactors(y)

design=model.matrix(~Diagnosis + Sex + Age + seqPC2 + seqPC3,data = colData(se_blocks)) # removing RIN as a covariate increases the number of significant genes
y <- estimateDisp(y, design,robust = TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit, coef=2)
dge.edger=topTags(qlf,n = Inf)
save(file = "./working_data/dge/dge.blocks.edger.RData",dge.edger)
dge.table=dge.edger$table
table(dge.table$FDR < 0.05)

# neuron
y<-DGEList(counts=assays(se_neuron)$counts)
y <- calcNormFactors(y)

design=model.matrix(~Diagnosis + Sex + Age,data = colData(se_neuron)) # removing RIN as a covariate increases the number of significant genes
y <- estimateDisp(y, design,robust = TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
coef(fit)
qlf <- glmQLFTest(fit, coef=2)
dge.edger=topTags(qlf,n = Inf)
save(file = "./working_data/dge/dge.neuron.edger.RData",dge.edger)
dge.table=dge.edger$table
table(dge.table$FDR < 0.05)


# oligo
y<-DGEList(counts=assays(se_oligo)$counts)
y <- calcNormFactors(y)

design=model.matrix(~Diagnosis + Sex + Age,data = colData(se_oligo)) # removing RIN as a covariate increases the number of significant genes
y <- estimateDisp(y, design,robust = TRUE)
fit <- glmQLFit(y, design, robust=TRUE)
coef(fit)
qlf <- glmQLFTest(fit, coef=2)
dge.edger=topTags(qlf,n = Inf)
save(file = "./working_data/dge/dge.oligo.edger.RData",dge.edger)
dge.table=dge.edger$table
table(dge.table$FDR < 0.05)
