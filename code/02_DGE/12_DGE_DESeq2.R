rm(list=ls());gc()
options(stringsAsFactors = F)

library(DESeq2)
library(ggplot2)
library(SummarizedExperiment)


load("working_data/summarizedExperiment/se_blocks_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/se_neuron_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/se_oligo_CPM_outlierRemoved.RData")

se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
se_oligo$Diagnosis = factor(se_oligo$Diagnosis, levels= c("Control", "Autism"))
se_neuron$Diagnosis = factor(se_neuron$Diagnosis, levels= c("Control", "Autism"))


# blocks
datMeta=colData(se_blocks)
datMeta$Age = datMeta$Age + 0.00000001
dds=DESeqDataSetFromMatrix(countData = round(assays(se_blocks)$counts),
                           colData = datMeta,
                           design = ~Diagnosis + Sex + Age + seqPC1 + seqPC2 + seqPC3)
dds <- DESeq(dds)
resultsNames(dds)
dge.deseq2 <- results(dds, name="Diagnosis_Autism_vs_Control")
save(file = "./working_data/dge/dge.blocks.deseq2.v2.moreCovariates.RData",dge.deseq2)
table(dge.deseq2$padj < 0.05)

# neuron
datMeta=colData(se_neuron)
datMeta$Age = datMeta$Age + 0.00000001
dds=DESeqDataSetFromMatrix(countData = round(assays(se_neuron)$counts),
                           colData = datMeta,
                           design = ~Diagnosis + Sex + Age + seqPC2 + seqPC3)
dds <- DESeq(dds)
resultsNames(dds)
dge.deseq2 <- results(dds, name="Diagnosis_Autism_vs_Control")
save(file = "./working_data/dge/dge.neuron.deseq2.v2.moreCovariates.RData",dge.deseq2)
table(dge.deseq2$padj < 0.05)
dat=data.frame(dge.deseq2)
annot=read.delim("D:/references/gencode.v19.gene.name.txt")
dat$genenames=annot$genename[match(rownames(dat),annot$geneid)]

# oligo
datMeta=colData(se_oligo)
datMeta$Age = datMeta$Age + 0.00000001
dds=DESeqDataSetFromMatrix(countData = round(assays(se_oligo)$counts),
                           colData = datMeta,
                           design = ~Diagnosis + Sex + Age + seqPC1 + seqPC2)
dds <- DESeq(dds)
resultsNames(dds)
dge.deseq2 <- results(dds, name="Diagnosis_Autism_vs_Control")
save(file = "./working_data/dge/dge.oligo.deseq2.v2.moreCovariates.RData",dge.deseq2)
table(dge.deseq2$padj < 0.05)
