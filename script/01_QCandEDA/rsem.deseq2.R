rm(list=ls());gc()
options(stringsAsFactors = F)

library(DESeq2)
library(ggplot2)
library(SummarizedExperiment)

load("working_data/summarizedExperiment/se_blocks.RData")
se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
datMeta=as.data.frame(colData(se_blocks))
datMeta=datMeta[,-c(16:82)]
datMeta$Age = datMeta$Age + 0.00000001
counts=read.delim("working_data/star_rsem/genes.results.txt",row.names = 1)
counts1=counts[rowSums(counts>0) > 0.8*ncol(counts),]
table(datMeta$Sample %in% colnames(counts1))
counts2=counts1[,match(datMeta$Sample,colnames(counts1))]
dds=DESeqDataSetFromMatrix(countData = round(assays(se_blocks)$counts),
                           colData = datMeta,
                           design = ~Diagnosis + Sex + Age + X260.280 + seqPC1 + seqPC2 + seqPC3 +seqPC4)
dds <- DESeq(dds)


resultsNames(dds)
dge.deseq2 <- results(dds, name="Diagnosis_Autism_vs_Control")
save(file = "./working_data/dge/dge.blocks.deseq2.rsem.RData",dge.deseq2)
table(dge.deseq2$padj < 0.05)


datMeta$Sample[!datMeta$Sample %in% colnames(counts1)]
