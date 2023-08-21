rm(list=ls());gc()
options(stringsAsFactors = F)

library(DESeq2)
library(tidyverse)
library(SummarizedExperiment)


load("se_blocks.RData")
load("se_neuron.RData")


se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))

se_neuron$Diagnosis = factor(se_neuron$Diagnosis, levels= c("Control", "Autism"))



# blocks
datMeta=colData(se_blocks)
datMeta$Age = datMeta$Age + 0.00000001
dds=DESeqDataSetFromMatrix(countData = round(assays(se_blocks)$counts),
                           colData = datMeta,
                           design = ~Diagnosis * Age + Sex + RIN + X260.280 + seqPC1 + seqPC2 + seqPC3)
dds <- DESeq(dds)
resultsNames(dds)
dge.deseq2 <- results(dds, name="DiagnosisAutism.Age", independentFiltering = FALSE)
dge.deseq2 <- results(dds, name="Diagnosis_Autism_vs_Control")
save(file = "dge.blocks.fc2pass.interaction.RData",dge.deseq2)
table(dge.deseq2$padj < 0.05)
dat=data.frame(dge.deseq2)
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
annot=read.delim("gencode.v29lift37.genename.genetype.txt")
dat$genenames=annot$genename[match(rownames(dat),annot$geneid)]
write.table(dat,file = "dge.blocks.fc2pass.interaction.txt",quote = F,sep = "\t",row.names = T)


# neuron
datMeta=colData(se_neuron)
datMeta$Age = datMeta$Age + 0.00000001
dds=DESeqDataSetFromMatrix(countData = round(assays(se_neuron)$counts),
                           colData = datMeta,
                           design = ~ Diagnosis * Age + Sex + Type_RNAseqRunNumber + seqPC1 + seqPC2 + seqPC3)
dds <- DESeq(dds)
resultsNames(dds)
dge.deseq2 <- results(dds, name="DiagnosisAutism.Age")
dge.deseq2 <- results(dds, name="Diagnosis_Autism_vs_Control" )
save(file = "dge.neuron.fc2pass.interaction.RData",dge.deseq2)
table(dge.deseq2$padj < 0.05)
dat=data.frame(dge.deseq2)
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
annot=read.delim("gencode.v29lift37.genename.genetype.txt")
dat$genenames=annot$genename[match(rownames(dat),annot$geneid)]
write.table(dat,file = "dge.neuron.fc2pass.interaction.txt",quote = F,sep = "\t",row.names = T)

