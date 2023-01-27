rm(list=ls());gc()
options(stringsAsFactors = F)

library(DESeq2)
library(ggplot2)
library(SummarizedExperiment)


load("working_data/summarizedExperiment/fc_2pass//se_blocks_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/fc_2pass//se_neuron_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/fc//se_oligo_CPM_outlierRemoved.RData")

se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
se_oligo$Diagnosis = factor(se_oligo$Diagnosis, levels= c("Control", "Autism"))
se_neuron$Diagnosis = factor(se_neuron$Diagnosis, levels= c("Control", "Autism"))

se_neuron=se_neuron[,!se_neuron$Sample %in% c("AO_B02","6221_N")]


# blocks
datMeta=colData(se_blocks)
datMeta$Age = datMeta$Age + 0.00000001
dds=DESeqDataSetFromMatrix(countData = round(assays(se_blocks)$counts),
                           colData = datMeta,
                           design = ~Diagnosis + Sex + Age + RIN + X260.280 + seqPC1 + seqPC2 + seqPC3)
dds <- DESeq(dds)
resultsNames(dds)
dge.deseq2 <- results(dds, name="Diagnosis_Autism_vs_Control")
save(file = "./working_data/dge/dge.blocks.fc.addCov.RData",dge.deseq2)
table(dge.deseq2$padj < 0.05)
dat=data.frame(dge.deseq2)
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
annot=read.delim("D:/references/gencode.v19.gene.name.txt")
dat$genenames=annot$genename[match(rownames(dat),annot$geneid)]
write.table(dat,file = "working_data/dge/dge.blocks.fc.addCov.txt",quote = F,sep = "\t",row.names = T)


# neuron
datMeta=colData(se_neuron)
datMeta$Age = datMeta$Age + 0.00000001
dds=DESeqDataSetFromMatrix(countData = round(assays(se_neuron)$counts),
                           colData = datMeta,
                           design = ~Diagnosis + Sex + Age + Type_RNAseqRunNumber + seqPC1 + seqPC2 + seqPC3)
dds <- DESeq(dds)
resultsNames(dds)
dge.deseq2 <- results(dds, name="Diagnosis_Autism_vs_Control")
save(file = "./working_data/dge/dge.neuron.fc.manualRemove.RData",dge.deseq2)
table(dge.deseq2$padj < 0.05)
dat=data.frame(dge.deseq2)
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
annot=read.delim("D:/references/gencode.v19.gene.name.txt")
dat$genenames=annot$genename[match(rownames(dat),annot$geneid)]
write.table(dat,file = "working_data/dge/dge.neuron.fc.manualRemove.txt",quote = F,sep = "\t",row.names = T)

# oligo
datMeta=colData(se_oligo)
datMeta$Age = datMeta$Age + 0.00000001
dds=DESeqDataSetFromMatrix(countData = round(assays(se_oligo)$counts),
                           colData = datMeta,
                           design = ~Diagnosis + Sex + Age + seqPC1 + seqPC2 + seqPC3)
dds <- DESeq(dds)
resultsNames(dds)
dge.deseq2 <- results(dds, name="Diagnosis_Autism_vs_Control")
save(file = "./working_data/dge/dge.oligo.fc_multimap.addCov.RData",dge.deseq2)
table(dge.deseq2$padj < 0.05)
dat=data.frame(dge.deseq2)
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
annot=read.delim("D:/references/gencode.v19.gene.name.txt")
dat$genenames=annot$genename[match(rownames(dat),annot$geneid)]
write.table(dat,file = "working_data/dge/dge.oligo.fc_multimap.addCov.txt",quote = F,sep = "\t",row.names = T)





# combine

rm(list=ls());gc()
options(stringsAsFactors = F)

library(DESeq2)
library(ggplot2)
library(SummarizedExperiment)


load("working_data/summarizedExperiment/fc_2pass/se_combine_CPM_outlierRemoved.RData")


se_combine$Diagnosis = factor(se_combine$Diagnosis, levels=c("Control", "Autism"))
se_combine$Type = factor(se_combine$Type, levels=c("Block", "Neuron"))

# blocks
datMeta=colData(se_combine)
datMeta$Age = datMeta$Age + 0.00000001
dds=DESeqDataSetFromMatrix(countData = round(assays(se_combine)$counts),
                           colData = datMeta,
                           design = ~Diagnosis + Sex + Age + Type + RNAseqPool + seqPC1 + seqPC2 + seqPC3 + Type:Diagnosis)
dds <- DESeq(dds)
resultsNames(dds)
block.ASDvCTL=results(dds,contrast = c("Diagnosis","Autism","Control"))
neuron.ASDvCTL=results(dds,contrast = list(c("Diagnosis_Autism_vs_Control","DiagnosisAutism.TypeNeuron")))
neuron.modify=results(dds,name = "DiagnosisAutism.TypeNeuron")
save(file = "./working_data/dge/dge.combine.fc2pass.RData",block.ASDvCTL,neuron.ASDvCTL,neuron.modify)
table(block.ASDvCTL$padj < 0.05)
dat=data.frame(neuron.ASDvCTL)
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")
dat=merge(dat,annot,by.x="row.names",by.y="geneid")

block_result=dat
neuron_result=dat
write.table(dat,file = "working_data/dge/dge.blocks.fc.addCov.txt",quote = F,sep = "\t",row.names = T)

