rm(list=ls());gc()
options(stringsAsFactors = F)

library(DESeq2)
library(tidyverse)
library(SummarizedExperiment)


load("working_data/summarizedExperiment/fc_2pass//se_blocks_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/fc_2pass//se_neuron_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/fc//se_oligo_CPM_outlierRemoved.RData")

se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
se_oligo$Diagnosis = factor(se_oligo$Diagnosis, levels= c("Control", "Autism"))
se_neuron$Diagnosis = factor(se_neuron$Diagnosis, levels= c("Control", "Autism"))

#table(se_neuron$Sample %in% c("AO_B02","6221_N"))
#se_neuron=se_neuron[,!se_neuron$Sample %in% c("AO_B02","6221_N")]


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
save(file = "./working_data/dge/dge.blocks.fc2pass.interaction.RData",dge.deseq2)
table(dge.deseq2$padj < 0.05)
dat=data.frame(dge.deseq2)
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")
dat$genenames=annot$genename[match(rownames(dat),annot$geneid)]
write.table(dat,file = "working_data/dge/dge.blocks.fc2pass.interaction.txt",quote = F,sep = "\t",row.names = T)

load("working_data/wgcna/voom.forWGCNA.input.blocks.fc.RData")
datMeta$SPEN=datExpr["ENSG00000065526.10_2",]
datMeta$CLUL1=datExpr["ENSG00000079101.16_3",]
datMeta$GAD2=datExpr["ENSG00000136750.12_4",]
datMeta$GAD1=datExpr["ENSG00000128683.13_3",]
ggplot(datMeta,aes(Age,GAD1,color=Diagnosis,group=Diagnosis))+
  geom_point(size=3)+
  stat_summary(fun = mean,geom = "line",size=1)
ggplot(datMeta,aes(Age,GAD2,color=Diagnosis,group=Diagnosis))+
  geom_point(size=3)+
  geom_smooth(method = "lm",alpha=0.3)+
#  xlab("GOT1 Intron chr10:101163631-101165513 PSI")+
  ylab("GAD2 expression\nin block tissue (log2CPM)")+
  theme_bw()+
  theme(text = element_text(family = "Arial",face = "bold"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=16))
tmp1=read.delim("working_data/dge/dge.blocks.fc2pass.txt")
tmp2=read.delim("working_data/dge/dge.blocks.fc2pass.interaction.txt")
dat=tmp2
rownames(dds)[grep("ENSG00000128683",rownames(dds))]
sigdat=dat[dat$padj < 0.05,]
sigdat=sigdat[complete.cases(sigdat),]
sigdat=sigdat[order(sigdat$log2FoldChange),]
godown = gprofiler2::gost(query=rownames(sigdat)[sigdat$log2FoldChange < 0], 
                 correction_method = "fdr",
                 exclude_iea = T, 
                 custom_bg = rownames(dat), 
                 sources = c("GO","KEGG"),
                 ordered_query = T)
tmp3=godown$result

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
save(file = "./working_data/dge/dge.neuron.fc2pass.interaction.RData",dge.deseq2)
table(dge.deseq2$padj < 0.05)
dat=data.frame(dge.deseq2)
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")
dat$genenames=annot$genename[match(rownames(dat),annot$geneid)]
write.table(dat,file = "working_data/dge/dge.neuron.fc2pass.interaction.txt",quote = F,sep = "\t",row.names = T)
plotinter=plotCounts(dds, "ENSG00000146001.5_4",  # "ENSG00000196417.12_3" , "ENSG00000146001.5_4","ENSG00000115317.11_2"
                     intgroup = c("Diagnosis","Age"), returnData = TRUE)
ggplot(plotinter,
       aes(x = Age, y = count, color = Diagnosis, group = Diagnosis)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()

load("working_data/wgcna/voom.forWGCNA.input.neuron.fc2pass.RData")
datMeta$Diagnosis=factor(datMeta$Diagnosis,levels = c("Control","Autism"))
datMeta$HTRA2=datExpr["ENSG00000115317.11_2",]
datMeta$ZNF765=datExpr["ENSG00000196417.12_3",]
datMeta$PCDHB18P=datExpr["ENSG00000146001.5_4",]
datMeta$GAD2=datExpr["ENSG00000136750.12_4",]
datMeta$GAD1=datExpr["ENSG00000128683.13_3",]
ggplot(datMeta,aes(Age,GAD2,color=Diagnosis,group=Diagnosis))+
  geom_point(size=3)+
  stat_summary(fun = mean,geom = "line",size=1)

ggplot(datMeta,aes(Age,HTRA2,color=Diagnosis,group=Diagnosis))+
  geom_point(size=3)+
  geom_smooth(method = "lm",alpha=0.3)+
  #  xlab("GOT1 Intron chr10:101163631-101165513 PSI")+
  ylab("HTRA2 expression\nin neuron")+
  theme_bw()+
  theme(text = element_text(family = "Arial",face = "bold"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=16))

ggplot(datMeta,aes(Age,GAD2,color=Diagnosis,group=Diagnosis))+
  geom_point(size=3)+
  geom_smooth(method = "lm",alpha=0.3)+
  #  xlab("GOT1 Intron chr10:101163631-101165513 PSI")+
  ylab("GAD2 expression\nin neuron")+
  scale_color_manual(values=c("Control"="#00BFC4","Autism"="#F8766D"))+
  theme_bw()+
  theme(text = element_text(family = "Arial",face = "bold"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=11),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=16))
ggplot_build(p)
