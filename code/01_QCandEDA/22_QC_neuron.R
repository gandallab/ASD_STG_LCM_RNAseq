#21_QC_neuron.R

options(stringsAsFactors = F)

library(readr)
library(rtracklayer)
library(pheatmap)
library(SummarizedExperiment)
library(edgeR)
library(WGCNA)
library(DESeq2)
# 
# library("BiocParallel")
# register(MulticoreParam(4))

load("working_data/summarizedExperiment/se_neuron.RData")


##### (1) Visualization on raw data  #####
boxplot(assays(se_neuron)$counts, range=0, main="Raw Counts")
boxplot(log2(1+assays(se_neuron)$counts), range=0, main = "log2 counts + 1")


i=1;  plot(density(log2(.001+assays(se_neuron)$counts[,i])),col = as.factor(se_neuron$Diagnosis)[i], main = "Counts from tximport, countsFromAbundance = 'lengthScaledTPM'", xlab="log2(raw_counts + 1)", xlim=c(-15,30), ylim=c(0,0.45))
for(i in 2:ncol(se_neuron)) {
  lines(density(log2(.001+assays(se_neuron)$counts[,i])),col = as.factor(se_neuron$Diagnosis)[i])
}
legend("topright", legend = unique(se_neuron$Diagnosis), fill = 1:3)


##### (2) Remove low expressed genes ####
genes_to_keep = filterByExpr(DGEList(assays(se_neuron)$counts))
table(genes_to_keep)
se_neuron = se_neuron[genes_to_keep,]


##### (3) Remove outlier samples ####
#calculate top 10 principle components
mds = cmdscale(dist(t(cpm(calcNormFactors(DGEList(assays(se_neuron)$counts), method = 'TMM'), log = TRUE))),   k = 5)
colnames(mds) = paste0("PC",1:ncol(mds))

pairs(mds, col=factor(se_neuron$Diagnosis), main="Group", pch=19)
par(xpd = TRUE,oma=c(1,1,1,1)); legend('bottomright', levels(factor(se_neuron$Diagnosis)),fill=1:2,cex=.5)

Zscore = apply(as.matrix(mds),2, scale)
rownames(Zscore) = rownames(mds)

Zscore = (abs(Zscore) > 3)
table(rowSums((Zscore)))
to_remove = rowSums(Zscore) >= 1
which(to_remove)

normadj <- 0.5 + 0.5*bicor(cpm(calcNormFactors(DGEList(assays(se_neuron)$counts), method = 'TMM'), log = TRUE))
netsummary <- fundamentalNetworkConcepts(normadj)
C <- netsummary$Connectivity
Z.C <- (C-mean(C))/sqrt(var(C))

# plot the connectivity z-score for the normalized data and set a cut-off at 3 STD
plot(1:length(Z.C),Z.C,main="Outlier Plot of log2(normalized_counts)",xlab = "Samples",ylab="Connectivity Z Score", col=factor(to_remove))
text(1:length(Z.C),Z.C - .2, labels = names(Z.C))
abline(h= -2, col="red")

# determine which samples fail the 'outlier threshold' via connectivity z-score
#to_remove = to_remove | Z.C < -2
to_remove = to_remove | se_neuron$Diagnosis==""
table(to_remove)

se_neuron = se_neuron[, !to_remove]


##### (4) Calculate normalized gene expression on outlier removed dataset
assays(se_neuron)$log2CPM = cpm(calcNormFactors(DGEList(assays(se_neuron)$counts), method = 'TMM'), log = TRUE)

boxplot(assays(se_neuron)$log2CPM,range=0, col=numbers2colors(se_neuron$TOTAL_READS))
boxplot(assays(se_neuron)$log2CPM,range=0, col=as.factor(se_neuron$Diagnosis))

i=1;  plot(density(assays(se_neuron)$log2CPM[,i]),col = as.factor(se_neuron$Diagnosis)[i], main = "Counts from tximport, countsFromAbundance = 'lengthScaledTPM'", xlab="log2(raw_counts + 1)", ylim=c(0,.5))
for(i in 2:ncol(assays(se_neuron)$log2CPM)){
  lines(density(assays(se_neuron)$log2CPM[,i]),col = as.factor(se_neuron$Diagnosis)[i])
}
legend("topright", legend = unique(se_neuron$Diagnosis), fill = 1:2)

#calculate top 10 principle components
mds = cmdscale(dist(t(assays(se_neuron)$log2CPM)),k = 10)
colnames(mds) = paste0("PC",1:ncol(mds))


pairs(mds, col=factor(se_neuron$Diagnosis), main="Group", pch=19)
pairs(mds, col=numbers2colors(se_neuron$Age), main="Age", pch=19)
pairs(mds, col=numbers2colors(se_neuron$PMI), main="PMI", pch=19)
pairs(mds, col=numbers2colors(se_neuron$seqPC1), main="seqPC1", pch=19)
pairs(mds, col=numbers2colors(se_neuron$seqPC3), main="seqPC2", pch=19)
pairs(mds, col=factor(se_neuron$Sex), main="Sex", pch=19)
pairs(mds, col=factor(se_neuron$RNAseqPool), main="RNAseqPool", pch=19)

save(file="working_data/summarizedExperiment/se_neuron_CPM_outlierRemoved.RData", se_neuron)

se_neuron$Diagnosis =factor(se_neuron$Diagnosis, levels=c("Control", "Autism"))


dds = DESeqDataSetFromMatrix(countData = round(assays(se_neuron)$counts),
                             colData = colData(se_neuron),
                             design = ~Diagnosis + Sex + Age + seqPC1 + seqPC2 + seqPC3)
dds = DESeq(dds)
res = results(dds, name='Diagnosis_Control_vs_Autism', filterFun = ihw)
res$gene = rowData(se_neuron)$gene_name
res$gene_type = rowData(se_neuron)$gene_type

#res <- lfcShrink(dds, coef="Diagnosis_Control_vs_Autism", type="apeglm")
sum(res$padj < 0.1, na.rm=TRUE)
res = as.data.frame(res)
res$log2FoldChange = -1*res$log2FoldChange

ggplot(as.data.frame(res),aes(x=log2FoldChange,y=-log10(pvalue), color=gene_type=="snoRNA")) + geom_point() + theme_bw() + 
  ggrepel::geom_text_repel(res[res$padj<.1 & res$gene_type=="snoRNA",],mapping=aes(x=log2FoldChange,y=-log10(pvalue), label=gene))
