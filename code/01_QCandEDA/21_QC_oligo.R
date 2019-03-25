#21_QC_oligo.R

options(stringsAsFactors = F)

library(readr)
library(rtracklayer)
library(pheatmap)
library(SummarizedExperiment)
library(edgeR)
library(WGCNA)

load("working_data/summarizedExperiment/se_oligo.RData")

se_oligo = se_oligo[,colnames(se_oligo) != "6469Oligo"]

##### (1) Visualization on raw data  #####
boxplot(assays(se_oligo)$counts, range=0, main="Raw Counts")
boxplot(log2(1+assays(se_oligo)$counts), range=0, main = "log2 counts + 1")


i=1;  plot(density(log2(.001+assays(se_oligo)$counts[,i])),col = as.factor(se_oligo$Diagnosis)[i], main = "Counts from tximport, countsFromAbundance = 'lengthScaledTPM'", xlab="log2(raw_counts + 1)", xlim=c(-15,30), ylim=c(0,0.45))
for(i in 2:ncol(se_oligo)) {
  lines(density(log2(.001+assays(se_oligo)$counts[,i])),col = as.factor(se_oligo$Diagnosis)[i])
}
legend("topright", legend = unique(se_oligo$Diagnosis), fill = 1:3)


##### (2) Remove low expressed genes ####
genes_to_keep = filterByExpr(calcNormFactors(DGEList(assays(se_oligo)$counts)))
table(genes_to_keep)
se_oligo = se_oligo[genes_to_keep,]


##### (3) Remove outlier samples ####
#calculate top 10 principle components
mds = cmdscale(dist(t(cpm(calcNormFactors(DGEList(assays(se_oligo)$counts), method = 'TMM'), log = TRUE))),   k = 10)
colnames(mds) = paste0("PC",1:ncol(mds))

pairs(mds, col=factor(se_oligo$Diagnosis), main="Group", pch=19)
par(xpd = TRUE,oma=c(1,1,1,1)); legend('bottomright', levels(factor(se_oligo$Diagnosis)),fill=1:2,cex=.5)

Zscore = apply(as.matrix(mds),2, scale)
rownames(Zscore) = rownames(mds)
Zscore = (abs(Zscore) > 3)
table(rowSums((Zscore)))
to_remove = rowSums(Zscore) >= 2


normadj <- 0.5 + 0.5*bicor(cpm(calcNormFactors(DGEList(assays(se_oligo)$counts), method = 'TMM'), log = TRUE))
netsummary <- fundamentalNetworkConcepts(normadj)
C <- netsummary$Connectivity
Z.C <- (C-mean(C))/sqrt(var(C))

# plot the connectivity z-score for the normalized data and set a cut-off at 3 STD
plot(1:length(Z.C),Z.C,main="Outlier Plot of log2(normalized_counts)",xlab = "Samples",ylab="Connectivity Z Score", col=factor(se_oligo$Diagnosis))
text(1:length(Z.C),Z.C - .2, labels = names(Z.C))
abline(h= -2, col="red")

# determine which samples fail the 'outlier threshold' via connectivity z-score
to_remove = to_remove | Z.C < -2
to_remove = to_remove | se_oligo$Diagnosis==""
table(to_remove)

se_oligo = se_oligo[, !to_remove]


##### (4) Calculate normalized gene expression on outlier removed dataset
assays(se_oligo)$log2CPM = cpm(calcNormFactors(DGEList(assays(se_oligo)$counts), method = 'TMM'), log = TRUE)

boxplot(assays(se_oligo)$log2CPM,range=0, col=numbers2colors(se_oligo$TOTAL_READS))
boxplot(assays(se_oligo)$log2CPM,range=0, col=as.factor(se_oligo$Diagnosis))

i=1;  plot(density(assays(se_oligo)$log2CPM[,i]),col = as.factor(se_oligo$Diagnosis)[i], main = "Counts from tximport, countsFromAbundance = 'lengthScaledTPM'", xlab="log2(raw_counts + 1)", ylim=c(0,.5))
for(i in 2:ncol(assays(se_oligo)$log2CPM)){
  lines(density(assays(se_oligo)$log2CPM[,i]),col = as.factor(se_oligo$Diagnosis)[i])
}
legend("topright", legend = unique(se_oligo$Diagnosis), fill = 1:2)

#calculate top 10 principle components
mds = cmdscale(dist(t(assays(se_oligo)$log2CPM)),k = 10)
colnames(mds) = paste0("PC",1:ncol(mds))


pairs(mds, col=factor(se_oligo$Diagnosis), main="Group", pch=19)
pairs(mds, col=numbers2colors(se_oligo$Age), main="Age", pch=19)
pairs(mds, col=numbers2colors(se_oligo$PMI), main="PMI", pch=19)
pairs(mds, col=numbers2colors(se_oligo$seqPC1), main="seqPC1", pch=19)
pairs(mds, col=numbers2colors(se_oligo$seqPC3), main="seqPC2", pch=19)
pairs(mds, col=factor(se_oligo$Sex), main="Sex", pch=19)
pairs(mds, col=factor(se_oligo$RNAseqPool), main="RNAseqPool", pch=19)


save(file="working_data/summarizedExperiment/se_oligo_CPM_outlierRemoved.RData", se_oligo)
