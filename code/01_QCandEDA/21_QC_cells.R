#21_QC_cells.R

options(stringsAsFactors = F)

library(readr)
library(rtracklayer)
library(pheatmap)
library(SummarizedExperiment)
library(edgeR)
library(WGCNA)

load("working_data/summarizedExperiment/se_cells.RData")


se_cells = se_cells[,colnames(se_cells) != "6469Oligo" & se_cells$Type != ""]

##### (1) Visualization on raw data  #####
boxplot(assays(se_cells)$counts, range=0, main="Raw Counts")
boxplot(log2(1+assays(se_cells)$counts), range=0, main = "log2 counts + 1")


i=1;  plot(density(log2(.001+assays(se_cells)$counts[,i])),col = as.factor(se_cells$Diagnosis)[i], main = "Counts from tximport, countsFromAbundance = 'lengthScaledTPM'", xlab="log2(raw_counts + 1)", xlim=c(-15,30), ylim=c(0,0.45))
for(i in 2:ncol(se_cells)) {
  lines(density(log2(.001+assays(se_cells)$counts[,i])),col = as.factor(se_cells$Diagnosis)[i])
}
legend("topright", legend = unique(se_cells$Diagnosis), fill = 1:3)


#calculate top 10 principle components
mds = cmdscale(dist(t(log2(assays(se_cells)$counts + .1))),k = 10)
colnames(mds) = paste0("PC",1:ncol(mds))

pairs(mds, col=factor(se_cells$Diagnosis), main="Group", pch=19)
par(xpd = TRUE,oma=c(1,1,1,1)); legend('bottomright', levels(factor(se_cells$Diagnosis)),fill=1:2,cex=.5)

pairs(mds, col=factor(se_cells$Type), main="CellType", pch=19)
par(xpd = TRUE,oma=c(1,1,1,1)); legend('bottomright', levels(factor(se_cells$Type)),fill=1:2,cex=.5)




##### (2) Remove low expressed genes and outlier samples ####
genes_to_keep = filterByExpr(DGEList(assays(se_cells)$counts))
table(genes_to_keep)


normadj <- 0.5 + 0.5*bicor(log2(.1+assays(se_cells)$counts))
netsummary <- fundamentalNetworkConcepts(normadj)
C <- netsummary$Connectivity
Z.C <- (C-mean(C))/sqrt(var(C))

# plot the connectivity z-score for the normalized data and set a cut-off at 3 STD
plot(1:length(Z.C),Z.C,main="Outlier Plot of log2(normalized_counts)",xlab = "Samples",ylab="Connectivity Z Score", col=factor(se_cells$Diagnosis))
abline(h= -2, col="red")

# determine which samples fail the 'outlier threshold' via connectivity z-score
to_remove = Z.C < -2
to_remove = to_remove | se_cells$Diagnosis==""
table(to_remove)

se_cells = se_cells[genes_to_keep, !to_remove]


##### (3) Calculate normalized gene expression on outlier removed dataset
assays(se_cells)$log2CPM = cpm(calcNormFactors(DGEList(assays(se_cells)$counts), method = 'TMM'), log = TRUE)

boxplot(assays(se_cells)$log2CPM,range=0, col=numbers2colors(se_cells$TOTAL_READS))

i=1;  plot(density(assays(se_cells)$log2CPM[,i]),col = as.factor(se_cells$Diagnosis)[i], main = "Counts from tximport, countsFromAbundance = 'lengthScaledTPM'", xlab="log2(raw_counts + 1)", ylim=c(0,.5))
for(i in 2:ncol(assays(se_cells)$log2CPM)){
  lines(density(assays(se_cells)$log2CPM[,i]),col = as.factor(se_cells$Diagnosis)[i])
}
legend("topright", legend = unique(se_cells$Diagnosis), fill = 1:2)

#calculate top 10 principle components
mds = cmdscale(dist(t(assays(se_cells)$log2CPM)),k = 10)
colnames(mds) = paste0("PC",1:ncol(mds))

pairs(mds, col=factor(se_cells$Diagnosis), main="Group", pch=19)
pairs(mds, col=numbers2colors(se_cells$Age), main="Age", pch=19)
pairs(mds, col=numbers2colors(se_cells$PMI), main="PMI", pch=19)
pairs(mds, col=numbers2colors(se_cells$seqPC1), main="seqPC1", pch=19)
pairs(mds, col=numbers2colors(se_cells$seqPC3), main="seqPC2", pch=19)
pairs(mds, col=factor(se_cells$Sex), main="Sex", pch=19)
pairs(mds, col=factor(se_cells$RNAseqPool), main="RNAseqPool", pch=19)


save(file="working_data/summarizedExperiment/se_cells_CPM_outlierRemoved.RData", se_cells)









