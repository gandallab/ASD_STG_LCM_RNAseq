rm(list=ls());gc()
options(stringsAsFactors = F)

library(limma)
library(edgeR)

load("./3months_Organoids_RSEM_Quant.genes.tpm.RData")

tpm1=tpm[rowSums(tpm>0.5) > 0.8*ncol(tpm),]

datMeta=read.delim("./3m_Organoids_metadata_pc_SVA.txt")
load("./3months_Organoids_RSEM_Quant.genes.counts.RData")
counts5=counts[rownames(tpm1),rownames(datMeta)]
dge.voom = voom(calcNormFactors(DGEList(counts = counts5),method = 'TMM'))

datExpr = as.matrix(dge.voom$E)
save(datMeta,datExpr,file = "./forContext/3mOrganoid.forContext.exprGenes.RData")


# use full dataset

rm(list=ls());gc()
options(stringsAsFactors = F)

library(limma)
library(edgeR)

datMeta=read.delim("./3m_metadata_pc_SVA_nohvg.txt")
gene=read.delim("3mOrg_Expressed_genes_tpm0.5_byGenotype.plus16p.txt",header = F)
gene.remove=read.delim("3mOrg_highly_variable_genes_betweenLabs.txt",header = F)
gene=gene[!gene$V1 %in% gene.remove$V1,]
load("3months_Organoids_RSEM_Quant.genes.counts.RData")
counts5=counts[gene,datMeta$Sample]
dge.voom = voom(calcNormFactors(DGEList(counts = counts5),method = 'TMM'))

datExpr = as.matrix(dge.voom$E)
table(colnames(datExpr)==datMeta$X)
save(datMeta,datExpr,file = "./forContext/3mOrganoid.forContext.exprGenes.fulldata.RData")
