rm(list=ls());gc()
options(stringsAsFactors = F)

library(VennDiagram)
library(export)

load("./working_data/dge/dge.neuron.deseq2.RData")
load("./working_data/dge/dge.neuron.edger.RData")

ds2=as.data.frame(dge.deseq2)
er=dge.edger$table

ds2.genes=rownames(ds2)[ds2$padj < 0.05 & !is.na(ds2$padj)]
er.genes=rownames(er)[er$FDR < 0.05]

draw.pairwise.venn(area1 = length(ds2.genes),
                   area2 = length(er.genes),
                   cross.area = length(intersect(ds2.genes,er.genes)),
                    category = c("DESeq2","EdgeR"),
                   fill = c("green","red"))

graph2ppt(file="./working_data/figures/deseq2.edger.venn.pptx")

ds2.er.overlap.fdr=ds2$padj[rownames(ds2) %in% intersect(er.genes,ds2.genes)]
ds2.only.fdr=ds2$padj[rownames(ds2) %in% setdiff(ds2.genes,er.genes)]
boxplot(ds2.er.overlap.fdr,ds2.only.fdr,
        names=c("DESeq2.EdgeR.Overlap","DESeq2.Only"),
        ylab="FDR Distribution")
graph2ppt(file="./working_data/figures/deseq2.edger.FDR.distributino.pptx")
