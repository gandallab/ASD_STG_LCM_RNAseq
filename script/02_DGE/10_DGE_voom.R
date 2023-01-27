#10_DGE_voom.R
options(stringsAsFactors = F)
library(limma)
library(ggplot2)
library(SummarizedExperiment)


load("working_data/summarizedExperiment/se_blocks_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/se_neuron_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/se_oligo_CPM_outlierRemoved.RData")

str(se_blocks)
se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
se_oligo$Diagnosis = factor(se_oligo$Diagnosis, levels= c("Control", "Autism"))
se_neuron$Diagnosis = factor(se_neuron$Diagnosis, levels= c("Control", "Autism"))


mod.blocks = model.matrix(~Diagnosis + Sex + Age + RIN + seqPC1 + seqPC2 + seqPC3 + seqPC4,data = colData(se_blocks))
v = voom(calcNormFactors(DGEList(assays(se_blocks)$counts)),design = mod.blocks, plot = T)
lm = lmFit(v, mod.blocks)
lm = eBayes(lm)

tt.blocks = topTable(lm, coef=2, number = Inf, sort.by = 'none', genelist = rowData(se_blocks)$gene_name)
table(tt.blocks$adj.P.Val < 0.05)
hist(tt.blocks$P.Value)


ggplot(tt.blocks, aes(x=logFC, y=-log10(P.Value), label=ID)) + geom_point() + geom_text(data = tt.blocks[tt.blocks$P.Value < .001,])





mod.neuron = model.matrix(~Diagnosis + Age + Sex + seqPC1 + seqPC2 + seqPC3, data = colData(se_neuron))

v = voom(calcNormFactors(DGEList(assays(se_neuron)$counts)),design = mod.neuron, plot = T)
lm = lmFit(v, mod.neuron)
lm = eBayes(lm)
tt.neuron = topTable(lm, coef=2, number = Inf, sort.by = 'none', genelist = rowData(se_neuron)$gene_name)
table(tt.neuron$adj.P.Val < 0.05)
hist(tt.neuron$P.Value)



mod.oligo= model.matrix(~Diagnosis + Age + Sex  + seqPC1 + seqPC2, data = colData(se_oligo))

v = voom(calcNormFactors(DGEList(assays(se_oligo)$counts)),design = mod.oligo, plot = T)
lm = lmFit(v, mod.oligo)
lm = eBayes(lm)
tt.oligo = topTable(lm, coef=2, number = Inf, sort.by = 'none', genelist = rowData(se_oligo)$gene_name)
table(tt.oligo$adj.P.Val < 0.05)
hist(tt.oligo$P.Value)
