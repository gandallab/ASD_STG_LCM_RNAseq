options(stringsAsFactors = F)

library(dplyr)
library(readr)
library(rtracklayer)
library(pheatmap)
library(SummarizedExperiment)

# Load Gencode v25 Annotations
if(!file.exists( './raw_data/gencode.v25lift37.annotation.gtf.gz')) download.file(url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz', destfile =  './raw_data/gencode.v25lift37.annotation.gtf.gz')
annot = import(gzfile( './raw_data/gencode.v25lift37.annotation.gtf.gz'))
annot = as.data.frame(annot[annot$type=='gene' & !grepl('PAR_Y', annot$gene_id),])
rownames(annot) = substr(annot$gene_id,0,15)
annot = annot[order(rownames(annot)),]


# Load expression data, match and order gene names
datExpr_cells = as.data.frame(read_csv("raw_data/salmon/salmon_counts_ASD_Cells.csv"))
rownames(datExpr_cells) = datExpr_cells$X1
datExpr_cells = datExpr_cells[,-1]

datExpr_blocks = as.data.frame(read_csv("raw_data/salmon/salmon_counts_ASD_Blocks.csv"))
datExpr_blocks = datExpr_blocks[!grepl('PAR_Y', datExpr_blocks$X1),]
rownames(datExpr_blocks) = substr(datExpr_blocks$X1,0,15)
datExpr_blocks = datExpr_blocks[,-1]
datExpr_blocks = datExpr_blocks[order(rownames(datExpr_blocks)),]

# Load Meta data
datMeta = read.csv("./raw_data/metaData/Metadata_ASD_CellsandBlocks.csv")
rownames(datMeta)= datMeta$Sample
datMeta_blocks = datMeta[match(colnames(datExpr_blocks), datMeta$Sample),]
datMeta_cells = datMeta[match(colnames(datExpr_cells), datMeta$Sample),]


# Split cells into Neurons and Oligos
datMeta_neuron = datMeta_cells %>% filter(Type=="Neuron")
datExpr_neuron = datExpr_cells[,match(datMeta_neuron$Sample, colnames(datExpr_cells))]
datMeta_oligo = datMeta_cells %>% filter(Type=="Oligo")
datExpr_oligo = datExpr_cells[,match(datMeta_oligo$Sample, colnames(datExpr_cells))]


# Load Seq Data for Blocks
datSeq_blocks = read.delim("./raw_data/QC/ASD_Blocks/multiqc_general_stats-salmon.txt")
datSeq_blocks$Sample = gsub("_abundance", "", datSeq_blocks$Sample)
datSeq_blocks = merge(datSeq_blocks, read.delim("./raw_data/QC/ASD_Blocks/multiqc_general_stats.txt"), by='Sample')
datSeq_blocks = merge(datSeq_blocks, read.delim("./raw_data/QC/ASD_Blocks/multiqc_picard_AlignmentSummaryMetrics.txt"), by='Sample')
datSeq_blocks = merge(datSeq_blocks, read.delim("./raw_data/QC/ASD_Blocks/multiqc_picard_dups.txt"), by='Sample')
datSeq_blocks = merge(datSeq_blocks, read.delim("./raw_data/QC/ASD_Blocks/multiqc_picard_gcbias.txt"), by='Sample')
datSeq_blocks = merge(datSeq_blocks, read.delim("./raw_data/QC/ASD_Blocks/multiqc_picard_RnaSeqMetrics.txt"), by='Sample')
rownames(datSeq_blocks) = datSeq_blocks$Sample
datSeq_blocks = datSeq_blocks[match(rownames(datMeta_blocks), rownames(datSeq_blocks)),]

# Remove columns with no variance or non-numeric data
v = apply(datSeq_blocks, 2, var)
datSeq_blocks = datSeq_blocks[,!is.na(v) & v>0]

library(pheatmap)
pheatmap(cor(datSeq_blocks, method='spearman'),fontsize = 6)

pcs = prcomp(datSeq_blocks,center = T, scale. = T)
screeplot(pcs,type = 'lines')
colnames(pcs$x) = paste0("seq", colnames(pcs$x))
datMeta_blocks = cbind(datMeta_blocks, pcs$x[,1:10])
datMeta_blocks = cbind(datMeta_blocks, datSeq_blocks)
rm(datSeq_blocks, pcs)



# Load Seq Data for Cells
datSeq_cells = read.delim("./raw_data/QC/ASD_Cells//multiqc_general_stats_salmon.txt")
datSeq_cells$Sample = gsub("_abundance", "", datSeq_cells$Sample)
datSeq_cells = merge(datSeq_cells, read.delim("./raw_data/QC/ASD_Cells/multiqc_general_stats.txt"), by='Sample')
datSeq_cells = merge(datSeq_cells, read.delim("./raw_data/QC/ASD_Cells/multiqc_picard_AlignmentSummaryMetrics.txt"), by='Sample')
datSeq_cells = merge(datSeq_cells, read.delim("./raw_data/QC/ASD_Cells/multiqc_picard_dups.txt"), by='Sample')
datSeq_cells = merge(datSeq_cells, read.delim("./raw_data/QC/ASD_Cells/multiqc_picard_gcbias.txt"), by='Sample')
datSeq_cells = merge(datSeq_cells, read.delim("./raw_data/QC/ASD_Cells/multiqc_picard_RnaSeqMetrics.txt"), by='Sample')
rownames(datSeq_cells) = datSeq_cells$Sample
datSeq_cells = datSeq_cells[match(rownames(datMeta_cells), rownames(datSeq_cells)),]


datSeq_neuron = datSeq_cells[match(datMeta_neuron$Sample, rownames(datSeq_cells)),]
datSeq_oligo = datSeq_cells[match(datMeta_oligo$Sample, rownames(datSeq_cells)),]

# Neuron seqPCs
# Remove columns with no variance or non-numeric data
v = apply(datSeq_neuron, 2, var)
datSeq_neuron = datSeq_neuron[,!is.na(v) & v>0]
pheatmap(cor(datSeq_neuron, method='spearman'),fontsize = 6)
pcs = prcomp(datSeq_neuron,center = T, scale. = T)
screeplot(pcs,type = 'lines')
colnames(pcs$x) = paste0("seq", colnames(pcs$x))
datMeta_neuron = cbind(datMeta_neuron, pcs$x[,1:10])
datMeta_neuron = cbind(datMeta_neuron, datSeq_neuron)
rm(datSeq_neuron, pcs)

# Oligo seqPCs
v = apply(datSeq_oligo, 2, var)
datSeq_oligo = datSeq_oligo[,!is.na(v) & v>0]
pheatmap(cor(datSeq_oligo, method='spearman'),fontsize = 6)
pcs = prcomp(datSeq_oligo,center = T, scale. = T)
screeplot(pcs,type = 'lines')
colnames(pcs$x) = paste0("seq", colnames(pcs$x))
datMeta_oligo = cbind(datMeta_oligo, pcs$x[,1:10])
datMeta_oligo = cbind(datMeta_oligo, datSeq_oligo)
rm(datSeq_oligo, pcs)
rm(datSeq_cells)




se_neuron = SummarizedExperiment(assays=list(counts=as.matrix(datExpr_neuron)),
                                rowData = annot, colData = datMeta_neuron)
save(file='./working_data/summarizedExperiment/se_neuron.RData',se_neuron)

se_oligo = SummarizedExperiment(assays=list(counts=as.matrix(datExpr_oligo)),
                                 rowData = annot, colData = datMeta_oligo)
save(file='./working_data/summarizedExperiment/se_oligo.RData',se_oligo)

se_blocks = SummarizedExperiment(assays=list(counts=as.matrix(datExpr_blocks)),
                                rowData = annot, colData = datMeta_blocks)
save(file='./working_data/summarizedExperiment/se_blocks.RData',se_blocks)

