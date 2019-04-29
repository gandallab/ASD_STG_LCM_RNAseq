rm(list=ls());gc()
options(stringsAsFactors = F)

library(ggplot2)
library(ggrepel)
library(scales)
library(export)
library(xlsx)
library(gProfileR)

asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

annot=read.delim("D:/references/gencode.v19.gene.name.txt")


# blocks
load("working_data/dge/dge.blocks.deseq2.v2.moreCovariates.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]
sigdat=dat[abs(dat$log2FoldChange) >= 0.58 & dat$padj < 0.05,]
nosigdat=dat[!rownames(dat) %in% rownames(sigdat),]
sigdat$genename=annot$genename[match(rownames(sigdat),annot$geneid)]
sigdat$genename[(sigdat$log2FoldChange > 0 & sigdat$log2FoldChange <= 1.2) & sigdat$padj > 0.006]=""
ggplot()+geom_point(data=sigdat,aes(log2FoldChange,-log10(padj)),size=2,color="black")+
  geom_point(data = nosigdat,aes(log2FoldChange,-log10(padj)),size=2,alpha=1,col="gray")+
  geom_text_repel(data = sigdat,aes(log2FoldChange,-log10(padj),label=genename),min.segment.length = unit(0, "lines"),col="red")+
  geom_hline(yintercept = 1.3)+
  geom_vline(xintercept = c(-0.58,0.58))+
  #  scale_color_manual(values = colors)+
  scale_x_continuous(trans = "asinh")+
  #scale_y_continuous(trans = "log10")+
  #  geom_vline(xintercept = c(-1,1))+                                                      
  theme_classic()+
  theme(text = element_text(family = "Arial"))

graph2ppt(file="working_data/figures/blocks.dge.plots.pptx",width=7.6,height=7.1,append=TRUE)

outwb=createWorkbook()
sigdat=sigdat[order(sigdat$log2FoldChange,decreasing = T),]
goup = gprofiler(query=rownames(sigdat)[sigdat$log2FoldChange > 0], 
               max_set_size = 1000,
               min_isect_size = 2,
               correction_method = "fdr",
               exclude_iea = T,
               hier_filtering = "moderate", 
               custom_bg = rownames(dge.deseq2), 
               src_filter = c("GO","KEGG"),
               ordered_query = T)
goup=goup[order(goup$p.value),]
sheet=createSheet(outwb,sheetName = "blocks.up")
addDataFrame(goup,sheet,row.names = F)

sigdat=sigdat[order(sigdat$log2FoldChange,decreasing = F),]
godown = gprofiler(query=rownames(sigdat)[sigdat$log2FoldChange < 0], 
                 max_set_size = 1000,
                 min_isect_size = 2,
                 correction_method = "fdr",
                 exclude_iea = T,
                 hier_filtering = "moderate", 
                 custom_bg = rownames(dge.deseq2), 
                 src_filter = c("GO","KEGG"),
                 ordered_query = F)
godown=godown[order(godown$p.value),]
sheet=createSheet(outwb,sheetName = "blocks.down")
addDataFrame(godown,sheet,row.names = F)

saveWorkbook(outwb,file = "working_data/dge/dge.go.blocks.v2.moreCovariates.xlsx")


goup=read.xlsx(file = "working_data/dge/dge.go.blocks.v2.moreCovariates.xlsx",sheetName = "blocks.up")
godown=read.xlsx(file = "working_data/dge/dge.go.blocks.v2.moreCovariates.xlsx",sheetName = "blocks.down")
goup$FDR=-log10(goup$p.value);goup=goup[1:10,];goup$dirc="up"
godown$FDR=log10(godown$p.value);godown=godown[1:10,];godown$dirc="down"
go=rbind(goup,godown);go=go[order(go$FDR,decreasing = T),]

ggplot(go, aes(x=reorder(term.name, FDR), y=FDR,fill=dirc)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("up" = "gold","down"="blue"))+
  coord_flip() + 
  xlab("") + 
  ylab("Signed -log10(FDR)") +
  geom_hline(yintercept=c(log10(0.05),-log10(0.05)), lty=2, color="red")+
  theme(text = element_text(family = "Arial",size = 15,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20))+
  guides(fill=F)
graph2ppt(file="working_data/figures/blocks.dge.plots.pptx",width=10.5,height=7.1,append=TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)
dat1=as.data.frame(dge.deseq2)
geneList=dat1$log2FoldChange
names(geneList)=rownames(dat1)
geneList=sort(geneList,decreasing = T)

ego <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
             keyType = "ENSEMBL",
              ont          = "BP",
              nPerm        = 500,
              minGSSize    = 2,
              maxGSSize    = 1000,
              pvalueCutoff = 1,
              verbose      = TRUE)

ego.table=ego@result
ego.table=ego.table[complete.cases(ego.table),]
ego.table=ego.table[order(ego.table$NES,decreasing = T),]
ego.plot=rbind(head(ego.table,10),tail(ego.table,10))

ggplot(ego.plot, aes(x=reorder(Description, NES), y=NES,fill=p.adjust)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  xlab("") + 
  ylab("NES") +
  theme(text = element_text(family = "Arial",size = 15,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20))
graph2ppt(file="working_data/figures/blocks.dge.plots.pptx",width=10.5,height=7.1,append=TRUE)




# neuron
load("working_data/dge/dge.neuron.deseq2.v2.moreCovariates.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]
sigdat=dat[abs(dat$log2FoldChange) >= 0.58 & dat$padj < 0.05,]
nosigdat=dat[!rownames(dat) %in% rownames(sigdat),]
sigdat$genename=annot$genename[match(rownames(sigdat),annot$geneid)]
sigdat$genename[(sigdat$log2FoldChange >= -1.8 & sigdat$log2FoldChange <= 1.5) & sigdat$padj > 0.006]=""
ggplot()+geom_point(data=sigdat,aes(log2FoldChange,-log10(padj)),size=2,color="black")+
  geom_point(data = nosigdat,aes(log2FoldChange,-log10(padj)),size=2,alpha=1,col="gray")+
  geom_text_repel(data = sigdat,aes(log2FoldChange,-log10(padj),label=genename),min.segment.length = unit(0, "lines"),col="red")+
  geom_hline(yintercept = 1.3)+
  geom_vline(xintercept = c(-0.58,0.58))+
  #  scale_color_manual(values = colors)+
  scale_x_continuous(trans = "asinh")+
  #scale_y_continuous(trans = "log10")+
  #  geom_vline(xintercept = c(-1,1))+                                                      
  theme_classic()+
  theme(text = element_text(family = "Arial"))

graph2ppt(file="working_data/figures/neuron.dge.plots.pptx",width=7.6,height=7.1,append=TRUE)

outwb=createWorkbook()
sigdat=sigdat[order(sigdat$log2FoldChange,decreasing = T),]
goup = gprofiler(query=rownames(sigdat)[sigdat$log2FoldChange > 0], 
                 max_set_size = 1000,
                 min_isect_size = 2,
                 correction_method = "fdr",
                 exclude_iea = T,
                 hier_filtering = "moderate", 
                 custom_bg = rownames(dge.deseq2), 
                 src_filter = c("GO","KEGG"),
                 ordered_query = T)
goup=goup[order(goup$p.value),]
sheet=createSheet(outwb,sheetName = "neuron.up")
addDataFrame(goup,sheet,row.names = F)

sigdat=sigdat[order(sigdat$log2FoldChange,decreasing = F),]
godown = gprofiler(query=rownames(sigdat)[sigdat$log2FoldChange < 0], 
                   max_set_size = 1000,
                   min_isect_size = 2,
                   correction_method = "fdr",
                   exclude_iea = T,
                   hier_filtering = "moderate", 
                   custom_bg = rownames(dge.deseq2), 
                   src_filter = c("GO","KEGG"),
                   ordered_query = T)
godown=godown[order(godown$p.value),]
sheet=createSheet(outwb,sheetName = "neuron.down")
addDataFrame(godown,sheet,row.names = F)

saveWorkbook(outwb,file = "working_data/dge/dge.go.neuron.v2.moreCovariates.xlsx")

goup=read.xlsx(file = "working_data/dge/dge.go.neuron.v2.moreCovariates.xlsx",sheetName = "neuron.up")
godown=read.xlsx(file = "working_data/dge/dge.go.neuron.v2.moreCovariates.xlsx",sheetName = "neuron.down")
goup$FDR=-log10(goup$p.value);goup=goup[1:10,];goup$dirc="up"
godown$FDR=log10(godown$p.value);godown=godown[1:min(10,nrow(godown)),];godown$dirc="down"
go=rbind(goup,godown);go=go[order(go$FDR,decreasing = T),]

ggplot(go, aes(x=reorder(term.name, FDR), y=FDR,fill=dirc)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("up" = "gold","down"="blue"))+
  coord_flip() + 
  xlab("") + 
  ylab("Signed -log10(FDR)") +
  geom_hline(yintercept=c(log10(0.05),-log10(0.05)), lty=2, color="red")+
  theme(text = element_text(family = "Arial",size = 15,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20))+
  guides(fill=F)
graph2ppt(file="working_data/figures/neuron.dge.plots.pptx",width=10.5,height=7.1,append=TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)
dat1=as.data.frame(dge.deseq2)
geneList=dat1$log2FoldChange
names(geneList)=rownames(dat1)
geneList=sort(geneList,decreasing = T)

ego <- gseGO(geneList     = geneList,
             OrgDb        = org.Hs.eg.db,
             keyType = "ENSEMBL",
             ont          = "BP",
             nPerm        = 500,
             minGSSize    = 2,
             maxGSSize    = 1000,
             pvalueCutoff = 1,
             verbose      = TRUE)

ego.table=ego@result
ego.table=ego.table[complete.cases(ego.table),]
ego.table=ego.table[order(ego.table$NES,decreasing = T),]
ego.plot=rbind(head(ego.table,10),tail(ego.table,10))

ggplot(ego.plot, aes(x=reorder(Description, NES), y=NES,fill=p.adjust)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  xlab("") + 
  ylab("NES") +
  theme(text = element_text(family = "Arial",size = 15,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20))
graph2ppt(file="working_data/figures/neuron.dge.plots.pptx",width=10.5,height=7.1,append=TRUE)



# oligo
load("working_data/dge/dge.oligo.deseq2.v2.moreCovariates.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]
sigdat=dat[abs(dat$log2FoldChange) >= 0.58 & dat$padj < 0.05,]
nosigdat=dat[!rownames(dat) %in% rownames(sigdat),]
sigdat$genename=annot$genename[match(rownames(sigdat),annot$geneid)]
sigdat$genename[(sigdat$log2FoldChange >= -1.8 & sigdat$log2FoldChange <= 1.5) & sigdat$padj > 0.006]=""
ggplot()+geom_point(data=sigdat,aes(log2FoldChange,-log10(padj)),size=2,color="black")+
  geom_point(data = nosigdat,aes(log2FoldChange,-log10(padj)),size=2,alpha=1,col="gray")+
  geom_text_repel(data = sigdat,aes(log2FoldChange,-log10(padj),label=genename),min.segment.length = unit(0, "lines"),col="red")+
  geom_hline(yintercept = 1.3)+
  geom_vline(xintercept = c(-0.58,0.58))+
  #  scale_color_manual(values = colors)+
  scale_x_continuous(trans = "asinh")+
  #scale_y_continuous(trans = "log10")+
  #  geom_vline(xintercept = c(-1,1))+                                                      
  theme_classic()+
  theme(text = element_text(family = "Arial"))

graph2ppt(file="working_data/figures/oligo.dge.plots.pptx",width=7.6,height=7.1,append=TRUE)

outwb=createWorkbook()
sigdat=sigdat[order(sigdat$log2FoldChange,decreasing = T),]
goup = gprofiler(query=rownames(sigdat)[sigdat$log2FoldChange > 0], 
                 max_set_size = 1000,
                 min_isect_size = 2,
                 correction_method = "fdr",
                 exclude_iea = T,
                 hier_filtering = "moderate", 
                 custom_bg = rownames(dge.deseq2), 
                 src_filter = c("GO","KEGG"),
                 ordered_query = T)
goup=goup[order(goup$p.value),]
sheet=createSheet(outwb,sheetName = "oligo.up")
addDataFrame(goup,sheet,row.names = F)
sigdat=sigdat[order(sigdat$log2FoldChange,decreasing = F),]
godown = gprofiler(query=rownames(sigdat)[sigdat$log2FoldChange < 0], 
                   max_set_size = 1000,
                   min_isect_size = 2,
                   correction_method = "fdr",
                   exclude_iea = T,
                   hier_filtering = "moderate", 
                   custom_bg = rownames(dge.deseq2), 
                   src_filter = c("GO","KEGG"),
                   ordered_query = T)
godown=godown[order(godown$p.value),]
sheet=createSheet(outwb,sheetName = "oligo.down")
addDataFrame(godown,sheet,row.names = F)

saveWorkbook(outwb,file = "working_data/dge/dge.go.oligo.v2.moreCovariates.xlsx")

goup=read.xlsx(file = "working_data/dge/dge.go.oligo.v2.moreCovariates.xlsx",sheetName = "oligo.up")
godown=read.xlsx(file = "working_data/dge/dge.go.oligo.v2.moreCovariates.xlsx",sheetName = "oligo.down")
goup$FDR=-log10(goup$p.value);goup=goup[1:min(10,nrow(goup)),];goup$dirc="up"
godown$FDR=log10(godown$p.value);godown=godown[1:10,];godown$dirc="down"
go=rbind(goup,godown);go=go[order(go$FDR,decreasing = T),]

ggplot(go[complete.cases(go),], aes(x=reorder(term.name, FDR), y=FDR,fill=dirc)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c("up" = "gold","down"="blue"))+
  coord_flip() + 
  xlab("") + 
  ylab("Signed -log10(FDR)") +
  geom_hline(yintercept=c(log10(0.05),-log10(0.05)), lty=2, color="red")+
  theme(text = element_text(family = "Arial",size = 15,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20))+
  guides(fill=F)
graph2ppt(file="working_data/figures/oligo.dge.plots.pptx",width=10.5,height=7.1,append=TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)
dat1=as.data.frame(dge.deseq2)
geneList=dat1$log2FoldChange
names(geneList)=rownames(dat1)
geneList=sort(geneList,decreasing = T)

ego <- gseGO(geneList     = geneList,
             OrgDb        = org.Hs.eg.db,
             keyType = "ENSEMBL",
             ont          = "BP",
             nPerm        = 500,
             minGSSize    = 2,
             maxGSSize    = 1000,
             pvalueCutoff = 1,
             verbose      = TRUE)

ego.table=ego@result
ego.table=ego.table[complete.cases(ego.table),]
ego.table=ego.table[order(ego.table$NES,decreasing = T),]
ego.plot=rbind(head(ego.table,10),tail(ego.table,10))

ggplot(ego.plot, aes(x=reorder(Description, NES), y=NES,fill=p.adjust)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  xlab("") + 
  ylab("NES") +
  theme(text = element_text(family = "Arial",size = 15,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20))
graph2ppt(file="working_data/figures/oligo.dge.plots.pptx",width=10.5,height=7.1,append=TRUE)
