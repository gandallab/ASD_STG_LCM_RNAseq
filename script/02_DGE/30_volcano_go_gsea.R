rm(list=ls());gc()
options(stringsAsFactors = F)

library(ggplot2)
library(ggrepel)
library(scales)
#library(export)
#library(xlsx)
library(gProfileR)
library(tidyverse)

asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")


# blocks
load("working_data/dge/dge.blocks.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
dat$genename=annot$genename[match(rownames(dat),annot$geneid)]
dat$genetype=annot$genetype[match(rownames(dat),annot$geneid)]
dat.sno=dat[dat$genetype == "snRNA",]
dat.sno=dat.sno[!is.na(dat.sno$genetype),]
dat.sno$genename[dat.sno$pvalue > 0.01]=""
dat.nonsno=dat[dat$genetype != "snRNA",]

sigdat=dat[abs(dat$log2FoldChange) >= 0.58 & dat$padj < 0.05,]
nosigdat=dat[!rownames(dat) %in% rownames(sigdat),]
sigdat$genename=annot$genename[match(rownames(sigdat),annot$geneid)]
sigdat$genename[(sigdat$log2FoldChange > 0 & sigdat$log2FoldChange <= 1.2) & sigdat$padj > 0.006]=""

fdrdat=dat[dat$padj < 0.05,]
pdat=dat[dat$padj > 0.05 & dat$pvalue < 0.05,]
restdat=dat[dat$pvalue > 0.05,]
labeldat=dat[dat$genename %in% c("HSPA1A","HSPA1B","DNAJB1","DNAJB4","HSPB1","HSPB8","BAG3","PTGES3","KCNH3","KCNIP1"),]


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
  theme(text = element_text(family = "Arial"),
        axis.title = element_text(size = 17,face = "bold"),
        axis.text = element_text(size=15))

ggplot()+
  geom_point(data=dat.nonsno,aes(log2FoldChange,-log10(pvalue)),color="gray",size=1.5)+
  geom_point(data=dat.sno,aes(log2FoldChange,-log10(pvalue)),color="red",size=2)+
  geom_text_repel(data = dat.sno,aes(log2FoldChange,-log10(pvalue),label=genename),min.segment.length = unit(0, "lines"),col="red")+
  scale_x_continuous(trans = "asinh")+
  theme_classic()+
  theme(text = element_text(family = "Arial"))

ggplot()+
  geom_point(data=fdrdat[fdrdat$log2FoldChange > 0,],aes(log2FoldChange,-log10(pvalue)),color="gold",size=2.5)+
  geom_point(data=fdrdat[fdrdat$log2FoldChange < 0,],aes(log2FoldChange,-log10(pvalue)),color="royalblue",size=2.5)+
  geom_point(data=pdat,aes(log2FoldChange,-log10(pvalue)),color="gray",size=2)+
  geom_point(data = restdat,aes(log2FoldChange,-log10(pvalue)),color="gray",size=1.5)+
  geom_text_repel(data = fdrdat[(!fdrdat$genename %in% labeldat$genename) & (fdrdat$padj < 0.025),],aes(log2FoldChange,-log10(pvalue),label=genename),segment.size = 0.1,col="black",size=2.5)+
  geom_text_repel(data = labeldat,aes(log2FoldChange,-log10(pvalue),label=genename),segment.size = 0.1,col="red",fontface="bold",size=3)+
  scale_x_continuous(trans = "asinh")+
  theme_classic(base_size = 16)+
  theme(text = element_text(family = "Arial"),
        axis.title = element_text(size = 17,face = "bold"),
        axis.text = element_text(size=15))


graph2ppt(file="working_data/figures/blocks.dge.plots.pptx",width=7.6,height=7.1,append=TRUE)

outwb=createWorkbook()
sigdat=dat[dat$padj < 0.05,]
sigdat=sigdat[order(sigdat$log2FoldChange,decreasing = T),]
goup = gprofiler(query=rownames(sigdat)[sigdat$log2FoldChange > 0], 
               max_set_size = 1000,
               min_isect_size = 2,
               correction_method = "fdr",
               exclude_iea = T,
               hier_filtering = "moderate", 
               custom_bg = rownames(dat), 
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
                 custom_bg = rownames(dat), 
                 src_filter = c("GO","KEGG"),
                 ordered_query = T)
godown=godown[order(godown$p.value),]
sheet=createSheet(outwb,sheetName = "blocks.down")
addDataFrame(godown,sheet,row.names = F)

saveWorkbook(outwb,file = "working_data/dge/dge.go.blocks.v3.xlsx")


goup=read.xlsx(file = "working_data/dge/dge.go.blocks.v3.xlsx",sheetName = "blocks.up")
godown=read.xlsx(file = "working_data/dge/dge.go.blocks.v3.xlsx",sheetName = "blocks.down")
goup$FDR=-log10(goup$p.value);goup=goup[1:15,];goup$dirc="up"
godown$FDR=log10(godown$p.value);godown=godown[1:15,];godown$dirc="down"
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
graph2ppt(file="working_data/figures/blocks.pptx",width=11.2,height=7.4,append=TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(caret)
annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")
dat1=as.data.frame(dge.deseq2)
rownames(dat1)=sapply(rownames(dat1),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
geneList=dat1$log2FoldChange
names(geneList)=rownames(dat1)
geneList=sort(geneList,decreasing = T)
ego.table.all=data.frame()
for (onto in c("BP","MF","CC")){
  ego <- gseGO(geneList     = geneList,
                OrgDb        = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont          = onto,
                nPerm        = 5000,
                minGSSize    = 50,
                maxGSSize    = 300,
                pvalueCutoff = 1,
                verbose      = TRUE)
  
  ego.table=ego@result
  ego.table=ego.table[complete.cases(ego.table),]
  ego.table=ego.table[ego.table$qvalues < 0.1,]
  if (nrow(ego.table) == 0){next}
  hsGO <- godata('org.Hs.eg.db', ont=onto,computeIC = F)
  go_cor=mgoSim(ego.table$ID,ego.table$ID,semData = hsGO,measure = "Wang",combine = NULL)
  to_remove=findCorrelation(go_cor,cutoff = 0.7,names = T)
  ego.table=ego.table[!ego.table$ID %in% to_remove,]
  ego.table$onto=onto
  ego.table.all=rbind(ego.table.all,ego.table)
}
  

ego.table.all=ego.table.all[order(ego.table.all$NES,decreasing = T),]
ego.table.all$core_enrichment_genename=apply(ego.table.all,1,function(x){
  a=str_split_fixed(x[11],"/",Inf)[1,]
  b=annot$genename[match(a,annot$geneid)]
  paste(b,collapse = ",")
})
write.table(ego.table.all,file = "working_data/dge/gsea/blocks.dge.gsea.txt",quote = F,sep = "\t",row.names = F)

ego.plot=rbind(head(ego.table.all,10),tail(ego.table.all,10))


ggplot(ego.plot, aes(x=reorder(Description, NES), y=NES,fill=(NES < 0))) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  xlab("") + 
  ylab("NES") +
  scale_fill_manual(values = c("gold","royalblue"))+
  theme(text = element_text(family = "Arial",size = 15,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20))
graph2ppt(file="working_data/figures/blocks.dge.plots.pptx",width=10.5,height=7.1,append=TRUE)




# neuron
rm(list=ls());gc()
options(stringsAsFactors = F)
load("working_data/dge/dge.neuron.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
dat$genename=annot$genename[match(rownames(dat),annot$geneid)]
dat$genetype=annot$genetype[match(rownames(dat),annot$geneid)]
dat.sno=dat[dat$genetype == "snRNA",]
dat.sno=dat.sno[!is.na(dat.sno$genetype),]
dat.sno$genename[dat.sno$pvalue > 0.05]=""
dat.nonsno=dat[dat$genetype != "snRNA",]

fdrdat=dat[dat$padj < 0.05,]
pdat=dat[dat$padj > 0.05 & dat$pvalue < 0.05,]
restdat=dat[dat$pvalue > 0.05,]
labeldat=dat[dat$genename %in% c("GAD1","GAD2","JUN","JUNB","FOS","SOX9","S1PR1","PPP1R16B","KCNJ2","NFKBID"),]

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
  theme(text = element_text(family = "Arial"),
        axis.title = element_text(size = 17,face = "bold"),
        axis.text = element_text(size=15))

ggplot()+
  geom_point(data=dat.nonsno,aes(log2FoldChange,-log10(pvalue)),color="gray",size=1.5)+
  geom_point(data=dat.sno,aes(log2FoldChange,-log10(pvalue)),color="red",size=2)+
  geom_hline(linetype=2,size = 1,col="black",yintercept = c(1.3,-log10(max(dat$pvalue[dat$padj < 0.05]))))+
  geom_text_repel(data = dat.sno,aes(log2FoldChange,-log10(pvalue),label=genename),min.segment.length = unit(0, "lines"),col="red")+
  scale_x_continuous(trans = "asinh")+
  theme_classic()+
  theme(text = element_text(family = "Arial"),
        axis.title = element_text(size = 17,face = "bold"),
        axis.text = element_text(size=15))

ggplot()+
  geom_point(data=fdrdat[fdrdat$log2FoldChange > 0,],aes(log2FoldChange,-log10(pvalue)),color="gold",size=2.5)+
  geom_point(data=fdrdat[fdrdat$log2FoldChange < 0,],aes(log2FoldChange,-log10(pvalue)),color="royalblue",size=2.5)+
  geom_point(data=pdat,aes(log2FoldChange,-log10(pvalue)),color="gray",size=2)+
  geom_point(data = restdat,aes(log2FoldChange,-log10(pvalue)),color="gray",size=1.5)+
  geom_text_repel(data = fdrdat[!fdrdat$genename %in% labeldat$genename,],aes(log2FoldChange,-log10(pvalue),label=genename),segment.size = 0.1,col="black",size=3)+
  geom_text_repel(data = labeldat,aes(log2FoldChange,-log10(pvalue),label=genename),segment.size = 0.1,col="red",fontface="bold",size=3.5)+
  scale_x_continuous(trans = "asinh")+
  theme_classic(base_size = 16)+
  theme(text = element_text(family = "Arial"),
        axis.title = element_text(size = 17,face = "bold"),
        axis.text = element_text(size=15))

graph2ppt(file="working_data/figures/neuron.dge.plots.pptx",width=7.6,height=7.1,append=TRUE)

outwb=createWorkbook()
sigdat=dat[dat$padj < 0.05,]
sigdat=sigdat[order(sigdat$log2FoldChange,decreasing = T),]
goup = gprofiler(query=rownames(sigdat)[sigdat$log2FoldChange > 0], 
                 max_set_size = 1000,
                 min_isect_size = 2,
                 correction_method = "fdr",
                 exclude_iea = T,
                 hier_filtering = "moderate", 
                 custom_bg = rownames(dat), 
                 src_filter = c("GO","KEGG"),
                 ordered_query = T)
goup=goup[order(goup$p.value),]
sheet=createSheet(outwb,sheetName = "neuron.up")
addDataFrame(goup,sheet,row.names = F)

sigdat=sigdat[order(sigdat$log2FoldChange,decreasing = F),]
sigdat$genename=annot$genename[match(rownames(sigdat),annot$geneid)]
godown = gprofiler(query=rownames(sigdat)[sigdat$log2FoldChange < 0], 
                   max_set_size = 1000,
                   min_isect_size = 2,
                   correction_method = "fdr",
                   exclude_iea = T,
                   hier_filtering = "moderate", 
                   custom_bg = rownames(dat), 
                   src_filter = c("GO","KEGG"),
                   ordered_query = T)
godown = gprofiler(query=sigdat$genename[sigdat$log2FoldChange < 0], 
                   max_set_size = 1000,
                   min_isect_size = 2,
                   correction_method = "fdr",
                   exclude_iea = T,
                   hier_filtering = "moderate", 
                   custom_bg = annot$genename[match(rownames(dge.deseq2),annot$geneid)], 
                   src_filter = c("GO","KEGG"),
                   ordered_query = T)
godown=godown[order(godown$p.value),]
sheet=createSheet(outwb,sheetName = "neuron.down")
addDataFrame(godown,sheet,row.names = F)

saveWorkbook(outwb,file = "working_data/dge/dge.go.neuron.v3.xlsx")

goup=read.xlsx(file = "working_data/dge/dge.go.neuron.v3.xlsx",sheetName = "neuron.up")
godown=read.xlsx(file = "working_data/dge/dge.go.neuron.v3.xlsx",sheetName = "neuron.down")
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
graph2ppt(file="working_data/figures/neuron.pptx",width=10.5,height=7.1,append=TRUE)

library(clusterProfiler)
library(org.Hs.eg.db)
rm(list=ls());gc()
options(stringsAsFactors = F)
load("working_data/dge/dge.neuron.fc2pass.RData")
dat1=as.data.frame(dge.deseq2)
rownames(dat1)=sapply(rownames(dat1),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
geneList=dat1$log2FoldChange
#geneList= 1 - dat1$pvalue
names(geneList)=rownames(dat1)
geneList=sort(geneList,decreasing = T)
ego.table.all=data.frame()
for (onto in c("BP","MF","CC")){
  ego <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               keyType = "ENSEMBL",
               ont          = onto,
               nPerm        = 5000,
               minGSSize    = 50,
               maxGSSize    = 300,
               pvalueCutoff = 1,
               verbose      = TRUE)
  
  ego.table=ego@result
  ego.table=ego.table[complete.cases(ego.table),]
  ego.table=ego.table[ego.table$qvalues < 0.1,]
  if (nrow(ego.table) == 0){next}
  # hsGO <- godata('org.Hs.eg.db', ont=onto, computeIC = F)
  # go_cor=mgoSim(ego.table$ID,ego.table$ID,semData = hsGO,measure = "Wang",combine = NULL)
  # to_remove=findCorrelation(go_cor,cutoff = 0.7,names = T)
  # ego.table=ego.table[!ego.table$ID %in% to_remove,]
  ego.table$onto=onto
  ego.table.all=rbind(ego.table.all,ego.table)
}


ego.table.all=ego.table.all[order(ego.table.all$NES,decreasing = T),]
annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")
ego.table.all$core_enrichment_genename=apply(ego.table.all,1,function(x){
  a=str_split_fixed(x[11],"/",Inf)[1,]
  b=annot$genename[match(a,annot$geneid)]
  paste(b,collapse = ",")
})
write.table(ego.table.all,file = "working_data/dge/gsea/neuron.dge.gsea.txt",quote = F,sep = "\t",row.names = F)

ego.plot=rbind(head(ego.table.all,15),tail(ego.table.all,15))


ggplot(ego.plot, aes(x=reorder(Description, NES), y=NES,fill=(NES < 0))) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  xlab("") + 
  ylab("NES") +
  scale_fill_manual(values = c("gold","royalblue"))+
  theme(text = element_text(family = "Arial",size = 15,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20))
graph2ppt(file="working_data/figures/neuron.dge.plots.pptx",width=10.5,height=7.1,append=TRUE)



# oligo
load("working_data/dge/dge.oligo.fc_2pass.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
dat$genename=annot$genename[match(rownames(dat),annot$geneid)]
dat$genetype=annot$genetype[match(rownames(dat),annot$geneid)]
dat.sno=dat[dat$genetype == "snoRNA",]
dat.sno=dat.sno[!is.na(dat.sno$genetype),]
dat.sno$genename[dat.sno$pvalue > 0.01]=""
dat.nonsno=dat[dat$genetype != "snoRNA",]


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


ggplot()+
  geom_point(data=dat.nonsno,aes(log2FoldChange,-log10(pvalue)),color="gray",size=1.5)+
  geom_point(data=dat.sno,aes(log2FoldChange,-log10(pvalue)),color="red",size=2)+
  geom_text_repel(data = dat.sno,aes(log2FoldChange,-log10(pvalue),label=genename),min.segment.length = unit(0, "lines"),col="red")+
  scale_x_continuous(trans = "asinh")+
  theme_classic()+
  theme(text = element_text(family = "Arial"))

graph2ppt(file="working_data/figures/oligo.dge.plots.pptx",width=7.6,height=7.1,append=TRUE)

outwb=createWorkbook()
sigdat=dat[dat$padj < 0.05,]
sigdat=sigdat[order(sigdat$log2FoldChange,decreasing = T),]
goup = gprofiler(query=rownames(sigdat)[sigdat$log2FoldChange > 0], 
                 max_set_size = 1000,
                 min_isect_size = 2,
                 correction_method = "fdr",
                 exclude_iea = T,
                 hier_filtering = "moderate", 
                 custom_bg = rownames(dat), 
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

saveWorkbook(outwb,file = "working_data/dge/dge.go.oligo.v3.xlsx")

goup=read.xlsx(file = "working_data/dge/dge.go.oligo.v3.xlsx",sheetName = "oligo.up")
godown=read.xlsx(file = "working_data/dge/dge.go.oligo.v3.xlsx",sheetName = "oligo.down")
goup$FDR=-log10(goup$p.value);goup=goup[1:min(10,nrow(goup)),];goup$dirc="up"
godown$FDR=log10(godown$p.value);godown=godown[1:3,];godown$dirc="down"
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
graph2ppt(file="working_data/figures/oligo.pptx",width=10.5,height=7.1,append=TRUE)

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









# blocks interaction
rm(list=ls());gc()
options(stringsAsFactors = F)

library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(caret)
annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")
dat1=read.delim("working_data/dge/dge.neuron.fc2pass.interaction.txt")
rownames(dat1)=sapply(rownames(dat1),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
geneList=dat1$log2FoldChange
names(geneList)=rownames(dat1)
geneList=sort(geneList,decreasing = T)
ego.table.all=data.frame()
for (onto in c("BP","MF","CC")){
  ego <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               keyType = "ENSEMBL",
               ont          = onto,
               nPerm        = 5000,
               minGSSize    = 50,
               maxGSSize    = 300,
               pvalueCutoff = 1,
               verbose      = TRUE)
  
  ego.table=ego@result
  ego.table=ego.table[complete.cases(ego.table),]
  ego.table=ego.table[ego.table$qvalues < 0.1,]
  if (nrow(ego.table) == 0){next}
  hsGO <- godata('org.Hs.eg.db', ont=onto,computeIC = F)
  go_cor=mgoSim(ego.table$ID,ego.table$ID,semData = hsGO,measure = "Wang",combine = NULL)
  to_remove=findCorrelation(go_cor,cutoff = 0.7,names = T)
  ego.table=ego.table[!ego.table$ID %in% to_remove,]
  ego.table$onto=onto
  ego.table.all=rbind(ego.table.all,ego.table)
}


ego.table.all=ego.table.all[order(ego.table.all$NES,decreasing = T),]
ego.table.all$core_enrichment_genename=apply(ego.table.all,1,function(x){
  a=str_split_fixed(x[11],"/",Inf)[1,]
  b=annot$genename[match(a,annot$geneid)]
  paste(b,collapse = ",")
})
write.table(ego.table.all,file = "working_data/dge/gsea/blocks_interaction.dge.gsea.txt",quote = F,sep = "\t",row.names = F)

ego.table.all=read.delim("working_data/dge/gsea/blocks_interaction.dge.gsea.txt")
ego.table.all=ego.table.all[ego.table.all$onto != "CC",]
ego.table.all=ego.table.all[ego.table.all$Description != "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",]
ego.plot=rbind(head(ego.table.all,5),tail(ego.table.all,5))


ggplot(ego.plot, aes(x=reorder(Description, NES), y=NES,fill=(NES < 0))) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  xlab("") + 
  ylab("NES") +
  scale_fill_manual(values = c("gold","royalblue"))+
  theme(text = element_text(family = "Arial",size = 15,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20))

ggplot(ego.plot, aes(x=reorder(Description, NES), y=NES,fill=qvalues)) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  xlab("") + 
  ylab("NES") +
  #scale_fill_manual(values = c("gold","royalblue"))+
  theme(text = element_text(family = "Arial",size = 15,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20))

# blocks interaction Gprofiler

dat=read.delim("working_data/dge/dge.blocks.fc2pass.interaction.txt")
dat1=dat[complete.cases(dat),]
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
dat$genename=annot$genename[match(rownames(dat),annot$geneid)]

sigdat=dat1[dat1$padj < 0.05,]
sigdat=sigdat[order(sigdat$log2FoldChange,decreasing = T),]
goup = gprofiler(query=rownames(sigdat)[sigdat$log2FoldChange > 0], 
                 max_set_size = 1000,
                 min_isect_size = 2,
                 correction_method = "fdr",
                 exclude_iea = T,
                 hier_filtering = "moderate", 
                 custom_bg = rownames(dat), 
                 src_filter = c("GO","KEGG"),
                 ordered_query = T)
goup=goup[order(goup$p.value),]


sigdat=sigdat[order(sigdat$log2FoldChange,decreasing = F),]
godown = gprofiler(query=rownames(sigdat)[sigdat$log2FoldChange < 0], 
                   max_set_size = 1000,
                   min_isect_size = 2,
                   correction_method = "fdr",
                   exclude_iea = T,
                   hier_filtering = "moderate", 
                   custom_bg = rownames(dat), 
                   src_filter = c("GO","KEGG"),
                   ordered_query = T)
godown=godown[order(godown$p.value),]
write.table(godown,file = "working_data/dge/block_interaction_down_go.txt",quote = F,sep = "\t",row.names = F)
