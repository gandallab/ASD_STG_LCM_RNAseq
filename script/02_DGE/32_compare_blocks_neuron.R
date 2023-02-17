rm(list=ls());gc()
options(stringsAsFactors = F)

load("working_data/dge/dge.blocks.fc2pass.RData")
block=as.data.frame(dge.deseq2)
colnames(block)=paste("block",colnames(block),sep = "_")
load("working_data/dge/dge.neuron.fc2pass.RData")
neuron=as.data.frame(dge.deseq2)
colnames(neuron)=paste("neuron",colnames(neuron),sep = "_")
dat=merge(block,neuron,by="row.names")
rownames(dat)=sapply(dat$Row.names,function(x) str_split_fixed(x,"\\.",Inf)[1,1])
dat$Row.names=NULL
dat$diff.log2FC = dat$neuron_log2FoldChange - dat$block_log2FoldChange
#dat$diff.log2FC = log2((2 ^ dat$neuron_log2FoldChange) - (2 ^ dat$block_log2FoldChange))

annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")

dat1=merge(annot,dat,by.x="geneid",by.y="row.names")

library(clusterProfiler)
library(org.Hs.eg.db)

geneList=dat1$diff.log2FC
names(geneList)=dat1$geneid
geneList=sort(geneList,decreasing = T)
ego.table.all=data.frame()
for (onto in c("BP","MF","CC")){
  ego <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               keyType = "ENSEMBL",
               ont          = onto,
               nPerm        = 500,
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
write.table(ego.table.all,file = "working_data/dge/gsea/contrast.blocks.neuron.txt",quote = F,sep = "\t",row.names = F)
ego.plot=ego.plot=rbind(head(ego.table.all,15),tail(ego.table.all,15))


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
