rm(list=ls());gc()
options(stringsAsFactors = F)

library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)

# mds on log2CPM
load("working_data/summarizedExperiment/fc_2pass/se_neuron_CPM_outlierRemoved.RData")
se_neuron=se_neuron[,!se_neuron$Sample %in% c("AO_B02","6221_N")]
mds=cmdscale(dist(t(assays(se_neuron)$log2CPM)))
colnames(mds)=c("PC1","PC2")
dat=cbind(mds,as.data.frame(colData(se_neuron)))
dat$Diagnosis=factor(dat$Diagnosis,levels = c("Control","Autism"))
ggplot(dat)+
  geom_point(aes(PC1,PC2,col=Diagnosis,shape=Sex),size=4)+
  #  geom_text_repel(aes(PC1,PC2,label=rownames(dat)),size=4)+
  theme_classic()+
  theme(text = element_text(family = "Arial"))

graph2ppt(file="working_data/figures/neuron.pptx",width=8.5,height=6.5,append=TRUE)

# process wgcna
rm(list=ls());gc()
options(stringsAsFactors = F)
load("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/neuron.fc2pass.recut.RData")
datExpr=networks$datExpr
datMeta=networks$datMeta

mds=cmdscale(dist(t(datExpr)))
colnames(mds)=c("PC1","PC2")
dat=cbind(mds,datMeta)

ggplot(dat)+
  geom_point(aes(PC1,PC2,col=Diagnosis,shape=Sex),size=4)+
  #  geom_text_repel(aes(PC1,PC2,label=rownames(dat)),size=4)+
  theme_classic()+
  theme(text = element_text(family = "Arial"))

graph2ppt(file="working_data/figures/neuron.pptx",width=8.5,height=6.5,append=TRUE)

# modTrait
library(WGCNA)

geneTree=networks$tree
merged=networks$merged
modules=merged$colors
print(length(unique(modules))-1)
MEs=networks$MEs
kMEtable=networks$kMEtable
table(datMeta$Sample == colnames(datExpr))
datMeta$Diagnosis=factor(datMeta$Diagnosis,levels = c("Control","Autism"))

modTrait=data.frame()
for(i in 2:length(unique(modules))) {
  me = MEs$eigengenes[,i]
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  moduleColor = labels2colors(moduleNumber)
  #s = tryCatch(summary(lme(me ~ Genotype, data=datMeta, random = ~1 | individual))$tTable,error=function(e){NA})
  s = summary(lm(me ~ Diagnosis + Age + Sex,data=datMeta))$coefficients
  if (!is.na(s)){
    for(grp in c("Autism")) {
      rowID = paste0("Diagnosis", grp)
      # modTrait = rbind(modTrait,
      #                  data.frame(Module=moduleColor, moduleNumber= moduleNumber, Group=grp,
      #                             beta = s[rowID, "Value"], SE = s[rowID, "Std.Error"], t=s[rowID, "t-value"], p=s[rowID, "p-value"]))
      modTrait = rbind(modTrait,
                       data.frame(Module=moduleColor, moduleNumber= moduleNumber, Group=grp,
                                  beta = s[rowID, "Estimate"], SE = s[rowID, "Std. Error"], t=s[rowID, "t value"], p=s[rowID, "Pr(>|t|)"]))
    }
  }
}
modTrait$fdr=p.adjust(modTrait$p)
modTrait$signedLog10fdr = -log10(modTrait$fdr) * sign(modTrait$beta)
modTrait$signedLog10fdr[modTrait$fdr > .05] = 0
modTrait$text = signif(modTrait$beta, 1)
modTrait$text[modTrait$fdr > 0.05] = ""
modTrait$Module=factor(modTrait$Module,levels = unique(modTrait$Module))
p=(ggplot(modTrait, aes(x=Module,y=Group, label=text)) +
     geom_tile(aes(fill=signedLog10fdr),color="grey60") + scale_fill_gradient2(low = "blue", mid = "white", high = "red","[beta]\nsigned\n-log10FDR\n") + 
     geom_text(size=3, color="black") + ylab("") + xlab("") + theme(axis.text.x = element_text(angle=30, hjust=1), axis.text.y = element_text(size=14))+
     geom_tile(aes(x=Module,y=0.1),fill=modTrait$Module,height=0.1))
print(p)
#graph2ppt(file="working_data/figures/neuron.pptx",width=11,height=2.7,append=TRUE)

write.table(modTrait,file = "working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds1mod40/modTrait.neuron.fc2pass.txt",
            quote = F,
            sep="\t",
            row.names = F)

# go
library(gProfileR)

genes=rownames(datExpr)
merged=networks$merged
modules=merged$colors
MEs = networks$MEs
kMEtable=networks$kMEtable

genes=sapply(genes,function(x) str_split_fixed(x,"\\.",Inf)[1,1])
names(genes)=NULL
rownames(kMEtable)=sapply(rownames(kMEtable),function(x) str_split_fixed(x,"\\.",Inf)[1,1])

annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")

for(i in 2:length(unique(modules))){
  
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  moduleColor = labels2colors(moduleNumber)
  #if (!moduleColor %in% modTrait$Module[modTrait$fdr < 0.05]){next}
  moduleGenes = genes[modules==moduleNumber]
  go = gprofiler(query=moduleGenes[order(kMEtable[moduleGenes,i], decreasing = T)], 
                 max_set_size = 500, 
                 #min_set_size = 2,
                 min_isect_size = 2,
                 correction_method = "fdr",
                 exclude_iea = T,
                 hier_filtering = "moderate", 
                 custom_bg = genes, 
                 src_filter = c("GO","KEGG"),
                 ordered_query = T)
  #	go = go[go$overlap.size > 2,]
  #	go = go[order(go$p.value)[1:min(10,nrow(go))],]
  go = go[order(go$p.value),]
  go = go[go$p.value < 0.05,]
  go$intersection_genename=apply(go,1,function(x){
    a=str_split_fixed(x[14],",",Inf)[1,]
    b=annot$genename[match(a,annot$geneid)]
    paste(b,collapse = ",")
  })
  write.table(go,file=paste("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds1mod40/go/go.M",moduleNumber,"_",moduleColor,".txt",sep=""),quote=F,sep="\t",row.names=F)
}

# gsea
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(caret)
library(WGCNA)


genes=rownames(datExpr)
merged=networks$merged
modules=merged$colors
MEs = networks$MEs
kMEtable=networks$kMEtable

genes=sapply(genes,function(x) str_split_fixed(x,"\\.",Inf)[1,1])
names(genes)=NULL
rownames(kMEtable)=sapply(rownames(kMEtable),function(x) str_split_fixed(x,"\\.",Inf)[1,1])

modTrait = read.delim("working_data/wgcna/neuron.fc/modTrait.neuron.fc.manualRemove.txt")

hsGO.bp <- godata('org.Hs.eg.db', ont="BP", computeIC = FALSE)
hsGO.mf <- godata('org.Hs.eg.db', ont="MF", computeIC = FALSE)
hsGO.cc <- godata('org.Hs.eg.db', ont="CC", computeIC = FALSE)

annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")

for(i in 2:length(unique(modules))){
  
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  moduleColor = labels2colors(moduleNumber)
 # if (!moduleColor %in% modTrait$Module[modTrait$fdr < 0.05]){next}
 # if (!moduleColor %in% c("blue","royalblue","steelblue")){next}
  geneList=kMEtable[,i]
  names(geneList)=rownames(kMEtable)
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
                 pvalueCutoff = 0.1,
                 verbose      = TRUE)
    
    ego.table=ego@result
    ego.table=ego.table[complete.cases(ego.table),]
    ego.table=ego.table[ego.table$qvalues < 0.1,]
    if (nrow(ego.table) == 0){next}
    #hsGO <- godata('org.Hs.eg.db', ont=onto, computeIC = FALSE)
    if (onto == "BP"){
      go_cor=mgoSim(ego.table$ID,ego.table$ID,semData = hsGO.bp,measure = "Wang",combine = NULL)
    }
    if (onto == "MF"){
      go_cor=mgoSim(ego.table$ID,ego.table$ID,semData = hsGO.mf,measure = "Wang",combine = NULL)
    }
    if (onto == "CC"){
      go_cor=mgoSim(ego.table$ID,ego.table$ID,semData = hsGO.cc,measure = "Wang",combine = NULL)
    }
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
  write.table(ego.table.all,file=paste("working_data/wgcna/neuron.fc/gsea.manualRemove/go.M",moduleNumber,"_",moduleColor,".txt",sep=""),quote=F,sep="\t",row.names=F)
}

# membership

membership=data.frame(geneid=genes,module=labels2colors(modules))
annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")

membership=merge(membership,annot,by="geneid",all.x=T)
write.table(membership,file = "working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/membership.fc2pass.txt",quote = F,sep = "\t",row.names = F)

