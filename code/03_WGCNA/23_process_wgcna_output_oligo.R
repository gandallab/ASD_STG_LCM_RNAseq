rm(list=ls());gc()
options(stringsAsFactors = F)

library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)

# mds on log2CPM
load("working_data/summarizedExperiment/se_oligo_CPM_outlierRemoved.RData")

mds=cmdscale(dist(t(assays(se_oligo)$log2CPM)))
colnames(mds)=c("PC1","PC2")
dat=cbind(mds,as.data.frame(colData(se_oligo)))

ggplot(dat)+
  geom_point(aes(PC1,PC2,col=Diagnosis,shape=Sex),size=4)+
  #  geom_text_repel(aes(PC1,PC2,label=rownames(dat)),size=4)+
  theme_classic()+
  theme(text = element_text(family = "Arial"))

graph2ppt(file="working_data/figures/oligo.wgcna.v2.moreCovar.pptx",width=8.5,height=6.5)

# process wgcna
rm(list=ls());gc()
options(stringsAsFactors = F)
load("working_data/wgcna/oligo.v2.moreCovar//oligo.bicor.recut.RData")
load("working_data/wgcna/voom.forWGCNA.input.oligo.v2.RData")

mds=cmdscale(dist(t(datExpr)))
colnames(mds)=c("PC1","PC2")
dat=cbind(mds,datMeta)

ggplot(dat)+
  geom_point(aes(PC1,PC2,col=Diagnosis,shape=Sex),size=4)+
  #  geom_text_repel(aes(PC1,PC2,label=rownames(dat)),size=4)+
  theme_classic()+
  theme(text = element_text(family = "Arial"))

graph2ppt(file="working_data/figures/oligo.wgcna.v2.moreCovar.pptx",width=8.5,height=6.5,append=TRUE)
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
graph2ppt(file="working_data/figures/oligo.wgcna.v2.moreCovar.pptx",width=15,height=2.7,append=TRUE)
write.table(modTrait,file = "working_data/wgcna/oligo.v2.moreCovar//modTrait.oligo.bicor.txt",
            quote = F,
            sep="\t",
            row.names = F)

# go
library(gProfileR)
datExpr=networks$datExpr
genes=rownames(datExpr)
merged=networks$merged
modules=merged$colors
MEs = networks$MEs
kMEtable=networks$kMEtable

for(i in 2:length(unique(modules))){
  
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  moduleColor = labels2colors(moduleNumber)
  moduleGenes = genes[modules==moduleNumber]
  go = gprofiler(query=moduleGenes[order(kMEtable[moduleGenes,i], decreasing = T)], 
                 max_set_size = 1000,
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
  write.table(go,file=paste("working_data/wgcna/oligo.v2.moreCovar//go/go.M",moduleNumber,"_",moduleColor,".txt",sep=""),quote=F,sep="\t",row.names=F)
}
