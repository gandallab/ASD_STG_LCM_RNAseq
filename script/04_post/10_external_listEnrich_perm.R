options(stringsAsFactors = F);rm(list=ls());gc()
library(tidyverse)
library(clusterProfiler)
library(WGCNA)

load("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/neuron.fc2pass.recut.RData")
datExpr=networks$datExpr
allGenes=rownames(datExpr)
allGenes=sapply(allGenes,function(x) str_split_fixed(x,"\\.",Inf)[1,1])
names(allGenes)=NULL
merged=networks$merged
modules=merged$colors
length(unique(modules))
MEs = networks$MEs

geneFeature=read.delim("D:/references/gencodev19/gene.length.gc.txt",header = T,sep = "\t")

geneFeature1=geneFeature[geneFeature$geneid %in% allGenes,]
colnames(geneFeature1)=c("V1","V2","V3")
rm(geneFeature);gc()

#gmt=read.gmt("F:/datasets/genesets/Mariani_Vaccarino_Cell_FOXG1_TableS2_WGCNAmodules.gmt")
#gmt=read.gmt("working_data/parikshak.nature.cortex.module.gmt")
#gmt=read.gmt("D:/datasets/genesets/voineagu.modules.gmt")
gmt=read.gmt("D:/datasets/genesets/houseLists03242020.gmt")
genename=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt",header = T,sep = "\t")
gmt$geneid=genename$geneid[match(gmt$gene,genename$genename)]
permTable=data.frame()


for(i in 2:length(unique(modules))) {
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  moduleColor = labels2colors(moduleNumber)
  moduleGenes = allGenes[modules==moduleNumber]
  
  permut=list()
  for (j in moduleGenes){
    if (!j %in% geneFeature1$V1){next}
    len=geneFeature1$V2[geneFeature1$V1 == j]
    gc=geneFeature1$V3[geneFeature1$V1 == j]
    current=geneFeature1 %>% filter(abs(V2-len)/len < 0.2 & abs(V3-gc)/gc < 0.2)
    permut[[j]]=current$V1
  }
  
  permut1=list();ii=0
  for (j in 1:1000){
    current=c()
    for (k in 1:length(permut)){
      current=c(current,permut[[k]][sample(1:length(permut[[k]]),1)])
    }
    #    if (length(unique(current)) != length(permut)){next}
    ii=ii+1
    permut1[[ii]]=unique(current)
  }
  
  
  
  permutRecorder=data.frame()
  for (j in unique(gmt$ont)){
    genesets=unique(gmt$geneid[gmt$ont==j])
    genesets=genesets[!is.na(genesets)]
    recorder=c()
    for (k in 1:length(permut1)){
      recorder=c(recorder,length(intersect(genesets,permut1[[k]])))
    }
    downsample=sample(moduleGenes,median(lengths(permut1)),replace = F)
    permutRecorder=rbind(permutRecorder,data.frame(category=j,mean=mean(recorder),sd=sd(recorder),actual=length(intersect(genesets,downsample))))
  }
  
  
  pvalue=apply(permutRecorder,1,function(x){
    1-pnorm(as.numeric(x[4]),mean = as.numeric(x[2]),sd=as.numeric(x[3]))
  })
  
  fdr=p.adjust(pvalue,method = "fdr")
  
  zscore=apply(permutRecorder,1,function(x){
    (as.numeric(x[4]) - as.numeric(x[2]))/ as.numeric(x[3])
  })
  
  permutRecorder$zscore=zscore
  permutRecorder$fdr=fdr
  permutRecorder$moduleNumber=moduleNumber
  permutRecorder$moduleColor=moduleColor
  
  permTable=rbind(permTable,permutRecorder)
}

permTable$Module=paste(permTable$moduleNumber,permTable$moduleColor,sep = "")
permTable$signedLog10fdr = -log10(permTable$fdr) * sign(permTable$zscore)
permTable$signedLog10fdr[permTable$fdr>0.05]=0
permTable$signedLog10fdr[permTable$signedLog10fdr == Inf]=max(permTable$signedLog10fdr[permTable$signedLog10fdr != Inf])

permTable$text = round(permTable$zscore, 1)
permTable$text[permTable$fdr>0.05] = ""
permTable$Module1=reorder(permTable$Module,permTable$moduleNumber)
ggplot(permTable, aes(x=Module1,y=category, label=text)) +
  geom_tile(aes(fill=signedLog10fdr),color="grey60") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red","signed\n-log10FDR\n") + 
  geom_text(size=3, color="black") + ylab("") + xlab("") + 
  theme(axis.text.x = element_text(angle=30, hjust=1), axis.text.y = element_text(size=7))+
  geom_tile(aes(x=Module,y=0.1),fill=permTable$moduleColor,height=0.1)

save(permTable,file = "working_data/listEnrich/permutTable.1000times.neuron_fc2passModule.newhouseList.RData")
permTable$fdr[permTable$fdr == 0]=min(permTable$fdr[permTable$fdr > 0])
