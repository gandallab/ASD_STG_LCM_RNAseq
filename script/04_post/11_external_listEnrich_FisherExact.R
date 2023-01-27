rm(list = ls());gc();options(stringsAsFactors = F)

OR <- function(q,k,m,t) {
  
  ## 2 x 2 table:
  
  ##         inTest   !inTest
  
  ## inRef     q        k
  
  ## !inRef    m        t
  
  
  
  fisher.out <- fisher.test(matrix(c(q, k-q, m-q, t-m-k+q), 2, 2),conf.int=TRUE)
  
  OR <- fisher.out$estimate
  
  pval <- fisher.out$p.value
  
  upCI <- fisher.out$conf.int[1]
  
  downCI <- fisher.out$conf.int[2]
  
  
  
  output <- c(OR,pval,upCI,downCI)
  
  names(output) <- c("OR","Fisher p","-95%CI","+95%CI")
  
  return(output)
  
}



## count overlaps and run the analysis

ORA <- function(testpath,refpath,testbackground,refbackground) {
  
  testpath = testpath[testpath %in% testbackground]
  
  refpath = refpath[refpath %in% refbackground]
  
  q <- length(intersect(testpath,refpath)) ## overlapped pathway size
  
  k <- length(intersect(refpath,testbackground))  ## input gene set
  
  m <- length(intersect(testpath,refbackground)) ## input module
  
  t <- length(intersect(testbackground,refbackground)) ## Total assessed background (intersect reference and test backgrounds)
  
  
  
  empvals <- OR(q,k,m,t)
  
  
  
  tmpnames <- names(empvals)
  
  empvals <- as.character(c(empvals,q,k,m,t,100*signif(q/k,3)))
  
  names(empvals) <- c(tmpnames,"Overlap","Reference List","Input List","Background","% List Overlap")
  
  return(empvals)
  
}

library(clusterProfiler)
#gmt=read.gmt("F:/datasets/genesets/Mariani_Vaccarino_Cell_FOXG1_TableS2_WGCNAmodules.gmt")
#gmt=read.gmt("F:/datasets/genesets/voineagu.modules.gmt")
gmt=read.gmt("D:/datasets/genesets/houseLists02062020.gmt")
genename=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt",header = T,sep = "\t")
gmt$geneid=genename$geneid[match(gmt$gene,genename$genename)]
pLI=read.delim("D:/datasets/exac/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt",header = T,sep = "\t")
pLI_geneid=unique(genename$geneid[match(pLI$gene,genename$genename)])


library(WGCNA)
load("working_data/wgcna/neuron.fc/neuron.fc.manualRemove.recut.RData")
datExpr=networks$datExpr
allGenes=rownames(datExpr)
allGenes=sapply(allGenes,function(x) str_split_fixed(x,"\\.",Inf)[1,1])
names(allGenes)=NULL
merged=networks$merged
modules=merged$colors
length(unique(modules))
MEs = networks$MEs

refpath=list();refbackgroud=list()
for (i in unique(gmt$ont)){
  genesets=gmt$geneid[gmt$ont==i]
  genesets=genesets[!is.na(genesets)]
  refpath[[i]]=genesets
  if (grepl("pLI",i)) {refbackgroud[[i]]=pLI_geneid}
  else{refbackgroud[[i]]=allGenes}
}



ORtable=data.frame(matrix(NA,nrow = length(unique(modules))-1,ncol=length(refpath)+1))
colnames(ORtable)=c("Module",names(refpath))
ptable=data.frame(matrix(NA,nrow = length(unique(modules))-1,ncol=length(refpath)+1))
colnames(ptable)=c("Module",names(refpath))

for(i in 2:length(unique(modules))) {
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  moduleColor = labels2colors(moduleNumber)
  moduleGenes = allGenes[modules==moduleNumber]
  ORtable[i-1,1]=paste(moduleNumber,moduleColor,sep = "")
  ptable[i-1,1]=paste(moduleNumber,moduleColor,sep = "")
  for (j in 1:length(refpath)){
    fitest=ORA(testpath = moduleGenes,refpath = refpath[[j]],testbackground = allGenes,refbackground = refbackgroud[[j]])
    ORtable[i-1,j+1]=fitest[1]
    ptable[i-1,j+1]=fitest[2]
    print(length(unique(c(names(refpath)[j],names(refbackgroud)[j],colnames(ORtable)[j+1],colnames(ptable)[j+1]))))
  }
}

ptable1=as.data.frame(t(apply(ptable,1,function(x){p.adjust(x[2:ncol(ptable)],method = "fdr")})))
ptable1$Module=ptable$Module

library(tidyverse)
ptable2=ptable1 %>%
  gather("category","fdr",1:(ncol(ptable1)-1))%>%
  unite(id,Module,category,sep="++")

ORtable2=ORtable %>%
  gather("category","OR",2:ncol(ORtable))%>%
  unite(id,Module,category,sep="++")

mergeTable=merge(ORtable2,ptable2,by="id")
mergeTable1=mergeTable %>%
  separate(id,c("Module","category"),sep="\\+\\+")
mergeTable1$OR=as.numeric(mergeTable1$OR)
mergeTable1$signedLog10fdr = -log10(mergeTable1$fdr) * sign(log(mergeTable1$OR))
mergeTable1$signedLog10fdr[mergeTable1$OR==0 | mergeTable1$fdr>0.05]=0

mergeTable1$text = round(mergeTable1$OR, 1)
mergeTable1$text[mergeTable1$OR==0 | mergeTable1$fdr>0.05] = ""
mergeTable1$moudleColor=gsub("[0-9]+?","",mergeTable1$Module)
mergeTable1$moudleNumber=sapply(mergeTable1$Module,function(x)str_match(x,"^[0-9]*")[1,1])
mergeTable1$moudleNumber=as.numeric(mergeTable1$moudleNumber)
mergeTable1$Module1=reorder(mergeTable1$Module,mergeTable1$moudleNumber)
ggplot(mergeTable1, aes(x=Module1,y=category, label=text)) +
  geom_tile(aes(fill=signedLog10fdr),color="grey60") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red","signed\n-log10FDR\n") + 
  geom_text(size=3, color="black") + ylab("") + xlab("") + 
  theme(axis.text.x = element_text(angle=30, hjust=1), axis.text.y = element_text(size=7))+
  geom_tile(aes(x=Module,y=0.1),fill=mergeTable1$moudleColor,height=0.1)
save(mergeTable1,file = "working_data/listEnrich/fisherExact.neuronModule.houseList.RData")

# combine with permutation results
library(Cairo)
load("working_data/listEnrich/fisherExact.blockModule.houseList.RData")
load("working_data/listEnrich/permutTable.1000times.blockModule.houseList.RData")
lists=read.delim("F:/datasets/genesets/houseLists6.1moOrganoid.txt",header = F)
permTable1=permTable %>%
  unite(id,Module,category,sep="++") %>%
  dplyr::select(id,zscore,fdr,moduleColor,moduleNumber)

mergeTable2=mergeTable1 %>%
  unite(id,Module,category,sep="++") %>%
  dplyr::select(id,OR,fdr)

combineTable=merge(mergeTable2,permTable1,by="id")

combineTable1=combineTable %>% 
  separate(id,c("Module","category"),sep="\\+\\+") %>%
  mutate(fdr=pmax(fdr.x,fdr.y))

combineTable1$signedLog10fdr = -log10(combineTable1$fdr) * sign(combineTable1$zscore)
combineTable1$signedLog10fdr[combineTable1$fdr>0.05]=0

combineTable1$text = round(combineTable1$zscore, 1)
combineTable1$text[combineTable1$fdr>0.05] = ""
combineTable1$Module1=reorder(combineTable1$Module,combineTable1$moduleNumber)
# combineTable1$category1=factor(combineTable1$category,levels = lists$V1)
# combineTable1=combineTable1 %>% filter(category %in% c("yellow","green","blue","magenta","brown","tan"))
# combineTable1$category=factor(combineTable1$category,levels = c("yellow","green","blue","magenta","brown","tan"))
#combineTable1=combineTable1 %>% filter(category %in% lists$V1[1:11])
cairo_pdf("working_data/figures/listEnrich.blockModules.fisherExact.permut.showZscore.forTile.pdf",width = 17.8,height = 9.2,family = "Arial")
ggplot(combineTable1, aes(x=Module1,y=category, label=text)) +
  geom_tile(aes(fill=log2(signedLog10fdr+1)),color="grey60") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red","signed\n-log10FDR\n") + 
  geom_text(size=3, color="black") +
  ylab("") + 
  xlab("") + 
  theme(text=element_text(family = "Arial"),
        axis.text.x = element_text(angle=30, hjust=1), 
        axis.text.y = element_text(size=10))+
  geom_tile(aes(x=Module,y=0.1),fill=combineTable1$moduleColor,height=0.1)
dev.off()

cairo_pdf("working_data/figures/listEnrich.blockModules.fisherExact.permut.showZscore.forLegend.pdf",width = 17.8,height = 9.2,family = "Arial")
ggplot(combineTable1, aes(x=Module1,y=category, label=text)) +
  geom_tile(aes(fill=signedLog10fdr),color="grey60") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red","signed\n-log10FDR\n") + 
  geom_text(size=3, color="black") +
  ylab("") + 
  xlab("") + 
  theme(text=element_text(family = "Arial"),
        axis.text.x = element_text(angle=30, hjust=1), 
        axis.text.y = element_text(size=10))+
  geom_tile(aes(x=Module,y=0.1),fill=combineTable1$moduleColor,height=0.1)
dev.off()
