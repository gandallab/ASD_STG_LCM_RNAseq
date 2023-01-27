rm(list = ls());gc()
options(stringsAsFactors = F)
library(WGCNA)
library(tidyverse)
library(igraph)
library(pSI)
library(gridExtra)
library(Cairo)


load("working_data/wgcna/blocks.fc/blocks.fc.recut.RData")
softpower=8

geneTree = networks$tree
datExpr=networks$datExpr
datMeta=networks$datMeta
merged = networks$merged
modules = merged$colors
MEs = networks$MEs
kMEtable = networks$kMEtable
rownames(datExpr)=sapply(rownames(datExpr),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
rownames(kMEtable)=sapply(rownames(kMEtable),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt",header = T,sep = "\t")
annot=annot[match(rownames(datExpr),annot$geneid),]
table(datMeta$Sample == colnames(datExpr))


i=11

me = MEs$eigengenes[,i]
moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
moduleColor = labels2colors(moduleNumber)
#if (!moduleColor %in% modTrait$Module[modTrait$fdr < 0.05]){next}
moduleGenes = annot$geneid[modules==moduleNumber]
moduleGenes_name=annot$genename[modules==moduleNumber]

gmt=read.gmt("D:/datasets/genesets/houseLists03242020.gmt")
genename=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt",header = T,sep = "\t")
gmt$geneid=genename$geneid[match(gmt$gene,genename$genename)]

asd102=gmt$gene[gmt$ont == "ASD_102genes"]
syndromic=gmt$gene[gmt$ont == "Syndromic ASD"]
length(intersect(asd102,syndromic))

pli099=gmt$gene[gmt$ont == "pLI>0.99"]
constrained=gmt$gene[gmt$ont == "Constrained"]
length(intersect(pli099,constrained))

pres=gmt$gene[gmt$ont == "Presynaptic"]
posts=gmt$gene[gmt$ont == "Postsynaptic"]

synaptic=union(pres,posts)

intersect(synaptic, moduleGenes_name)
intersect(syndromic, moduleGenes_name)
intersect(asd102, moduleGenes_name)

filter1=intersect(synaptic, moduleGenes_name)
gmt1=gmt[gmt$gene %in% filter1,]
filter2=gmt1$gene[gmt1$ont == "pLI>0.99"]

final=union(filter2,intersect(syndromic, moduleGenes_name))
node_property=data.frame(node=final,class="LOF intolerant synaptic genes")
node_property$class[node_property$node %in% syndromic]="ASD genes"
write.table(node_property,file = "working_data/wgcna/blocks.fc/cytoscape.nodeProperty.txt",quote = F,sep = "\t",row.names = F)


#Hub gene network
hubGene.symbols = final
hubGenes=annot$geneid[match(hubGene.symbols,annot$genename)]

adjMat = adjacency(t(datExpr[hubGenes,]),type = "signed",corFnc = "bicor", power=softpower)
adjMat[adjMat < quantile(adjMat,0.1)]=0
graph <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=F)
el=cbind(get.edgelist(graph),round(E(graph)$weight,3))
el1=data.frame(V1=annot$genename[match(el[,1],annot$geneid)],
               V2=annot$genename[match(el[,2],annot$geneid)],
               weight=as.numeric(el[,3]))
el2=el1[order(el1$weight,decreasing = T),]
for (i in 1:nrow(el2)){
  current=el2[1:i,1:2]
  print(c(i,length(unique(c(current$V1,current$V2)))))
}
write.table(el2[1:233,],file = "working_data/wgcna/blocks.fc/cytoscape.edgelist.v2.txt",quote = F,sep = "\t",row.names = F)






