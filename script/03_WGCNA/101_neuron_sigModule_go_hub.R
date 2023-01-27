rm(list = ls());gc()
options(stringsAsFactors = F)
library(WGCNA)
library(tidyverse)
library(igraph)
library(Cairo)


# load PPI
ppi=read.delim("D:/datasets/ppi/allVidal.lit17.mvp.hprdInvivo.hintBinary.txt")

load("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/neuron.fc2pass.recut.RData")
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
datMeta$Diagnosis=factor(datMeta$Diagnosis,levels = c("Control","Autism"))
datMeta=datMeta[order(datMeta$Diagnosis),]
datMeta$seq=c(1:nrow(datMeta))
datMeta=datMeta[match(colnames(datExpr),datMeta$Sample),]
modTrait=read.delim("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/modTrait.neuron.fc2pass.txt")

# loop through each module except gray

for(i in 2:ncol(MEs$eigengenes)) {
  
  me = MEs$eigengenes[,i]
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  moduleColor = labels2colors(moduleNumber)
  if (!moduleColor %in% modTrait$Module[modTrait$fdr < 0.05]){next}
  moduleGenes = annot$geneid[modules==moduleNumber]
  
  
  
  #Gene ontology  
  go = read.delim(paste("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100//go//go.M",moduleNumber,"_",moduleColor,".txt",sep = ""),header = T,sep = "\t")
  go=go[complete.cases(go),]
  go=go[go$domain != "keg",]
  go=go[1:min(7,nrow(go)),]
  if(nrow(go)>0){
    cairo_pdf(file=paste0("working_data//figures/moduleAnnot.neuron.fc2pass/", moduleNumber, moduleColor, ".go.pdf"), width=9.7, height=3.8,family = "Arial")
    
    p=ggplot(go, aes(x=reorder(term.name, -log10(p.value)), y=-log10(p.value))) + geom_bar(stat="identity", fill="grey",width = 0.3) + 
      coord_flip() + xlab("") + ylab("-log10(FDR)")+
      geom_hline(yintercept=-log10(0.05), lty=2, color="red")+
      theme_bw()+
      theme(text=element_text(family = "Arial",face = "bold"),
            axis.text.y = element_text(size=18),
            axis.text.x = element_text(size = 18),
            axis.title.x = element_text(size=20))
    print(p)
    dev.off()
  }
  
  #Hub gene network
  hubGenes = moduleGenes[order(kMEtable[moduleGenes,i], decreasing = T)[1:50]]
  hubGene.symbols = annot$genename[match(hubGenes, annot$geneid)]
  adjMat = adjacency(t(datExpr[hubGenes,]),type = "signed",corFnc = "bicor", power=softpower)
  #adjMat[adjMat < quantile(adjMat,0.1)]=0
  graph <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=F)
  plotcord= data.frame(layout_with_fr(graph))
  colnames(plotcord) = c("X1","X2")
  edgelist <- get.edgelist(graph,names = F)
  edgelist1=as.data.frame(edgelist)
  
  edgelist1$V3=apply(edgelist1,1,function(x){
    a=c(hubGenes[x[1]],hubGenes[x[2]])
    b=a[order(a)]
    paste(b[1],b[2],sep = "_")
  })
  edgelist1$V4="N"
  edgelist1$V4[edgelist1$V3 %in% ppi$V3]="Y"
  table(edgelist1$V4)
  edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
  colnames(edges) <- c("X1","Y1","X2","Y2")
  edges$ppi=edgelist1$V4
  plotcord = cbind(plotcord, data.frame(gene=hubGene.symbols))
  
  cairo_pdf(file=paste0("working_data//figures/moduleAnnot.neuron.fc2pass/", moduleNumber, moduleColor, ".hub.pdf"), width=9.7, height=8,family = "Arial")
  
  p=ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.5, colour="grey",alpha=0.3) +
    geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges[edges$ppi == "Y",], size = 1, colour="black") +
    geom_point(aes(X1, X2), data=plotcord,color=moduleColor,size=4) + geom_text(aes(x=X1, y=X2+.2, label=gene),fontface="bold",size=6,data=plotcord) +
    theme_classic() + labs(x="", y="") + 
    theme(text=element_text(family = "Arial",face = "bold",size=12),
          axis.line = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  print(p)
  dev.off()
}