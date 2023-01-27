rm(list = ls());gc()
options(stringsAsFactors = F)
library(WGCNA)
library(tidyverse)
library(igraph)
library(pSI)
library(gridExtra)
library(Cairo)



# load blocks data first
load("working_data/wgcna/cross_datasets_me.fc.RData")
## blocks data loaded

# load cell type enrichment
load("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/ewce.neuron.level2.RData")

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
  

  cairo_pdf(file=paste0("working_data//figures/moduleAnnot.neuron.fc2pass/", moduleNumber, moduleColor, ".pdf"), width=10, height=10,family = "Arial")
  
  # # neruon MEs
  # 
  # 
  # dat2=rbind(modTrait[modTrait$moduleNumber == moduleNumber , c(1:8)],modTrait.exprBlocks.moduleNeuron[modTrait.exprBlocks.moduleNeuron$moduleNumber == moduleNumber , ])
  # dat2$Group=c("Neuron","Blocks")
  # 
  # dat2$p.symbol = ""; dat2$p.symbol[dat2$fdr<.05] ="*"
  # g1=ggplot(dat2, aes(x=Group, y=beta, label=p.symbol)) + geom_bar(stat="identity",fill=moduleColor,width = 0.6) +  
  #   geom_errorbar(aes(ymin=(beta - SE), ymax=(beta + SE)), position=position_dodge(width=0.8), width=0.25,size=0.25) +
  #   geom_text(color="red",size=8,aes(y=beta+ 1.3*sign(beta)*SE ))+
  #   ggtitle("Module-Trait Association")+
  #   xlab("")+
  #   ylab("Linear Regression Beta")+
  #   theme_classic(base_family = "Arial")+
  #   theme(panel.background = element_rect(colour = "black"),axis.text.x = element_text(size = 18),axis.text.y = element_text(size=18),axis.title.y = element_text(size=20))
  # 
  # ggplot(final, aes(x=Group, y=beta, label=p.symbol,fill=Module)) + geom_bar(stat="identity",width = 0.6,position = "dodge") +  
  #   geom_errorbar(aes(ymin=(beta - SE), ymax=(beta + SE)), position=position_dodge(width=0.5), width=0.25,size=0.25) +
  #   geom_text(color="red",size=8,aes(y=beta+ 1.3*sign(beta)*SE ))+
  #   scale_fill_manual(values = c("purple"="purple","royalblue"="royalblue"))+
  #   ggtitle("Module-Trait Association")+
  #   xlab("")+
  #   ylab("Linear Regression Beta")+
  #   theme_classic(base_family = "Arial")+
  #   theme(panel.background = element_rect(colour = "black"),axis.text.x = element_text(size = 18),axis.text.y = element_text(size=18),axis.title.y = element_text(size=20))
  # 
  # 
  
  
  # #Module eigengene vs Age trajectory
  # dat=data.frame(Eigengene=me, Age=datMeta$Age,Group=datMeta$Diagnosis)
  # g2=ggplot(dat,aes(x=Age,y=Eigengene,color=Group)) + geom_point(size=3)+ geom_smooth(method = "lm")+
  #   xlab("Age")+ylab("Eigengene")+ggtitle("Eigengene Expression")+
  #   #    scale_color_manual(values = c("black","red","blue"))+
  #   theme(axis.text = element_text(size=15),axis.title.y = element_text(size=20),legend.text = element_text(size = 15),
  #         legend.title = element_text(size=17))
  
  #Gene ontology  
  go = read.delim(paste("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100//go//go.M",moduleNumber,"_",moduleColor,".txt",sep = ""),header = T,sep = "\t")
  go=go[complete.cases(go),]
  go=go[go$domain != "keg",]
  #go=go[go$overlap.size > 1,]
  #go=go[go$NES > 0,]
  go=go[1:min(7,nrow(go)),]
  #go$term.name=gsub("leukocyte cytokine","leukocyte \n cytokine",go$term.name)
  if(nrow(go)==0){
    df <- data.frame()
    g3=ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10)
  }else{
    g3=ggplot(go, aes(x=reorder(term.name, -log10(p.value)), y=-log10(p.value))) + geom_bar(stat="identity", fill="grey",width = 0.3) + 
      coord_flip() + xlab("") + ylab("-log10(FDR)")+
      geom_hline(yintercept=-log10(0.05), lty=2, color="red")+
      theme_classic(base_family = "Arial")+
      theme(axis.text.y = element_text(size=16,face = "bold"),axis.text.x = element_text(size = 18),axis.title.x = element_text(size=20),panel.background = element_rect(colour = "black"))+
      ggtitle("Functional Enrichment")
  }
  
  #Hub gene network
  # allpairs=as.data.frame(t(combn(moduleGenes,2)))
  # allpairs$V3=apply(allpairs,1,function(x){
  #   a=c(x[1],x[2])
  #   b=a[order(a)]
  #   paste(b[1],b[2],sep = "_")
  # })
  # 
  # ppinode=c()
  # for (eachpair in intersect(allpairs$V3,ppi$V3)){
  #   ppinode=c(ppinode,str_split_fixed(eachpair,"_",Inf)[1,])
  # }
  # ppinode=unique(ppinode)
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
  
  g4=ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.5, colour="grey",alpha=0.3) +
    geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges[edges$ppi == "Y",], size = 1, colour="black") +
    geom_point(aes(X1, X2), data=plotcord,color=moduleColor,size=4) + geom_text(aes(x=X1, y=X2+.2, label=gene),fontface="bold",size=6,data=plotcord) +
    theme_classic() + labs(x="", y="") + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())+
    ggtitle("Hub Genes Network")
  
  
  #Cell-Type enrichment
  
  ctdat=adult_all[adult_all$moduleNumber == moduleNumber,]
  ctdat=ctdat[grepl("Ex|In",ctdat$CellType),]
  ctdat$neurontype=str_sub(ctdat$CellType,1,2)
  ctdat$sd_from_mean[ctdat$sd_from_mean < 0]=0.01
  ctdat$symbol=""
  ctdat$symbol[ctdat$fdr < 0.05]="*"
  g5=ggplot(ctdat)+
    geom_col(aes(CellType,sd_from_mean,fill=neurontype),width = 0.6)+
    geom_text(aes(CellType,sd_from_mean + 0.2,label=symbol),color="red",size = 8)+
    xlab("")+ylab("Z-score")+
    theme_classic(base_family = "Arial")+
    theme(panel.background = element_rect(colour = "black"),
          axis.text = element_text(size=18),axis.title = element_text(size=20),
          axis.text.x = element_text(angle = 90,hjust = 0.5,vjust = 0.5))+
    ggtitle("Cell-type enrichment")
  

  
  
  
  grid.arrange(grobs=list(g3,g4), layout_matrix=rbind(c(1),c(1),c(2),c(2),c(2),c(2),c(2)))
  # if (i != 2){
  #   ind=TRUE
  # }else{
  #     ind=FALSE
  #     }
  # graph2ppt(file="working_data/figures/module.annotation.blocks.v2.moreCovar.pptx",width=18,height=15,append=ind)
  dev.off()
}

joint_ct=adult_all[adult_all$moduleColor %in% modTrait$Module[modTrait$fdr < 0.05],]
joint_ct=joint_ct[grepl("Ex|In",joint_ct$CellType),]
joint_ct$neurontype=str_sub(joint_ct$CellType,1,2)
joint_ct$sd_from_mean[joint_ct$sd_from_mean < 0]=0.01
joint_ct$symbol=""
joint_ct$symbol[joint_ct$fdr < 0.05]="*"
ggplot(joint_ct)+
  geom_col(aes(CellType,sd_from_mean,fill=moduleColor),width = 1,position = "dodge",color="black")+
  geom_text(aes(CellType,sd_from_mean + 0.2,label=symbol),color="red",size = 8)+
  scale_fill_manual(values = c("green"="green","red"="red","greenyellow"="greenyellow","lightcyan"="lightcyan","grey60"="grey60"))+
  facet_wrap(~neurontype,scales = "free_x")+
  xlab("")+ylab("Z-score")+
  theme_classic(base_family = "Arial")+
  theme(panel.background = element_rect(colour = "black"),
        axis.text = element_text(size=18),axis.title = element_text(size=20),
        axis.text.x = element_text(angle = 90,hjust = 0.5,vjust = 0.5))+
  ggtitle("Cell-type enrichment")
