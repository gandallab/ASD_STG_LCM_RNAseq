rm(list = ls());gc()
options(stringsAsFactors = F)
library(WGCNA)
library(tidyverse)
library(igraph)
library(pSI)
library(gridExtra)
library(Cairo)

# define functions
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
ORA <- function(testpath,refpath,testbackground,refbackground) {
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

# functions defined

load("working_data/wgcna/neuron.v2.moreCovar//neuron.bicor.recut.RData")
load("working_data/wgcna/voom.forWGCNA.input.neuron.v2.RData")
softpower=7

geneTree = networks$tree
datExpr=networks$datExpr
merged = networks$merged
modules = merged$colors
MEs = networks$MEs
kMEtable = networks$kMEtable
annot=read.delim("D:/references/gencode.v19.gene.name.txt",header = T,sep = "\t")
annot=annot[match(rownames(datExpr),annot$geneid),]
table(datMeta$Sample == colnames(datExpr))
datMeta$Diagnosis=factor(datMeta$Diagnosis,levels = c("Control","Autism"))
datMeta=datMeta[order(datMeta$Diagnosis),]
datMeta$seq=c(1:nrow(datMeta))
datMeta=datMeta[match(colnames(datExpr),datMeta$Sample),]
modTrait=read.delim("working_data/wgcna/neuron.v2.moreCovar//modTrait.neuron.bicor.txt")

#Load cell-type expression signatures for enrichment
pSI.zhang = read.csv(file="D:/datasets/psychencode/cellTypes/Zhang_Human_Neuro2015_pSI_GSE21653.csv",row.names=1,stringsAsFactors = F)

pSI.zeisel = read.csv(file="D:/datasets/psychencode/cellTypes/Zeisel_level1_Mouse_Science2015_pSI.ensg.csv",row.names=1,stringsAsFactors = F)

pSI.goldmann = read.csv(file="D:/datasets/psychencode/cellTypes/Goldman_levelHybrid_Mouse_NatImmunol2015_pSI.ensg.csv",row.names=1,stringsAsFactors = F)

markers.Habib = read.csv("D:/datasets/psychencode/cellTypes/Habib_NatMeth2017_Human_DroncSeq_clusters.csv",row.names=1)

# loop through each module except gray

for(i in 2:ncol(MEs$eigengenes)) {
  
  me = MEs$eigengenes[,i]
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  moduleColor = labels2colors(moduleNumber)
  moduleGenes = annot$geneid[modules==moduleNumber]
  
  #cairo_pdf(file=paste0("working_data//figures/moduleAnnot.neuron/", moduleNumber, moduleColor, ".pdf"), width=21, height=11,family = "Arial")
  
  s = summary(lm(me ~ Diagnosis + Age + Sex,data=datMeta))$coefficients
  #  s = tryCatch(summary(lme(me ~ Group, data=datMeta, random = ~1 | Individual))$tTable,error=function(e){NA})
  dat2=data.frame()
  if (!is.na(s)){
    for(grp in c("Autism")) {
      rowID = paste0("Diagnosis", grp)
      dat2 = rbind(dat2,
                   data.frame(Module=moduleColor, moduleNumber= moduleNumber, Group=grp,
                              beta = s[rowID, "Estimate"], SE = s[rowID, "Std. Error"], t=s[rowID, "t value"], p=s[rowID, "Pr(>|t|)"]))
      # dat2 = rbind(dat2,
      #              data.frame(Module=moduleColor, moduleNumber= moduleNumber, Group=grp,
      #                         beta = s[rowID, "Value"], SE = s[rowID, "Std.Error"], t=s[rowID, "t-value"], p=s[rowID, "p-value"]))
    }
    currentModTrait=modTrait %>%
      filter(moduleNumber== !!moduleNumber)
    
    dat2$fdr = currentModTrait$fdr[match(dat2$Group,currentModTrait$Group)]
    
    dat2$p.symbol = ""; dat2$p.symbol[dat2$fdr<.05] ="*"
    g1=ggplot(dat2, aes(x=Group, y=beta, label=p.symbol)) + geom_bar(stat="identity",fill=moduleColor) +  
      geom_errorbar(aes(ymin=(beta - SE), ymax=(beta + SE)), position=position_dodge(width=0.8), width=0.25,size=0.25) +
      geom_text(color="red",size=8,aes(y=beta+ 1.3*sign(beta)*SE ))+
      ggtitle("Module-Trait Association")+
      xlab("")+
      ylab("Linear Regression Beta")+
      theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size=15),axis.title.y = element_text(size=20))
  }else{
    df <- data.frame()
    g1=ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10)
  }
  
  
  
  #Module eigengene vs Age trajectory
  dat=data.frame(Eigengene=me, Age=datMeta$Age,Group=datMeta$Diagnosis)
  g2=ggplot(dat,aes(x=Age,y=Eigengene,color=Group)) + geom_point(size=3)+ geom_smooth(method = "lm")+
    xlab("Age")+ylab("Eigengene")+ggtitle("Eigengene Expression")+
    #    scale_color_manual(values = c("black","red","blue"))+
    theme(axis.text.x = element_blank(),axis.text.y = element_text(size=15),axis.title.y = element_text(size=20),legend.text = element_text(size = 15),
          legend.title = element_text(size=17))
  
  #Gene ontology  
  go = read.delim(paste("working_data/wgcna/neuron.v2.moreCovar//go//go.M",moduleNumber,"_",moduleColor,".txt",sep = ""),header = T,sep = "\t")
  go=go[complete.cases(go),]
  go=go[go$overlap.size > 1,]
  go=go[1:min(10,nrow(go)),]
  if(nrow(go)==0){
    df <- data.frame()
    g3=ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10)
  }else{
    g3=ggplot(go, aes(x=reorder(term.name, -log10(p.value)), y=-log10(p.value))) + geom_bar(stat="identity", fill="royalblue") + 
      coord_flip() + xlab("") + geom_hline(yintercept=-log10(0.05), lty=2, color="red")+
      theme(axis.text.y = element_text(size=12,face = "bold"),axis.text.x = element_text(size = 14),axis.title.x = element_text(size=20))+
      ggtitle("Functional Enrichment")
  }
  
  #Hub gene network
  hubGenes = moduleGenes[order(kMEtable[moduleGenes,i], decreasing = T)[1:25]]
  hubGene.symbols = annot$genename[match(hubGenes, annot$geneid)]
  adjMat = adjacency(t(datExpr[hubGenes,]),type = "signed",corFnc = "bicor", power=softpower)
  adjMat[adjMat < quantile(adjMat,0.1)]=0
  graph <- graph.adjacency(as.matrix(adjMat),mode="undirected",weighted=T,diag=F)
  plotcord= data.frame(layout_with_fr(graph))
  colnames(plotcord) = c("X1","X2")
  edgelist <- get.edgelist(graph,names = F)
  edges <- data.frame(plotcord[edgelist[,1],], plotcord[edgelist[,2],])
  colnames(edges) <- c("X1","Y1","X2","Y2")
  plotcord = cbind(plotcord, data.frame(gene=hubGene.symbols))
  
  g4=ggplot() + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.5, colour="grey") + 
    geom_point(aes(X1, X2), data=plotcord,color=moduleColor,size=4) + geom_text(aes(x=X1, y=X2+.2, label=gene),fontface="bold",size=6,data=plotcord) +
    theme_classic() + labs(x="", y="") + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())+
    ggtitle("Hub Genes Network")
  
  
  #Cell-Type enrichment
  f.zhang = fisher.iteration(pSI.zhang, moduleGenes,p.adjust = T)
  
  f.zeisel = fisher.iteration(pSI.zeisel, moduleGenes,p.adjust = T)
  
  f.goldman = fisher.iteration(pSI.goldmann, moduleGenes,p.adjust = T)
  
  
  
  f.habib = data.frame(P=matrix(NA,nrow=length(unique(na.omit(markers.Habib$Cluster.ID))), ncol=1), row.names=unique(na.omit(markers.Habib$Cluster.ID)))
  
  for(j in 1:nrow(f.habib)) {
    
    f = ORA(testpath = annot$genename[modules==moduleNumber], 
            
            refpath = markers.Habib$Gene.ID[na.omit(markers.Habib$Cluster.ID==rownames(f.habib)[j])],
            
            testbackground = annot$genename, refbackground = markers.Habib$Gene.ID)
    
    if(as.numeric(f[[1]]) > 1) { 
      
      f.habib$P[j] = as.numeric(f[[2]])
      
    } else {
      
      f.habib$P[j] = 1
      
    }
    
  }
  modCellType = rbind(data.frame(Dataset="Zhang", CellType=rownames(f.zhang), log10fdr=-log10(f.zhang[,2])),
                      
                      data.frame(Dataset="Zeisel", CellType=rownames(f.zeisel), log10fdr=-log10(f.zeisel[,2])),
                      
                      data.frame(Dataset="Goldman", CellType=rownames(f.goldman), log10fdr=-log10(f.goldman[,2])),
                      
                      data.frame(Dataset="Habib",  CellType=rownames(f.habib), log10fdr=-log10(f.habib[,1])))
  
  
  
  modCellType$CellType = factor(modCellType$CellType, levels=unique(modCellType$CellType))
  g5=ggplot(modCellType,aes(x=CellType,y=log10fdr, fill=Dataset)) + geom_bar(stat="identity") + coord_flip() + geom_hline(yintercept = -log10(0.05),lty=2) + xlab("") + 
    ggtitle("Cell-Type Enrichment") +theme(axis.text.y = element_text(size=14,face = "bold"),axis.text.x = element_text(size=9),axis.title.x = element_text(size=15))
  
  
  
  
  grid.arrange(grobs=list(g1,g2,g3,g4,g5), layout_matrix=rbind(c(1,2,2,2,5,5,5),c(3,3,3,4,4,4,4)),
               top=paste0("Module ", moduleNumber, " (", moduleColor, ")"),padding=unit(2, "line"))
  if (i != 2){
    ind=TRUE
  }else{
    ind=FALSE
  }
  graph2ppt(file="working_data/figures/module.annotation.neuron.v2.pptx",width=18,height=15,append=ind)
  #  dev.off()
}


