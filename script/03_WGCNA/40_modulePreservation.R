rm(list=ls());gc()
options(stringsAsFactors = F)
library(WGCNA)

# blocks
load("working_data/wgcna/blocks.fc/blocks.fc.recut.RData")
blocks=networks;rm(networks)
datExpr.blocks=blocks$datExpr
merged.blocks = blocks$merged
modules.blocks = merged.blocks$colors
length(unique(modules.blocks))
MEs.blocks = blocks$MEs
kMEtable.blocks = blocks$kMEtable

MEblocks=MEs.blocks$eigengenes[,-1]
MEblocks1=hclust(dist(t(MEblocks)))
colors.blocks=labels2colors(as.numeric(gsub("ME","",MEblocks1$labels[MEblocks1$order])))

member.blocks=labels2colors(modules.blocks)
membership.blocks=data.frame(id=rownames(datExpr.blocks),module=member.blocks)

# neuron
load("working_data/wgcna/neuron.fc/neuron.fc.manualRemove.recut.RData")
neuron=networks;rm(networks)
datExpr.neuron=neuron$datExpr
merged.neuron = neuron$merged
modules.neuron = merged.neuron$colors
length(unique(modules.neuron))
MEs.neuron = neuron$MEs
kMEtable.neuron = neuron$kMEtable

MEneuron=MEs.neuron$eigengenes[,-1]
MEneuron1=hclust(dist(t(MEneuron)))
colors.neuron=labels2colors(as.numeric(gsub("ME","",MEneuron1$labels[MEneuron1$order])))

member.neuron=labels2colors(modules.neuron)
membership.neuron=data.frame(id=rownames(datExpr.neuron),module=member.neuron)

# oligo
load("working_data/wgcna/oligo.fc/oligo.fc.recut.RData")
oligo=networks;rm(networks)
datExpr.oligo=oligo$datExpr
merged.oligo = oligo$merged
modules.oligo = merged.oligo$colors
length(unique(modules.oligo))
MEs.oligo = oligo$MEs
kMEtable.oligo = oligo$kMEtable

MEoligo=MEs.oligo$eigengenes[,-1]
MEoligo1=hclust(dist(t(MEoligo)))
colors.oligo=labels2colors(as.numeric(gsub("ME","",MEoligo1$labels[MEoligo1$order])))

member.oligo=labels2colors(modules.oligo)
membership.oligo=data.frame(id=rownames(datExpr.oligo),module=member.oligo)


# perservation
multiExpr=list(blocks=list(data=t(datExpr.blocks)),neuron=list(data=t(datExpr.neuron)),oligo=list(data=t(datExpr.oligo)))
multiColor=list(blocks=membership.blocks$module,neuron=membership.neuron$module,oligo=membership.oligo$module)
save(file = "working_data/wgcna/preservation.fc/mp.input.RData",multiColor,multiExpr)

load("fetalBrain/wgcna/preservation/mp.input.RData")
mp = modulePreservation(multiExpr, multiColor,
                        networkType = "signed",
                        corFnc = "bicor",
                        referenceNetworks = c(1,2),
                        nPermutations = 100,
                        randomSeed = 1,
                        quickCor = 0,
                        useInterpolation=F,
                        verbose = 3)
save(file = "./preservation/mp.output.3asRef.RData",mp)

# plot
rm(list=ls());gc()
options(stringsAsFactors = F)
library(tidyverse)
library(ggrepel)
library(WGCNA)

colTrans=data.frame(number=seq(1,100),color=labels2colors(seq(1,100)))

load("working_data/wgcna/preservation.fc/mp.output.RData")
ref = 3
test = 2
plotData=cbind(mp$preservation$observed[[ref]][[test]][, 1:2], mp$preservation$Z[[ref]][[test]][, 1:2])
plotData[,3]=NULL
plotData$color=rownames(plotData)
plotData1=plotData[!plotData$color %in% c("grey","gold"),]

cairo_pdf("../modulePreservations/modulePreservation.zsummary.3moRef.1moTest.pdf",width = 9.6,height = 9.1,family = "Arial")
ggplot(plotData1,aes(moduleSize,Zsummary.pres))+
  geom_point(size=6,color=plotData1$color)+
  scale_y_continuous(breaks = pretty(plotData1$Zsummary.pres,n=7))+
  scale_x_continuous(trans = "log1p",breaks = c(50,100,200,300,500,700,900))+
  xlab("Module size")+
  ylab("Preservation Zsummary")+
  geom_text_repel(aes(moduleSize,Zsummary.pres,label=color))+
  geom_hline(yintercept = c(2,10),linetype=2,color=c("blue","green"))+
  theme_classic()+
  theme(axis.title=element_text(size=17,face = "bold"),
        axis.text = element_text(face = "bold",size = 16),
        text = element_text(family = "Arial"))

dev.off()

counts=as.data.frame(mp$accuracy$observedCounts[[ref]][[test]])
counts$testModule=rownames(counts)
counts1=counts %>%
  filter(testModule != "grey") %>%
  select(-grey) %>%
  gather("refModule","overlap",1:(ncol(counts)-2)) %>%
  unite("modulePair",testModule,refModule)

fishp=as.data.frame(mp$accuracy$observedFisherPvalues[[ref]][[test]])
fishp$testModule=rownames(fishp)
fishp1=fishp %>%
  filter(testModule != "grey") %>%
  select(-grey) %>%
  gather("refModule","fishp",1:(ncol(counts)-2)) %>%
  unite("modulePair",testModule,refModule)

plotData2=merge(counts1,fishp1,by="modulePair") %>%
  separate(modulePair,c("testModule","refModule"),sep="_")
plotData2$overlap[plotData2$fishp > 0.01]=""
plotData2$testModule=reorder(plotData2$testModule,colTrans$number[match(plotData2$testModule,colTrans$color)])
plotData2$refModule=reorder(plotData2$refModule,colTrans$number[match(plotData2$refModule,colTrans$color)])
# ggplot(plotData2,aes(x=rnaModule,y=proteinModule,label=overlap))+
#   geom_tile(aes(fill=log2(-log10(fishp)+1)))+
#   scale_fill_gradient2(low = "white",mid = "white", high = "red",midpoint = 1.58)+
#   geom_text(size=3,color="black")

# plotData3=plotData2 %>%
#   filter(rnaModule %in% c("violet","sienna3","mediumpurple3","bisque4","skyblue2","plum","greenyellow","darkgrey","lightcyan1","lightpink4","coral1"))%>%
#   filter(proteinModule %in% c("blue","red","magenta","tan","cyan","grey60","lightyellow","turquoise","brown"))
cairo_pdf("../modulePreservations/modulePreservation.crossTab.1moRef.3moTest.pdf",width = 9.6,height = 10.3,family = "Arial")
ggplot(plotData2,aes(x=testModule,y=refModule,label=overlap))+
  geom_tile(aes(fill=-log10(fishp)),color="grey90")+
  scale_fill_gradient2(low = "white",mid = "white", high = "red",midpoint = 2,"-log10\np-value")+
  geom_text(size=3,color="black")+
  xlab("test module")+
  ylab("ref module")+
  #coord_flip()+
  theme(axis.text.x = element_text(angle=30, hjust=1), 
        axis.text.y = element_text(size=10),
        axis.title = element_text(size = 12),
        text = element_text(family = "Arial",face = "bold"))+
  geom_tile(aes(x=testModule,y=0.1),fill=plotData2$testModule,height=0.2)+
  geom_tile(aes(x=0.1,y=refModule),fill=plotData2$refModule,width=0.3)
dev.off()

