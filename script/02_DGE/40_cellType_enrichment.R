options(stringsAsFactors = F)
library(tidyverse)
library(EWCE)

rm(list = ls());gc()

load("working_data/dge/dge.neuron.fc2pass.RData")
dat=as.data.frame(dge.deseq2)
dat=dat[complete.cases(dat),]
rownames(dat)=sapply(rownames(dat),function(x) str_split_fixed(x,"\\.",Inf)[1,1])

annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")
dat$genename=annot$genename[match(rownames(dat),annot$geneid)]

bg=dat$genename
dge_up=dat$genename[dat$padj < 0.1 & dat$log2FoldChange > 0]
dge_down=dat$genename[dat$padj < 0.1 & dat$log2FoldChange < 0]

reps=10000
level=1
load("D:/projects/organoid/revision/CellTypeData_DamonNeuralFetalOnly.rda")
ctd_fetal = ctd
load("D:/projects/organoid/revision/CellTypeData_LakeNeuralAdultOnly.rda")
ctd_adult = ctd
rm(ctd)

up_adult_cells=bootstrap.enrichment.test(sct_data=ctd_adult,hits=dge_up,bg=bg,reps=reps,annotLevel=level,genelistSpecies="human",sctSpecies="human")
up_fetal_cells=bootstrap.enrichment.test(sct_data=ctd_fetal,hits=dge_up,bg=bg,reps=reps,annotLevel=level,genelistSpecies="human",sctSpecies="human")

down_adult_cells=bootstrap.enrichment.test(sct_data=ctd_adult,hits=dge_down,bg=bg,reps=reps,annotLevel=level,genelistSpecies="human",sctSpecies="human")
down_fetal_cells=bootstrap.enrichment.test(sct_data=ctd_fetal,hits=dge_down,bg=bg,reps=reps,annotLevel=level,genelistSpecies="human",sctSpecies="human")

up_adult=up_adult_cells$results
down_adult=down_adult_cells$results
up_fetal=up_fetal_cells$results
down_fetal=down_fetal_cells$results

up_adult$fdr=p.adjust(up_adult$p)
down_adult$fdr=p.adjust(down_adult$p)
up_fetal$fdr=p.adjust(up_fetal$p)
down_fetal$fdr=p.adjust(down_fetal$p)


up_adult$module = up_fetal$module ="DGE_UP"
down_adult$module = down_fetal$module ="DGE_DOWN"

up_adult$stage = down_adult$stage ="Adult"
up_fetal$stage = down_fetal$stage ="Fetal"


adult_dge=rbind(up_adult,down_adult)
fetal_dge=rbind(up_fetal,down_fetal)

dat = rbind(up_adult,down_adult,up_fetal,down_fetal)
dat$sd_from_mean[dat$sd_from_mean < 0]=0
dat$fdr[dat$fdr == 0]=min(dat$fdr[dat$fdr>0]) * 0.1
dat$sd_from_mean=round(dat$sd_from_mean,1)
dat$sd_from_mean[dat$fdr > 0.05]=""
dat$CellType=factor(dat$CellType,levels = dat$CellType[order(dat$stage)])
ggplot(dat,aes(CellType,module,label=sd_from_mean))+
  geom_tile(aes(fill=-log10(fdr)),color="grey60")+
  facet_wrap(~stage,scales = "free_x")+
  geom_text(size=3, color="black")+
  scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = 0.1)+
  xlab("")+ylab("")+
  theme(text=element_text(family="Arial",size=16),axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))

markers = read.delim("clipboard")
astro=markers$GeneName[grepl("Adult-In",markers$CellType)]
intersect(dge_up,markers$GeneName)
