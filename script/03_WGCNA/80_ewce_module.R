rm(list = ls());gc()
options(stringsAsFactors = F)
library(WGCNA)
library(tidyverse)
library(EWCE)


#load("working_data/wgcna/blocks.fc/blocks.fc.recut.RData")
load("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/neuron.fc2pass.recut.RData")

geneTree = networks$tree
datExpr=networks$datExpr
rownames(datExpr)=sapply(rownames(datExpr),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
merged = networks$merged
modules = merged$colors
MEs = networks$MEs
kMEtable = networks$kMEtable

annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt",header = T,sep = "\t")
annot=annot[match(rownames(datExpr),annot$geneid),]


#modTrait=read.delim("./sv1.6/modTrait.ds2mod70.txt")

#Load cell-type expression signatures for enrichment
reps=10000
level=2
load("D:/projects/organoid/revision/CellTypeData_DamonNeuralFetalOnly.rda")
ctd_fetal = ctd
load("D:/projects/organoid/revision/CellTypeData_LakeNeuralAdultOnly.rda")
ctd_adult = ctd
rm(ctd)
bg=annot$genename
background_adult = intersect(bg,rownames(ctd_adult[[level]]$specificity))
background_fetal = intersect(bg,rownames(ctd_fetal[[level]]$specificity))

# loop through each module except gray
adult_all=data.frame()
fetal_all=data.frame()
for(i in 2:ncol(MEs$eigengenes)) {
  
  me = MEs$eigengenes[,i]
  moduleNumber = as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
  #if (!moduleNumber %in% c(5,8,10,20,30)) {next}
  moduleColor = labels2colors(moduleNumber)
  moduleGenes = annot$genename[modules==moduleNumber]
  
  input_adult = intersect(moduleGenes,rownames(ctd_adult[[level]]$specificity))
  input_fetal = intersect(moduleGenes,rownames(ctd_fetal[[level]]$specificity))
  
  
  #Cell-Type enrichment
  adult_cells=bootstrap_enrichment_test(sct_data=ctd_adult,hits=input_adult,bg=background_adult,reps=reps,annotLevel=level,genelistSpecies="human",sctSpecies="human")
  fetal_cells=bootstrap_enrichment_test(sct_data=ctd_fetal,hits=input_fetal,bg=background_fetal,reps=reps,annotLevel=level,genelistSpecies="human",sctSpecies="human")
  
  # print(adult_cells$results[order(adult_cells$results$p),3:5][1:6,])
  # print(ewce.plot(adult_cells$results,mtc_method="BH")$plain)
  # print(fetal_cells$results[order(fetal_cells$results$p),3:5][1:6,])
  # print(ewce.plot(fetal_cells$results,mtc_method="BH")$plain)
  
  adult_results=adult_cells$results
  adult_results$fdr=p.adjust(adult_results$p)
  adult_results$moduleNumber=moduleNumber
  adult_results$moduleColor=moduleColor
  fetal_results=fetal_cells$results
  fetal_results$fdr=p.adjust(fetal_results$p)
  fetal_results$moduleNumber=moduleNumber
  fetal_results$moduleColor=moduleColor
  
  adult_all=rbind(adult_all,adult_results)
  fetal_all=rbind(fetal_all,fetal_results)
  
}

save(adult_all,fetal_all,file = "working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/ewce.neuron.level1.RData")

tmp=adult_all[grepl("In",adult_all$CellType),]

# plot adult

flag=unique(adult_all$module[adult_all$fdr < 0.05])
dat=adult_all[adult_all$module %in% flag,]
#dat$realp=2*pnorm(-abs(dat$sd_from_mean))
dat$sd_from_mean[dat$sd_from_mean < 0]=0
dat$fdr[dat$fdr == 0]=min(dat$fdr[dat$fdr>0]) * 0.1
dat$sd_from_mean=round(dat$sd_from_mean,1)
dat$sd_from_mean[dat$fdr > 0.05]=""
dat$moduleNumber=as.numeric(str_sub(str_extract(dat$module,"_[1-9]*"),2,-1))
dat$module=reorder(dat$module,dat$moduleNumber)
ggplot(dat,aes(CellType,module,label=sd_from_mean))+
  geom_tile(aes(fill=-log10(fdr)))+
  geom_text(size=3, color="black")+
  scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = 1.3)

# plot fetal

flag=unique(fetal_all$module[fetal_all$fdr < 0.05])
dat=fetal_all[fetal_all$module %in% flag,]
#dat$realp=2*pnorm(-abs(dat$sd_from_mean))
dat$sd_from_mean[dat$sd_from_mean < 0]=0
dat$fdr[dat$fdr == 0]=min(dat$fdr[dat$fdr>0]) * 0.1
dat$sd_from_mean=round(dat$sd_from_mean,1)
dat$sd_from_mean[dat$fdr > 0.05]=""
dat$moduleNumber=as.numeric(str_sub(str_extract(dat$module,"_[0-9]*"),2,-1))
dat$module=reorder(dat$module,dat$moduleNumber)
ggplot(dat,aes(CellType,module,label=sd_from_mean))+
  geom_tile(aes(fill=-log10(fdr)))+
  geom_text(size=3, color="black")+
  scale_fill_gradient2(low = "white",mid = "white",high = "red",midpoint = 1.3)
