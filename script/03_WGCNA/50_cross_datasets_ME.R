rm(list = ls());gc()
options(stringsAsFactors = F)
library(WGCNA)

load("working_data/wgcna/blocks.fc/blocks.fc.recut.RData")
blocks=networks;rm(networks)
blocks.datExpr=blocks$datExpr
blocks.datMeta=blocks$datMeta
blocks.modules=blocks$merged$colors
blocks.genes=rownames(blocks.datExpr)
blocks.power=8

load("working_data/wgcna/neuron.fc/neuron.fc.manualRemove.recut.RData")
neuron=networks;rm(networks)
neuron.datExpr=neuron$datExpr
neuron.datMeta=neuron$datMeta
neuron.modules=neuron$merged$colors
neuron.genes=rownames(neuron.datExpr)  
neuron.power=8
  
shared.genes=intersect(blocks.genes,neuron.genes)
blocks.index=match(shared.genes,blocks.genes)
neuron.index=match(shared.genes,neuron.genes)

blocks.datExpr1=blocks.datExpr[blocks.index,]
blocks.modules1=blocks.modules[blocks.index]
neuron.datExpr1=neuron.datExpr[neuron.index,]
neuron.modules1=neuron.modules[neuron.index]

ME.exprBlock.moduleNeuron=moduleEigengenes(t(blocks.datExpr1),colors = neuron.modules1,softPower = blocks.power)
ME.exprNeuron.moduleBlocks=moduleEigengenes(t(neuron.datExpr1),colors = blocks.modules1,softPower = neuron.power)

# pseudo modTrait

blocks.datMeta$Diagnosis=factor(blocks.datMeta$Diagnosis,levels = c("Control","Autism"))
neuron.datMeta$Diagnosis=factor(neuron.datMeta$Diagnosis,levels = c("Control","Autism"))

modules=neuron.modules1
MEs=ME.exprBlock.moduleNeuron
datMeta=blocks.datMeta

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
modTrait.exprBlocks.moduleNeuron=modTrait


modules=blocks.modules1
MEs=ME.exprNeuron.moduleBlocks
datMeta=neuron.datMeta

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
modTrait.exprNeuron.moduleBlocks=modTrait

save(file = "working_data/wgcna/cross_datasets_me.fc.RData",
     ME.exprBlock.moduleNeuron,ME.exprNeuron.moduleBlocks,
     modTrait.exprBlocks.moduleNeuron,modTrait.exprNeuron.moduleBlocks)
