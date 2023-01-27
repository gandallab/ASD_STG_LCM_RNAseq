options(stringsAsFactors = F)
rm(list = ls());gc()

library(tidyverse)

# blocks
input=read.delim("working_data/rmats/block.out/filter.combined.events.block.txt")
datMeta=read.delim("working_data/summarizedExperiment/fc/datMeta.blocks.txt")
ctr_samples=read.delim("working_data/rmats/block.out/block.ctr.samples.txt")
asd_samples=read.delim("working_data/rmats/block.out/block.asd.samples.txt")
sample_seq=c(ctr_samples$V1,asd_samples$V1)
datMeta1=datMeta[match(sample_seq,datMeta$Sample),]
output=data.frame()
for (i in 1:nrow(input)){
  inc_ctr=as.numeric(str_split_fixed(input[i,6],",",Inf)[1,])
  inc_asd=as.numeric(str_split_fixed(input[i,7],",",Inf)[1,])
  inc_all=c(inc_ctr,inc_asd)
  current=cbind(inc_all,datMeta1)
  current=current[!is.na(current$inc_all),]
  if (length(unique(current$Sex))==1){
    fml=as.formula("inc_all ~ Diagnosis + Age + RIN + X260.280 + seqPC1 + seqPC2 + seqPC3")
    unisex="yes"
  }else{
    fml=as.formula("inc_all ~ Diagnosis + Sex + Age + RIN + X260.280 + seqPC1 + seqPC2 + seqPC3")
    unisex="no"
  }
  current$Diagnosis=factor(current$Diagnosis,levels = c("Control","Autism"))
  rm(inc_all)
  mylm=summary(lm(fml,data = current))$coefficients
  output=rbind(output,mylm["DiagnosisAutism",])
}
colnames(output)=colnames(mylm)
result_block=cbind(input,output)
#result_block$fdr=p.adjust(result_block$`Pr(>|t|)`,method = "fdr")
#write.table(result_block,file = "working_data/rmats/block.results.txt",quote = F,sep = "\t",row.names = F)
result_block_fdr=data.frame()
for (i in unique(result_block$event)){
  current=result_block[result_block$event == i,]
  current$fdr=p.adjust(current$Pr...t..,method = "fdr")
  result_block_fdr=rbind(result_block_fdr,current)
}
write.table(result_block_fdr,file = "working_data/rmats/block.results.txt",quote = F,sep = "\t",row.names = F)

# neuron
input=read.delim("working_data/rmats/neuron.out/filter.combined.events.neuron.txt")
datMeta=read.delim("working_data/summarizedExperiment/fc/datMeta.neuron.txt")
ctr_samples=read.delim("working_data/rmats/neuron.out/neuron.ctr.samples.txt")
asd_samples=read.delim("working_data/rmats/neuron.out/neuron.asd.samples.txt")
sample_seq=c(ctr_samples$V1,asd_samples$V1)
datMeta1=datMeta[match(sample_seq,datMeta$Sample),]
output=data.frame()
for (i in 1:nrow(input)){
  inc_ctr=as.numeric(str_split_fixed(input[i,6],",",Inf)[1,])
  inc_asd=as.numeric(str_split_fixed(input[i,7],",",Inf)[1,])
  inc_all=c(inc_ctr,inc_asd)
  current=cbind(inc_all,datMeta1)
  current=current[!current$Sample %in% c("AO_B02","6221_N"),]
  current=current[!is.na(current$inc_all),]
  if (length(unique(current$Sex))==1){
    fml=as.formula("inc_all ~ Diagnosis + Age + Type_RNAseqRunNumber + seqPC1 + seqPC2 + seqPC3")
    unisex="yes"
  }else{
    fml=as.formula("inc_all ~ Diagnosis + Sex + Age + Type_RNAseqRunNumber + seqPC1 + seqPC2 + seqPC3")
    unisex="no"
  }
  current$Diagnosis=factor(current$Diagnosis,levels = c("Control","Autism"))
  rm(inc_all)
  mylm=summary(lm(fml,data = current))$coefficients
  output=rbind(output,mylm["DiagnosisAutism",])
}

colnames(output)=colnames(mylm)
result_neuron=cbind(input,output)
#result_block$fdr=p.adjust(result_block$`Pr(>|t|)`,method = "fdr")
#write.table(result_neuron,file = "working_data/rmats/neuron.results.txt",quote = F,sep = "\t",row.names = F)

result_neuron_fdr=data.frame()
for (i in unique(result_neuron$event)){
  current=result_neuron[result_neuron$event == i,]
  current$fdr=p.adjust(current$`Pr(>|t|)`,method = "fdr")
  result_neuron_fdr=rbind(result_neuron_fdr,current)
}
write.table(result_neuron_fdr,file = "working_data/rmats/neuron.results.txt",quote = F,sep = "\t",row.names = F)

