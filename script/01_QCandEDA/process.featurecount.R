options(stringsAsFactors = F)

rm(list = ls());gc()

counts=list.files("C:/Users/Pan Zhang/Downloads/fc/countsfiles_block_unique/","*counts.txt")
for (count in counts){
  dat=read.delim(paste("C:/Users/Pan Zhang/Downloads/fc/countsfiles_block_unique/",count,sep = ""),skip = 1)
  dat=dat[,c(1,7,8)]
  colnames(dat)[3]=gsub(".counts.txt","",count)
  if (exists("output")){
    if (output$Geneid != dat$Geneid){
      print(count)
      next()
    }
    output=cbind(output,dat[,3,drop=F])
  }
  else{
    output=dat
  }
}

rownames(output)=output$Geneid;output$Geneid=NULL;output$gene_name=NULL
output1=as.matrix(output)



library(SummarizedExperiment)
load("working_data/summarizedExperiment/newrsem/se_neuron.RData")
datMeta=as.data.frame(colData(se_neuron))
table(datMeta$Sample %in% colnames(output1))
output2=output1[,datMeta$Sample]
output3=output2[rowSums(output2>0)>0.8*ncol(output2),]
se_neuron=SummarizedExperiment(assays = list(counts=output3),colData=datMeta)
save(file = "working_data/summarizedExperiment/fc_multimap//se_neuron.RData",se_neuron)

load("working_data/summarizedExperiment/newrsem/se_oligo.RData")
datMeta=as.data.frame(colData(se_oligo))
table(datMeta$Sample %in% colnames(output1))
output2=output1[,datMeta$Sample]
output3=output2[rowSums(output2>0)>0.8*ncol(output2),]
se_oligo=SummarizedExperiment(assays = list(counts=output3),colData=datMeta)
save(file = "working_data/summarizedExperiment/fc_multimap//se_oligo.RData",se_oligo)

load("working_data/summarizedExperiment/newrsem/se_blocks.RData")
datMeta=as.data.frame(colData(se_blocks))
table(datMeta$Sample %in% colnames(output1))
output2=output1[,datMeta$Sample]
output3=output2[rowSums(output2>0)>0.8*ncol(output2),]
se_blocks=SummarizedExperiment(assays = list(counts=output3),colData=datMeta)
save(file = "working_data/summarizedExperiment/fc//se_blocks.RData",se_blocks)
