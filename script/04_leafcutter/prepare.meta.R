rm(list = ls());gc()
options(stringsAsFactors = F)

load("./working_data/wgcna/voom.forWGCNA.input.blocks.v3.RData")
block=c("Sample","Diagnosis", "Sex" , "Age" , "X260.280" , "seqPC1" , "seqPC2" , "seqPC3" , "seqPC4")

datMeta1=datMeta[,block]
datMeta1$Sample=paste(datMeta1$Sample,".merged",sep = "")
write.table(datMeta1,file = "./working_data/leafcutter/metadata/leafcutter.block.meta.txt",
            quote = F, sep = "\t", row.names = F)

rm(list = ls());gc()
options(stringsAsFactors = F)

load("./working_data/wgcna/voom.forWGCNA.input.neuron.v3.RData")
neuron=c("Sample","Diagnosis", "Sex" , "Age"  , "seqPC1" , "seqPC2" , "seqPC3")

datMeta1=datMeta[,neuron]
datMeta1$Sample=paste(datMeta1$Sample,".merged",sep = "")
write.table(datMeta1,file = "./working_data/leafcutter/metadata/leafcutter.neuron.meta.txt",
            quote = F, sep = "\t", row.names = F)


rm(list = ls());gc()
options(stringsAsFactors = F)

load("./working_data/wgcna/voom.forWGCNA.input.oligo.v3.RData")
oligo=c("Sample","Diagnosis", "Sex" , "Age"  , "seqPC1" , "seqPC2" , "seqPC3")

datMeta1=datMeta[,oligo]
datMeta1$Sample=paste(datMeta1$Sample,".merged",sep = "")
write.table(datMeta1,file = "./working_data/leafcutter/metadata/leafcutter.oligo.meta.txt",
            quote = F, sep = "\t", row.names = F)
