rm(list = ls());gc()
options(stringsAsFactors = F)
library(stringr)

f=read.delim("genes.results.list")
for (i in 1:nrow(f)){
  dat=read.delim(f[i,1])[,c(1,5)]
  if (exists("results")){
    dat=dat[match(results$gene_id,dat$gene_id),]
    results=cbind(results,dat[,2])
    colnames(results)[i+1]=str_split_fixed(f[i,1],"/",Inf)[1,8]
  }
  else{
    results=dat
    colnames(results)[i+1]=str_split_fixed(f[i,1],"/",Inf)[1,8]
  }
}
write.table(results,file = "genes.results.txt",quote = F,sep = "\t",row.names = F)