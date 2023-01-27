rm(list = ls());gc()
options(stringsAsFactors = F)

library(tidyverse)
library(mixOmics)

genes=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")
clusters=read.delim("working_data/leafcutter/v2/neuron_cluster_significance.txt")
cca_files=list.files("working_data/cca/",pattern = glob2rx("chr*rcca.RData"))

for (i in cca_files) {
  load(paste("working_data/cca/",i,sep = ""))
  
  loading_x=cca_shrink$loadings$X
  loading_y=cca_shrink$loadings$Y
  #plot(cca_shrink$cor) # top 42 correlations
  #cca_shrink$cor
  
  
  
  names_x = genes$genename[match(str_split_fixed(rownames(loading_x),"\\.",Inf)[,1],genes$geneid)]
  
  
  names_y=clusters$genes[match(str_split_fixed(rownames(loading_y),":",Inf)[,4] , str_split_fixed(clusters$cluster,":",Inf)[,2])]
  
  #plotVar(cca_shrink, comp = 1:2, cutoff = 0.7, var.names = c(TRUE,TRUE))
  
  loading_x1=as.data.frame(loading_x)
  loading_x1$name=names_x
  
  loading_y1=as.data.frame(loading_y)
  loading_y1$name=names_y
  
  ind=grepl("^SNOR",loading_x1$name)
  print(c(i,sum(abs(apply(loading_x[ind,],1,max)) > 0.0025)))
}
write.table(loading_x1,file = "working_data/cca/chr20_genes.txt",quote = F,sep = "\t",row.names = T)
write.table(loading_y1,file = "working_data/cca/chr20_clusters.txt",quote = F,sep = "\t",row.names = T)
