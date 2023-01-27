rm(list = ls());gc()
options(stringsAsFactors = F)

library(tidyverse)
library(mixOmics)
library(igraph)
library(RCy3)

genes=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")
clusters=read.delim("working_data/leafcutter/v2/neuron_cluster_significance.txt")
cca_files=list.files("working_data/cca/",pattern = glob2rx("chr*rcca.RData"))

for (i in cca_files) {
  load(paste("working_data/cca/",i,sep = ""))
  
  names_x = genes$genename[match(str_split_fixed(cca_shrink$names$colnames$X,"\\.",Inf)[,1],genes$geneid)]
  types_x = genes$genetype[match(str_split_fixed(cca_shrink$names$colnames$X,"\\.",Inf)[,1],genes$geneid)]
  names_x[duplicated(names_x)]=paste(names_x[duplicated(names_x)],"v2",sep="_")
  
  names_y=clusters$genes[match(str_split_fixed(cca_shrink$names$colnames$Y,":",Inf)[,4] , str_split_fixed(clusters$cluster,":",Inf)[,2])]
  names_y=paste(names_y,sample(1:20000,length(names_y),replace = F),sep = "_")
  
  nw=network(cca_shrink, row.names = names_x , col.names = names_y ,comp = 1:35, cutoff = 0.8)
  save(nw,file = paste("working_data/cca/network/",gsub("_rcca.RData","",i),"_network.RData",sep = ""))
  #table(dimnames(nw$M)[[1]] == cca_shrink$names$colnames$X)
  rownames(nw$M)=names_x
  colnames(nw$M)=names_y
  sub_mat=nw$M[types_x == "snoRNA",]
  print(c(i,sum(rowSums(sub_mat != 0) > 0)))
  vc=rownames(sub_mat)[rowSums(sub_mat != 0) > 0]
  eg=E(nw$gR)[adj(vc)]
  sg=subgraph.edges(nw$gR,eg)
  save(sg,file = paste("working_data/cca/network/",gsub("_rcca.RData","",i),"_subnetwork.RData",sep = ""))
  #createNetworkFromIgraph(sg,"tmp")
  #plot.igraph(sg,layout=layout_as_bipartite(sg,types = grepl("^ENSG",V(sg)$name)))
  
}
write.table(loading_x1,file = "working_data/cca/chr20_genes.txt",quote = F,sep = "\t",row.names = T)
write.table(loading_y1,file = "working_data/cca/chr20_clusters.txt",quote = F,sep = "\t",row.names = T)
