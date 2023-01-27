rm(list = ls());gc()
options(stringsAsFactors = F)

library(tidyverse)
library(mixOmics)

nComp = 30
genes=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt")
clusters=read.delim("working_data/leafcutter/v2/neuron_cluster_significance.txt")
cca_files=list.files("working_data/cca/",pattern = glob2rx("chr*rcca.RData"))
final_output=data.frame()
for (i in cca_files) {
  # chr=gsub("_rcca.RData","",i)
  load(paste("working_data/cca/",i,sep = ""))
  
  expr_x=cca_shrink$X
  expr_y=cca_shrink$Y
  #plot(cca_shrink$cor) # top 42 correlations
  print(cca_shrink$cor[nComp])
  
  
  
  names_x = genes$genename[match(str_split_fixed(colnames(expr_x),"\\.",Inf)[,1],genes$geneid)]
  
  
  names_y=clusters$genes[match(str_split_fixed(colnames(expr_y),":",Inf)[,4] , str_split_fixed(clusters$cluster,":",Inf)[,2])]
  
  ind=grepl("^SNOR",names_x)
  expr_x1=expr_x[,ind]
  colnames(expr_x1)=names_x[ind]
  variate_x=cca_shrink$variates$X
  #variate_y=cca_shrink$variates$Y
  #cor(variate_x[,1],variate_y[,1])
  
  output=matrix(NA,ncol = min(ncol(variate_x),nComp),nrow = ncol(expr_x1))
  for (j in 1:ncol(expr_x1)){
    for (k in 1:min(ncol(variate_x),nComp)){
      current_cor=cor(expr_x1[,j],variate_x[,k])
      output[j,k]=current_cor
    }
  }
  rownames(output)=colnames(expr_x1)
  
  output1=matrix(NA,ncol = min(ncol(variate_x),nComp),nrow = ncol(expr_y))
  for (j in 1:ncol(expr_y)){
    for (k in 1:min(ncol(variate_x),nComp)){
      current_cor=cor(expr_y[,j],variate_x[,k])
      output1[j,k]=current_cor
    }
  }
  rownames(output1)=colnames(expr_y)
  
  output2=as.data.frame(cor(t(output),t(output1)))
  output3=output2 %>% 
    tibble::rownames_to_column("gene") %>%
    tidyr::gather("cluster","correlation",2:(ncol(output2)+1))
  output3$cluster_gene=clusters$genes[match(str_split_fixed(output3$cluster,":",Inf)[,4] , str_split_fixed(clusters$cluster,":",Inf)[,2])]
  
  final_output=rbind(final_output,output3)
}
save(final_output,file = "working_data/cca/variable_cor.RData")
final_output1=final_output %>% filter(abs(correlation) > 0.7)
#final_output1=final_output %>% filter(grepl("^14:",cluster))
write.table(final_output1,file = "working_data/cca/top07correlation.txt",quote = F,sep = "\t",row.names = F)
