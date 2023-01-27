rm(list = ls());gc()
options(stringsAsFactors = F)

library(tidyverse)

ef=read.delim("working_data/leafcutter/v2/neuron_effect_sizes.txt")
sg=read.delim("working_data/leafcutter/v2/neuron_cluster_significance.txt")
ef=ef[abs(ef$deltapsi)>0.1,]
sg=sg[sg$status == "Success",]
ef$cluster=sapply(ef$intron,function(x){
  a=str_split_fixed(x,":",Inf)[1,]
  paste(a[1],a[length(a)],sep = ":")
})

sg=sg[sg$cluster %in% ef$cluster,]

sg1=sg %>%
  separate_rows(genes,sep = ",") %>%
  group_by(genes) %>%
  filter(p == min(p))
sg1=sg1[!duplicated(sg1$genes),]
sg2=sg1[order(sg1$p.adjust,decreasing = F),]

library(gProfileR)

go = gprofiler(query=sg2$genes[sg2$p.adjust < 0.01], 
               max_set_size = 500,
               min_isect_size = 2,
               correction_method = "fdr",
               exclude_iea = T,
               hier_filtering = "moderate", 
               custom_bg = sg2$genes, 
               src_filter = c("GO","KEGG"),
               ordered_query = F)

go=go[order(go$p.value),]
sg1=sg1[complete.cases(sg1),]
