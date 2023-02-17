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

sg1=sg1[complete.cases(sg1),]

# ef1=ef %>%
#   group_by(cluster) %>%
#   filter(abs(deltapsi) == max(abs(deltapsi))) %>%
#   top_n(1)

library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
library(caret)

sg2=sg1[order(sg1$p,decreasing = T),]
sg2$rank=seq(1,nrow(sg2))
geneList=sg2$rank 
names(geneList)=sg2$genes
geneList=sort(geneList,decreasing = T)
ego.table.all=data.frame()
for (onto in c("BP","MF","CC")){
  ego <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               keyType = "SYMBOL",
               ont          = onto,
               nPerm        = 500,
               minGSSize    = 50,
               maxGSSize    = 300,
               pvalueCutoff = 1,
               verbose      = TRUE)
  
  ego.table=ego@result
  ego.table=ego.table[complete.cases(ego.table),]
  ego.table=ego.table[ego.table$qvalues < 0.2,]
  if (nrow(ego.table) == 0){next}
  # hsGO <- godata('org.Hs.eg.db', ont=onto)
  # go_cor=mgoSim(ego.table$ID,ego.table$ID,semData = hsGO,measure = "Wang",combine = NULL)
  # to_remove=findCorrelation(go_cor,cutoff = 0.1,names = T)
  # ego.table=ego.table[!ego.table$ID %in% to_remove,]
  ego.table$onto=onto
  ego.table.all=rbind(ego.table.all,ego.table)
}


ego.table.all=ego.table.all[order(ego.table.all$NES,decreasing = T),]

write.table(ego.table.all,file = "working_data/leafcutter/v2/gsea/blocks.leafcutter.gsea.txt",quote = F,sep = "\t",row.names = F)
ego.plot=ego.table.all


ggplot(ego.plot, aes(x=reorder(Description, NES), y=NES,fill=(NES < 0))) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  xlab("") + 
  ylab("NES") +
  scale_fill_manual(values = c("steelblue","royalblue"))+
  theme(text = element_text(family = "Arial",size = 15,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20))
