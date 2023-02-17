rm(list=ls());gc()
options(stringsAsFactors = F)

ego.table.all=read.delim("working_data/dge/gsea/blocks.dge.gsea.txt")
ego.table.all=read.delim("working_data/dge/gsea/neuron.dge.gsea.txt")
ego.plot=rbind(head(ego.table.all,10),tail(ego.table.all,10))

ego.plot$overlap=apply(ego.plot,1,function(x){
  dim(str_split_fixed(x[11],"/",Inf))[2]
})


ggplot(ego.plot,aes(x=reorder(Description, NES), y=NES,size=overlap,color=qvalues))+
  geom_point()+
#  geom_hline(yintercept = 0)+
  coord_flip() + 
  xlab("") + 
  ylab("Signed NES") +
  scale_color_gradient(high = "pink",low = "purple")+
  theme_light()+
  theme(text = element_text(family = "Arial",size = 20,face = "bold"),
        axis.text.y = element_text(size=16,face = "bold"),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size=18),
        legend.title = element_text(size = 14,face = "plain"),
        legend.text = element_text(size=12,face = "plain"))



# revision block interaction
ego.table.all=read.delim("working_data/dge/gsea/blocks_interaction.dge.gsea.txt")
ego.table.all=ego.table.all[ego.table.all$onto != "CC",]
ego.table.all=ego.table.all[ego.table.all$Description != "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",]
ego.plot=rbind(head(ego.table.all,5),tail(ego.table.all,5))
ego.plot$NES=abs(ego.plot$NES)

ego.plot$overlap=sapply(ego.plot$core_enrichment,function(x){
  a=str_count(x,pattern = "/")+1
  return(a)
})

ego.plot$qvalues=as.numeric(ego.plot$qvalues)
ggplot(ego.plot,aes(x=reorder(Description, NES), y=NES,size=overlap,color=qvalues))+
  geom_point()+
  #  geom_hline(yintercept = 0)+
  coord_flip() + 
  xlab("") + 
  ylab("NES") +
  scale_y_continuous(limits = c(2.2,2.4),breaks = c(2.2,2.3,2.4))+
  scale_color_gradient2(high = "#D8BFD8",low = "#EE82EE",mid = "#DDA0DD",midpoint = 0.00251)+
  theme_light()+
  theme(text = element_text(family = "Arial",size = 20,face = "bold"),
        axis.text.y = element_text(size=14,face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size=18),
        legend.title = element_text(size = 14,face = "plain"),
        legend.text = element_text(size=12,face = "plain"))
