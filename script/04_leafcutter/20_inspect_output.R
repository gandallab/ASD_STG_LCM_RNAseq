rm(list = ls());gc()
options(stringsAsFactors = F)

load("manual_leafviz/leafviz_data/v2/blocks.leafviz.RData")
blocks_intron=introns
blocks_cluster=clusters
load("manual_leafviz/leafviz_data/v2/neuron.leafviz.RData")
neuron_intron=introns
neuron_cluster=clusters


intersect(blocks_cluster$clusterID,neuron_cluster$clusterID)
blocks_intron$id=paste(blocks_intron$clusterID,blocks_intron$chr,blocks_intron$start,blocks_intron$end,sep = ":")
neuron_intron$id=paste(neuron_intron$clusterID,neuron_intron$chr,neuron_intron$start,neuron_intron$end,sep = ":")
dat=merge(blocks_intron,neuron_intron,by="id")
ggplot(dat,aes(deltapsi.x,deltapsi.y))+
  geom_point()+
  xlim(-0.2,0.2)+
  ylim(-0.2,0.2)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)


rm(list = ls());gc()
options(stringsAsFactors = F)

blocks=read.delim("working_data/leafcutter/v2/blocks_effect_sizes.txt")
tmp=blocks[1:10,]
neuron=read.delim("working_data/leafcutter/v2/neuron_effect_sizes.txt")
table(neuron$intron %in% blocks$intron)
table(blocks$intron %in% neuron$intron)
dat=merge(blocks,neuron,by="intron")
ggplot(dat,aes(deltapsi.x,deltapsi.y))+
  geom_point()+
  xlim(-0.2,0.2)+
  ylim(-0.75,0.75)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)

blocks=read.delim("working_data/leafcutter/v2/blocks_cluster_significance.txt")
neuron=read.delim("working_data/leafcutter/v2/neuron_cluster_significance.txt")
blocks=blocks[blocks$status == "Success",]
neuron=neuron[neuron$status == "Success",]
table(blocks$p.adjust < 0.05)
table(neuron$p.adjust < 0.05)
