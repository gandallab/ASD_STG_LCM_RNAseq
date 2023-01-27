rm(list = ls());gc()
options(stringsAsFactors = F)

library(tidyverse)

# load cell type enrichment
load("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/ewce.neuron.level2.RData")
modTrait=read.delim("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/modTrait.neuron.fc2pass.txt")

joint_ct=adult_all[adult_all$moduleColor %in% modTrait$Module[modTrait$fdr < 0.05],]
joint_ct=joint_ct[grepl("Ex|In",joint_ct$CellType),]
joint_ct$neurontype=str_sub(joint_ct$CellType,1,2)
joint_ct$neurontype[joint_ct$neurontype == "Ex"]="Excitatory neurons"
joint_ct$neurontype[joint_ct$neurontype == "In"]="Inhibitory neurons"
joint_ct$sd_from_mean[joint_ct$sd_from_mean < 0]=0.01
joint_ct$symbol=""
joint_ct$symbol[joint_ct$fdr < 0.05]="*"
ggplot(joint_ct)+
  geom_col(aes(CellType,sd_from_mean,fill=moduleColor),width = 1,position = "dodge",color="black")+
  geom_text(aes(CellType,sd_from_mean + 0.2,label=symbol),color="red",size = 8)+
  scale_fill_manual(values = c("green"="green","red"="red","greenyellow"="greenyellow","lightcyan"="lightcyan","grey60"="grey60"),
                    labels=c("Neu-M5","Neu-M6","Neu-M11","Neu-M16","Neu-M17"))+
  facet_wrap(~neurontype,scales = "free_x")+
  xlab("")+ylab("Z-score")+
  guides(fill=guide_legend(title="Neuron Modules"))+
  theme_light()+
  theme(text=element_text(family="Arial",face = "bold"),
        panel.background = element_rect(colour = "black"),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle = 90,hjust = 0.5,vjust = 0.5))
