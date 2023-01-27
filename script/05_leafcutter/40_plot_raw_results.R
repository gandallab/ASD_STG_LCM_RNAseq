rm(list = ls());gc()
options(stringsAsFactors = F)
library(tidyverse)
library(ggrepel)

# block
dat=read.delim("working_data/final/tables/leafcutter.blocks_cluster_significance.txt")
#dat=read.delim("working_data/final/tables/leafcutter.neuron_cluster_significance.txt")
dat1=dat[dat$status == "Success",]

dat2=dat1[order(dat1$p,decreasing = F),]
dat2$expected = -log10(ppoints(nrow(dat1)))
dat2$observed = -log10(dat2$p)

dat2$genes[dat2$p.adjust > 0.01]=""

# ggplot(dat1,aes(sample = -log10(p.adjust)))+
#   geom_qq()+
# #  geom_qq_line()+
#   geom_abline(slope = 1)+
#   xlim(0,5)+
#   ylim(0,8)
# 
# df <- data.frame(
#   observed=-log10(sort(dat1$p)),
#   expected = -log10(ppoints(nrow(dat1)))
# )
# observed=-log10(sort(dat1$p))
# expected = -log10(ppoints(nrow(dat1)))

ggplot(dat2) +
  geom_point(aes(expected, observed, col = (p.adjust < 0.05)), shape = 1, size = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"))+
  geom_abline(intercept = 0, slope = 1, alpha = 0.5)+
  geom_text_repel(data = dat2[c(1:17),],aes(expected, observed, label=genes),fontface="bold",size=4)+
  theme_classic()+
  theme(text=element_text(family = "Arial",face = "bold"),
        panel.background = element_rect(colour = "black"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size=20))+
  guides(color = F)



# plot block and neuron together
rm(list = ls());gc()
options(stringsAsFactors = F)
library(tidyverse)
library(ggrepel)

dat=read.delim("working_data/final/tables/leafcutter.neuron_cluster_significance.txt")
dat1=dat[dat$status == "Success",]

dat2=dat1[order(dat1$p,decreasing = F),]
dat2$expected = -log10(ppoints(nrow(dat1)))
dat2$observed = -log10(dat2$p)

dat2$genes[dat2$p.adjust > 0.01]=""

p=ggplot(dat2) +
  geom_point(aes(expected, observed, col = (p.adjust < 0.05)), shape = 1, size = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray"))+
  geom_abline(intercept = 0, slope = 1, alpha = 0.5)+
  geom_text_repel(data = dat2[c(1:17),],aes(expected, observed, label=genes),fontface="bold",size=4)+
  theme_classic()+
  theme(text=element_text(family = "Arial",face = "bold"),
        panel.background = element_rect(colour = "black"),
        axis.text = element_text(size = 18),
        axis.title = element_text(size=20))+
  guides(color = F)

dat=read.delim("working_data/final/tables/leafcutter.blocks_cluster_significance.txt")
dat1=dat[dat$status == "Success",]

dat2=dat1[order(dat1$p,decreasing = F),]
dat2$expected = -log10(ppoints(nrow(dat1)))
dat2$observed = -log10(dat2$p)

dat2$genes[dat2$p.adjust > 0.01]=""

p+geom_point(data=dat2,aes(expected, observed, col = (p.adjust < 0.05)), shape = 1, size = 1,alpha=0.08)
neuron=dat2[,c("p.adjust","genes","expected","observed")]
neuron$dataset="neuron"