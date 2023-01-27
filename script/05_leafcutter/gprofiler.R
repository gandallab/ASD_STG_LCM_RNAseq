rm(list = ls());gc()
options(stringsAsFactors = F)

library(gProfileR)
library(tidyverse)

dat=read.delim("working_data/leafcutter/output/block_ds_cluster_significance.txt")
# table(complete.cases(dat))
# table(dat$status)
dat1=dat[dat$status == "Success",]
# table(complete.cases(dat1))
dat1=dat[complete.cases(dat),]
table(dat1$status)

dat2=dat1[order(dat1$p),]
table(duplicated(dat2$genes))
dat3=dat2[!duplicated(dat2$genes),]

dat4=dat3[dat3$p.adjust < 0.05,]

go = gprofiler(query=dat4$genes, 
                 max_set_size = 500,
                 min_isect_size = 2,
                 correction_method = "fdr",
                 exclude_iea = T,
                 hier_filtering = "moderate", 
                 custom_bg = dat3$genes, 
                 src_filter = c("GO","KEGG"),
                 ordered_query = T)

go=go[order(go$p.value),]


ggplot(go[1:1,], aes(x=reorder(term.name, -p.value), y=-log10(p.value))) + 
  geom_bar(stat="identity") + 
  coord_flip() + 
  xlab("") + 
  ylab("-log10(FDR)") +
  geom_hline(yintercept=c(-log10(0.05)), lty=2, color="red")+
  theme(text = element_text(family = "Arial",size = 15,face = "bold"),
        axis.text.y = element_text(size=12,face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size=20))+
  guides(fill=F)
