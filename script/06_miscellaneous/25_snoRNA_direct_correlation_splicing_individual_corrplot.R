options(stringsAsFactors = F)
library(tidyverse)
library(ggcorrplot)
rm(list = ls());gc()


load("working_data/snoRNA/individual_correlation/snorna_direct_cor_splicing_individual_results.RData")


dat = allp[allp$newfdr < 0.05,]
annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt",header = T,sep = "\t")
dat$snorna_name=annot$genename[match(dat$snorna_id,annot$geneid)]

dat$cluster=sapply(dat$event,function(x){
  a=str_split_fixed(x,":",Inf)[1,]
  paste("chr",a[1],":",a[length(a)],sep = "")
})

dat$intron=paste("chr",dat$event,sep = "")


sg=read.delim("working_data/leafcutter/v2/neuron_cluster_significance.txt")

#table(dat$cluster %in% sg$cluster)

sg=sg[sg$cluster %in% dat$cluster,]

table(sg$status)

dat$clsuter_gene=sg$genes[match(dat$cluster,sg$cluster)]
dat$cluster_fdr=sg$p.adjust[match(dat$cluster,sg$cluster)]

ef=read.delim("working_data/leafcutter/v2/neuron_effect_sizes.txt")
table(dat$intron %in% ef$intron)


dat1=dat[dat$cluster_fdr < 0.1 & !is.na(dat$cluster_fdr) & dat$intron %in% ef$intron,]
length(unique(dat1$clsuter_gene))

dat2=dat1[,c("event","cluster","intron","clsuter_gene","cluster_fdr","snorna_id","snorna_name","rho","newfdr")]
write.table(dat2,file = "working_data/snoRNA/individual_correlation/snoRNA_intron_significant_correlations.txt",
            quote = F,sep = "\t",row.names = F)
dat3=dat1[,c("snorna_name","clsuter_gene","rho")]
dat3=dat3[complete.cases(dat3),]
dat4=dat3 %>% 
  group_by(snorna_name,clsuter_gene) %>%
  slice_max(abs(rho),1)

dat5 = dat4 %>%
  ungroup() %>%
  group_by(clsuter_gene) %>%
  summarise(n=n())

dat6=dat4[dat4$clsuter_gene %in% dat5$clsuter_gene[dat5$n > 3],] %>%
  spread(clsuter_gene,rho)
dat7=as.data.frame(dat6);rownames(dat7)=dat7$snorna_name;dat7$snorna_name=NULL;dat7=as.matrix(dat7)
dat8=dat7;dat8[is.na(dat8)]=0
hc_snorna=hclust(dist(dat8),method = "ward.D2")
hc_gene=hclust(dist(t(dat8)),method = "ward.D2")
dat9=dat7[hc_snorna$order,hc_gene$order]

ggcorrplot(dat9,
           method = "circle",
           tl.srt = 90,
           tl.cex = 12,
           ggtheme = ggplot2::theme_light(base_family = "Arial") + theme(text = element_text(face = "bold"))
           )

