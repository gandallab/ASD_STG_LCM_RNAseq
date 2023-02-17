options(stringsAsFactors = F)
library(tidyverse)

rm(list = ls());gc()

load("working_data/snoRNA/individual_correlation/snorna_direct_cor_splicing_individual_results.RData")

#dat = allp[allp$newfdr < 0.05,]
dat = allp

dat$intron=paste("chr",dat$event,sep = "")

dat$cluster=sapply(dat$event,function(x){
  a=str_split_fixed(x,":",Inf)[1,]
  paste("chr",a[1],":",a[length(a)],sep = "")
})

sg=read.delim("working_data/leafcutter/v2/neuron_cluster_significance.txt")

table(dat$cluster %in% sg$cluster)

dat$clsuter_gene=sg$genes[match(dat$cluster,sg$cluster)]
dat$cluster_fdr=sg$p.adjust[match(dat$cluster,sg$cluster)]

ef=read.delim("working_data/leafcutter/v2/neuron_effect_sizes.txt")
table(dat$intron %in% ef$intron)
table(ef$intron %in% dat$intron)

sigdat=dat[dat$fdr < 0.05,]
sigdat=sigdat[!duplicated(sigdat$clsuter_gene),]
sigdat=sigdat[!is.na(sigdat$clsuter_gene),]

nosigdat=dat[dat$fdr >= 0.05,]
nosigdat=nosigdat[!duplicated(nosigdat$clsuter_gene),]
nosigdat=nosigdat[!is.na(nosigdat$clsuter_gene),]

length(intersect(sigdat$clsuter_gene,nosigdat$clsuter_gene))

dat1 = dat[!is.na(dat$clsuter_gene),]

library(gProfileR)
go = gprofiler(query=sigdat$clsuter_gene, 
               max_set_size = 500,
               min_isect_size = 2,
               correction_method = "fdr",
               exclude_iea = T,
               hier_filtering = "moderate", 
               custom_bg = unique(dat1$clsuter_gene), 
               src_filter = c("GO","KEGG"),
               ordered_query = F)
