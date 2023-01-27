options(stringsAsFactors = F)
library(tidyverse)

rm(list = ls());gc()

rhos=list.files("working_data/snoRNA/individual_correlation/",glob2rx("chr*rho*txt"),full.names = TRUE)
allrho=data.frame()

for (rho in rhos){
  current = read.delim(rho,header = T,row.names = 1)
  current1 = current %>%
    rownames_to_column("event") %>%
    gather("snorna_id","rho",2:41)
  allrho=rbind(allrho,current1)
}

ps=list.files("working_data/snoRNA/individual_correlation/",glob2rx("chr*cor_p*txt"),full.names = TRUE)
allp=data.frame()

for (p in ps){
  current = read.delim(p,header = T,row.names = 1)
  current1 = current %>%
    rownames_to_column("event") %>%
    gather("snorna_id","p",2:41)
  allp=rbind(allp,current1)
}

allp$fdr=1
for (snorna in unique(allp$snorna_id)){
  allp$fdr[allp$snorna_id == snorna]=p.adjust(allp$p[allp$snorna_id == snorna],method = "fdr")
}
allp$newfdr = p.adjust(allp$p,method = "fdr")

table(allrho$event == allp$event)
table(allrho$snorna_id == allp$snorna_id)

allp$rho=allrho$rho
save(allp,file = "working_data/snoRNA/individual_correlation/snorna_direct_cor_splicing_individual_results.RData")


dat = allp[allp$newfdr < 0.05,]
annot=read.delim("D:/references/gencode.v29lift37/gencode.v29lift37.genename.genetype.txt",header = T,sep = "\t")
dat$snorna_name=annot$genename[match(dat$snorna_id,annot$geneid)]

dat$cluster=sapply(dat$event,function(x){
  a=str_split_fixed(x,":",Inf)[1,]
  paste("chr",a[1],":",a[length(a)],sep = "")
})

dat$intron=paste("chr",dat$event,sep = "")


sg=read.delim("working_data/leafcutter/v2/neuron_cluster_significance.txt")

table(dat$cluster %in% sg$cluster)

sg=sg[sg$cluster %in% dat$cluster,]

table(sg$status)

dat$clsuter_gene=sg$genes[match(dat$cluster,sg$cluster)]
dat$cluster_fdr=sg$p.adjust[match(dat$cluster,sg$cluster)]

ef=read.delim("working_data/leafcutter/v2/neuron_effect_sizes.txt")
table(dat$intron %in% ef$intron)

dat1=dat[dat$cluster_fdr < 0.1 & !is.na(dat$cluster_fdr) & dat$intron %in% ef$intron,]
length(unique(dat1$clsuter_gene))

gmt=read.gmt("D:/datasets/genesets/houseLists03242020.gmt")
gmt=gmt[gmt$ont %in% c("Syndromic ASD","Presynaptic","Postsynaptic","pLI>0.99"),]

dat2=merge(dat1,gmt,by.x="clsuter_gene",by.y="gene")
length(unique(dat2$clsuter_gene))

explorer=data.frame()
for (gene in unique(dat2$clsuter_gene)){
  current=dat2[dat2$clsuter_gene == gene,]
  explorer=rbind(explorer,data.frame(gene=gene,nSNORNA=length(unique(current$snorna_id)),nCategory=length(unique(current$ont))))
}

dat3=dat1[,c("snorna_name","clsuter_gene","rho")]
dat3=dat3[complete.cases(dat3),]
dat3$direction = "pos"
dat3$direction[dat3$rho < 0]="neg"
dat4=data.frame(nodename=unique(c(dat3$snorna_name,dat3$clsuter_gene)),snorna="N",showlabel="",direction="pos")
dat4$snorna[dat4$nodename %in% dat3$snorna_name]="Y"
dat4$showlabel[dat4$nodename %in% explorer$gene[explorer$nSNORNA >1 & explorer$nCategory >1]]=dat4$nodename[dat4$nodename %in% explorer$gene[explorer$nSNORNA >1 & explorer$nCategory >1]]
dat4$showlabel[dat4$nodename %in% dat3$snorna_name]=dat4$nodename[dat4$nodename %in% dat3$snorna_name]
dat4$direction[dat4$nodename %in% dat3$clsuter_gene[dat3$rho < 0]]="neg"
write.table(dat3,"working_data/snoRNA/individual_correlation/cytoscape.snornaCor.edgeList.txt",quote = F,sep = "\t",row.names = F)
write.table(dat4,"working_data/snoRNA/individual_correlation/cytoscape.snornaCor.nodeProperty.txt",quote = F,sep = "\t",row.names = F)




# new workflow
options(stringsAsFactors = F)
library(tidyverse)

rm(list = ls());gc()
load("working_data/snoRNA/individual_correlation/snorna_direct_cor_splicing_individual_results.RData")

dat = allp



dat$intron=paste("chr",dat$event,sep = "")

ef=read.delim("working_data/leafcutter/v2/neuron_effect_sizes.txt")
table(dat$intron %in% ef$intron)

dat = dat[dat$intron %in% ef$intron,]

length(unique(dat$intron))

dat$cluster=sapply(dat$event,function(x){
  a=str_split_fixed(x,":",Inf)[1,]
  paste("chr",a[1],":",a[length(a)],sep = "")
})

sg=read.delim("working_data/leafcutter/v2/neuron_cluster_significance.txt")

table(dat$cluster %in% sg$cluster)

dat$clsuter_gene=sg$genes[match(dat$cluster,sg$cluster)]
dat$cluster_fdr=sg$p.adjust[match(dat$cluster,sg$cluster)]

save(dat,file = "working_data/snoRNA/individual_correlation/snorna_direct_cor_splicing_individual_results_matchToIntrons.RData")

sigdat=dat[dat$newfdr < 0.05,]
sigdat=sigdat[!duplicated(sigdat$clsuter_gene),]
sigdat=sigdat[!is.na(sigdat$clsuter_gene),]

nosigdat=dat[dat$newfdr >= 0.05,]
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
dat1=dat[dat$cluster_fdr < 0.1 & !is.na(dat$cluster_fdr) & dat$intron %in% ef$intron,]
length(unique(dat1$clsuter_gene))