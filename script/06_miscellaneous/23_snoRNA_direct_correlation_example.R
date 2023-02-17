options(stringsAsFactors = F)
library(tidyverse)

rm(list = ls());gc()
dat=read.delim("working_data/snoRNA/individual_correlation/got1.psi.vector.txt")

rownames(dat)=dat$ID
dat1=dat[,-c(1:4)]

datMeta=read.delim("working_data/snoRNA/individual_correlation/neuron.group.txt",header = F)
sum(colnames(dat1) == datMeta$V1) 

X = model.matrix(~ V2 + V3 + V4 + V5 + V6 + V7 + V8,datMeta)
Y = dat1
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
#datExpr = Y - t(as.matrix(X[,5]) %*% t(as.matrix(beta[5,])))
datSp = Y - t(X[,c(5:ncol(X))] %*% beta[c(5:nrow(beta)),])

load("working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/neuron.fc2pass.recut.RData")
datExpr=networks$datExpr
rownames(datExpr)=sapply(rownames(datExpr),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
rm(networks)

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


dat$clsuter_gene=sg$genes[match(dat$cluster,sg$cluster)]
dat$cluster_fdr=sg$p.adjust[match(dat$cluster,sg$cluster)]

ef=read.delim("working_data/leafcutter/v2/neuron_effect_sizes.txt")
table(dat$intron %in% ef$intron)

dat1=dat[dat$cluster_fdr < 0.1 & !is.na(dat$cluster_fdr) & dat$intron %in% ef$intron,]

targetsnorna=dat1[dat1$clsuter_gene == "GOT1" & !is.na(dat1$clsuter_gene) & dat1$rho >0,]

datExpr1=datExpr[targetsnorna$snorna_id,]
datExpr2 = as.data.frame(t(datExpr1)) %>%
  rownames_to_column("sample") %>%
  gather("snorna_id","RNA_expression",2:12)
datExpr2$sample = make.names(datExpr2$sample)

datSp1=t(datSp)

plotdata=merge(datExpr2,datSp1,by.x="sample",by.y="row.names")
colnames(plotdata)[4]="GOT1_clu_16816"
plotdata$snoRNAs=annot$genename[match(plotdata$snorna_id,annot$geneid)]

ggplot(plotdata,aes(GOT1_clu_16816,RNA_expression,color=snoRNAs,fill=snoRNAs))+
  geom_point()+
  geom_smooth(method = "lm",alpha=0.1)+
  xlab("GOT1 Intron chr10:101163631-101165513 PSI")+
  ylab("snoRNA expression")+
  theme_bw()+
  theme(text = element_text(family = "Arial",face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=16))

plotdata1=datMeta
plotdata1$V2=factor(plotdata1$V2,levels = c("Control","Autism"))
datSp2=as.data.frame(datSp1)
plotdata1$got1=datSp2$`10:101163631:101165513:clu_16816_NA`[match(datMeta$V1,row.names(datSp2))]
ggplot(plotdata1,aes(V2,got1,fill=V2))+
  geom_violin(width=0.6)+
  geom_jitter(width = 0.1,size=2)+
  guides(fill=FALSE)+
  xlab("")+
  ylab("GOT1 Intron\nchr10:101163631-101165513 PSI")+
  theme_bw()+
  theme(text = element_text(family = "Arial",face = "bold"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size=16))
t.test(plotdata1$got1[plotdata1$V2=="Control"],plotdata1$got1[plotdata1$V2=="Autism"])
