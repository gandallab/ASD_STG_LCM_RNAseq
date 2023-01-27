rm(list = ls());gc()
options(stringsAsFactors = F)

library(tidyverse)

dat=read.delim("working_data/wgcna/neuron.fc/neuron.gwas.enrich.txt")

for (i in unique(dat$SET)){
  current=dat[dat$SET == i,]
  current$fdr=p.adjust(current$Enrichment_p,method = "fdr")
  pdf(paste("working_data/figures/wgcna.neuron.module.gwas/",i,"_gwas.pdf",sep = ""),width=7.8,height=7.1)
  gg=ggplot(current, aes(x=reorder(GWAS, -log10(Enrichment_p)), y=-log10(fdr))) + geom_bar(stat="identity", fill="red") + 
    coord_flip() + xlab("") + 
    geom_hline(yintercept=-log10(0.05), lty=2, color="black")+
    theme(axis.text.y = element_text(size=12,face = "bold"),axis.text.x = element_text(size = 14),axis.title.x = element_text(size=20))+
    ggtitle("GWAS Enrichment")
  print(gg)
  dev.off()
}
