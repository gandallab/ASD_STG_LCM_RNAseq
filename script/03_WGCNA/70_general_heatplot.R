options(stringsAsFactors = F)
library(tidyverse)
library(WGCNA)
library(ggdendro)

rm(list=ls());gc()

to.keep=c("M1","M3","M10","M14","M19","M31","M32")

load("working_data/wgcna/blocks.fc/blocks.fc.recut.RData")

MEs=networks$MEs
me=MEs$eigengenes[,-1]
colnames(me)=gsub("ME","M",colnames(me))
me=me[,colnames(me) %in% to.keep]

hc=hclust(as.dist(1 - cor(me)))
#hc=hclust(dist(t(me)))
plot(hc,hang = -1)
dend=as.dendrogram(hc)
dend_data <- dendro_data(dend)

# heatmap1

# gene_pos_table = data.frame(x_center = seq(1,length(to.keep)),moduleNumber=to.keep,width=1)
gene_pos_table <- with(
  dend_data$labels, 
  data.frame(x_center = x, moduleNumber = as.character(label), width = 1))

sample_pos_table <- data.frame(Group = c("ASD")) %>%
  dplyr::mutate(y_center = (1:n()), 
                height = 1)

gene_axis_limits <- with(
  gene_pos_table, 
  c(min(x_center - 0.5 * width), max(x_center + 0.5 * width))
) + 
  0.1 * c(-1, 1)

heatdata=read.delim("working_data/wgcna/blocks.fc/modTrait.blocks.fc.txt")
heatdata$text[!is.na(heatdata$text)]="*"
heatdata$moduleNumber=paste("M",heatdata$moduleNumber,sep = "")

heatdata=heatdata[heatdata$moduleNumber %in% to.keep,]
heatdata$Group = "ASD"
heatmap_data <- heatdata %>% 
  left_join(gene_pos_table) %>%
  left_join(sample_pos_table)


plt_hmap1 <- ggplot(heatmap_data, 
                       aes(x = x_center, y = y_center, fill = beta, 
                           height = height, width = width,label=text)) + 
  geom_tile(aes(fill = beta),color="grey60") +
  geom_text(aes(y=y_center - 0.3),size=10)+
  scale_fill_gradient2("beta", high = "#EF8961",mid="white", low = "#8FBDDA") +
  scale_x_continuous(expand = c(0,0),
                     limits = gene_axis_limits,
                     breaks = gene_pos_table$x_center,
                     sec.axis = dup_axis())+
  scale_y_continuous(breaks = sample_pos_table$y_center, 
                     labels = sample_pos_table$Group,
                     #limits = c(0,2.5),
                     expand = c(0,0)) + 
  labs(x = "", y = "")+ 
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  coord_equal()+
  guides(fill=F)


# heatmap2
load("working_data/wgcna/blocks.fc/ewce.block.RData")
heatdata=adult_all
heatdata$fdr[heatdata$fdr == 0]=min(heatdata$fdr[heatdata$fdr > 0])
heatdata$text=""
heatdata$text[heatdata$fdr < 0.05]="*"
heatdata$moduleNumber=paste("M",heatdata$moduleNumber,sep = "")
heatdata=heatdata[heatdata$moduleNumber %in% to.keep,]


#gene_pos_table = data.frame(x_center = seq(1,length(to.keep)),moduleNumber=to.keep,width=1)
gene_pos_table <- with(
  dend_data$labels, 
  data.frame(x_center = x, moduleNumber = as.character(label), width = 1))

sample_pos_table <- data.frame(CellType = unique(heatdata$CellType)) %>%
  dplyr::mutate(y_center = (1:n()), 
                height = 1)

gene_axis_limits <- with(
  gene_pos_table, 
  c(min(x_center - 0.5 * width), max(x_center + 0.5 * width))
) + 
  0.1 * c(-1, 1)







heatmap_data <- heatdata %>% 
  left_join(gene_pos_table) %>%
  left_join(sample_pos_table)


plt_hmap2 <- ggplot(heatmap_data, 
                    aes(x = x_center, y = y_center, fill = -log10(fdr), 
                        height = height, width = width,label=text)) + 
  geom_tile(color="grey60") +
  geom_text(aes(y=y_center - 0.3),size=10)+
  scale_fill_gradient2("-log10FDR", low ="white", mid="white", high = "darkgreen",midpoint = 1.3) +
  scale_x_continuous(expand = c(0,0),
                     limits = gene_axis_limits,
                     breaks = gene_pos_table$x_center,
                     sec.axis = dup_axis())+
  scale_y_continuous(breaks = sample_pos_table$y_center, 
                     labels = sample_pos_table$CellType,
                     expand = c(0,0)) + 
  labs(x = "", y = "")+ 
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  coord_equal()+
  guides(fill=F)


# heatmap3
load("working_data/listEnrich/permutTable.1000times.blockModule.newhouseList.RData")
heatdata=permTable
heatdata$text[!is.na(heatdata$text) & heatdata$text != ""]="*"
heatdata$moduleNumber=paste("M",heatdata$moduleNumber,sep = "")
heatdata=heatdata[heatdata$moduleNumber %in% to.keep,]
# heatdata=heatdata[heatdata$category %in% c("Sfari_S_1_2",
#                                            "presynaptic",
#                                            "postsynaptic",
#                                            "pLI_099",
#                                            "HighlyMutationalConstraintGenes_SamochaNatGen2014",
#                                            "FMRP_bindingTarget_Darnell 2011",
#                                            "CHD8_targets"),]

#gene_pos_table = data.frame(x_center = seq(1,length(to.keep)),moduleNumber=to.keep,width=1)

gene_pos_table <- with(
  dend_data$labels, 
  data.frame(x_center = x, moduleNumber = as.character(label), width = 1))

sample_pos_table <- data.frame(category = unique(heatdata$category)) %>%
  dplyr::mutate(y_center = (1:n()), 
                height = 1)

gene_axis_limits <- with(
  gene_pos_table, 
  c(min(x_center - 0.5 * width), max(x_center + 0.5 * width))
) + 
  0.1 * c(-1, 1)





heatmap_data <- heatdata %>% 
  left_join(gene_pos_table) %>%
  left_join(sample_pos_table)


plt_hmap3 <- ggplot(heatmap_data, 
                       aes(x = x_center, y = y_center, fill = signedLog10fdr, 
                           height = height, width = width,label=text)) + 
  geom_tile(color="grey60") +
  geom_text(aes(y=y_center - 0.3),size=10)+
  scale_fill_gradient2("-log10FDR", high = "purple",low = "white") +
  scale_x_continuous(expand = c(0,0),
                     limits = gene_axis_limits,
                     breaks = gene_pos_table$x_center,
                     sec.axis = dup_axis())+
  scale_y_continuous(breaks = sample_pos_table$y_center, 
                     labels = sample_pos_table$category,
                     expand = c(0,0)) + 
  labs(x = "", y = "")+ 
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  coord_equal()+
  guides(fill=F)


# heatmap4



link=data.frame(moduleNumber=seq(1,50),moduleColor=labels2colors(seq(1,50)))
heatdata=read.delim("working_data/wgcna/blocks.fc/block.gwas.enrich.txt")
heatdata=heatdata[heatdata$GWAS == "ASD.iPSYCHPGC.2018",]
heatdata$fdr=p.adjust(heatdata$Enrichment_p,method = "fdr")
heatdata$text=""
heatdata$text[heatdata$fdr < 0.1]="*"
heatdata$moduleNumber=link$moduleNumber[match(heatdata$SET,link$moduleColor)]
heatdata$moduleNumber=paste("M",heatdata$moduleNumber,sep = "")

heatdata=heatdata[heatdata$moduleNumber %in% to.keep,]
heatdata=heatdata[,c("GWAS","fdr" , "text"  , "moduleNumber")]
heatdata$GWAS = "ASD_GWAS"
heatdata=rbind(heatdata,data.frame(GWAS = "ASD_GWAS", fdr =1, text = "", moduleNumber="M3"))

#gene_pos_table = data.frame(x_center = seq(1,length(to.keep)),moduleNumber=to.keep,width=1)

gene_pos_table <- with(
  dend_data$labels, 
  data.frame(x_center = x, moduleNumber = as.character(label), width = 1))

sample_pos_table <- data.frame(GWAS = c("ASD_GWAS")) %>%
  dplyr::mutate(y_center = (1:n()), 
                height = 1)

gene_axis_limits <- with(
  gene_pos_table, 
  c(min(x_center - 0.5 * width), max(x_center + 0.5 * width))
) + 
  0.1 * c(-1, 1)

heatmap_data <- heatdata %>% 
  left_join(gene_pos_table) %>%
  left_join(sample_pos_table)


plt_hmap4 <- ggplot(heatmap_data, 
                    aes(x = x_center, y = y_center, fill = -log10(fdr), 
                        height = height, width = width,label=text)) + 
  geom_tile(aes(fill = -log10(fdr)),color="grey60") +
  geom_text(aes(y=y_center - 0.3),size=10)+
  scale_fill_gradient( high = "pink",low = "white") +
  scale_x_continuous(expand = c(0,0),
                     limits = gene_axis_limits,
                     breaks = gene_pos_table$x_center,
                     sec.axis = dup_axis())+
  scale_y_continuous(breaks = sample_pos_table$y_center, 
                     labels = sample_pos_table$GWAS,
                     #limits = c(0,2.5),
                     expand = c(0,0)) + 
  labs(x = "", y = "")+ 
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  coord_equal()+
  guides(fill=F)


# plot dendrogram
dend_label=label(dend_data)
dend_label$color=lapply(dend_label$label,function(x){
  labels2colors(as.numeric(gsub("M","",x)))
})
plt_dendr <- ggplot(segment(dend_data)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_continuous(expand = c(0, 0),
                     limits = gene_axis_limits) + 
  scale_y_continuous(expand = c(0, 1)) + 
  labs(x = "", y = "") +
  geom_point(data=dend_label,aes(x=x,y=y),color=dend_label$color,size=5)+
  geom_text(data=dend_label,aes(x=x,y=y,label=label),angle=-90,hjust=-0.3,size=6)+
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

library(cowplot)
library(Cairo)


cairo_pdf("working_data/final/panels/blocks.module.heatplot.070520.pdf",width = 7,height = 20,family = "Arial")
# plot_grid(plt_hmap1, plt_hmap4, plt_hmap2, plt_hmap3,  
#           ncol = 1,align = "v",rel_heights = c(2.2,2.2,8,11))

grid.draw(rbind(ggplotGrob(plt_dendr), 
                ggplotGrob(plt_hmap1), 
                ggplotGrob(plt_hmap4),
                ggplotGrob(plt_hmap2),
                ggplotGrob(plt_hmap3),
                size = "last"))

dev.off()
save(file="../co.expr.general/3mOrg.RData",m3_plt_dendr, m3_plt_hmap1, m3_plt_hmap3)



