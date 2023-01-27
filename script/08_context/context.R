# clean count matrix and normalize by limma voom
datMeta=read.delim("./3m_metadata_pc_SVA_nohvg.txt")
gene=read.delim("3mOrg_Expressed_genes_tpm0.5_byGenotype.plus16p.txt",header = F)
gene.remove=read.delim("3mOrg_highly_variable_genes_betweenLabs.txt",header = F)
gene=gene[!gene$V1 %in% gene.remove$V1,]
load("3months_Organoids_RSEM_Quant.genes.counts.RData")
counts5=counts[gene,datMeta$Sample]
dge.voom = voom(calcNormFactors(DGEList(counts = counts5),method = 'TMM'))

datExpr = as.matrix(dge.voom$E)
table(colnames(datExpr)==datMeta$X)
save(datMeta,datExpr,file = "./forContext/3mOrganoid.forContext.exprGenes.fulldata.RData")

# match with brainspan data
load("F:/projects/organoid/rna3month/forContext/3mOrganoid.forContext.exprGenes.fulldata.RData")
link=read.delim("F:/references/idconvert/entrez.ensemble.unambigous.txt")
table(rownames(datExpr) %in% link$ensembl_gene)
datExpr1=as.data.frame(datExpr)
datExpr1$entrezid=link$entrez_gene[match(rownames(datExpr1),link$ensembl_gene)]
datExpr2=datExpr1[complete.cases(datExpr1),]
rownames(datExpr2)=datExpr2$entrezid
datExpr2$entrezid=NULL
datExpr.HNP=datExpr2

load("./MEKAfiles/annot_1e-05_0.001.Rdata")

ind = match(annotnew$ENTREZ_GENE_ID,rownames(datExpr.HNP))
datExpr.HNP = datExpr.HNP[ind,]

datExpr = read.csv("InVivoData/Kang17565genes1340samples_fromGEO.txt")
annot = read.csv("InVivoData/annot.csv")
ind = match(rownames(datExpr),rownames(annot))
annot = annot[ind,]

ind = match(annotnew$ENTREZ_GENE_ID,annot$ENTREZ_GENE_ID);
datExprnew = datExpr[ind,];
datExpr.HNP = datExpr.HNP - median(apply(datExpr.HNP,1,median,na.rm=TRUE),na.rm=TRUE) + median(apply(datExprnew,1,median,na.rm=TRUE));

ARFFfnamefull = "G:/context/YourData.arff"



allstages = separatebycommas(seq(1,15));
allregions = "AMY,CBC,Cortex,HIP,STR,THAL";


for (i in 1:length(annotnew$ENTREZ_GENE_ID)) {
  cat("@ATTRIBUTE",annotnew$ENTREZ_GENE_ID[i],"numeric\n",file=ARFFfnamefull,append=TRUE);
}
cat("\n@DATA\n",file=ARFFfnamefull,append=TRUE);
for (i in 1:ncol(datExpr.HNP)) {
  cat(separatebycommas(c(1,'AMY',t(datExpr.HNP[,i]))),"\n",file=ARFFfnamefull,append=TRUE,sep="");
}

datTraits = read.table("../1mOrganoid//metadata.forContext.txt");
nsamples = nrow(datTraits);

dist = list(Stage=matrix(NA,nrow=15,ncol=nsamples),Region=matrix(NA,nrow=6,ncol=nsamples));
distunparsed = read.table("../YourData.output",sep="}")$V1;

overallcount = 0;
count = 0;
subcount = 0;
subjectcount = 1;
while (overallcount <= nsamples*24) {

  overallcount = overallcount + 1;
  count = count + 1;

  if (count >= 2 & count <= 16) {
    subcount = subcount + 1;
    dist$Stage[subcount,subjectcount] = as.numeric(unlist(strsplit(distunparsed[overallcount],":",fixed=TRUE))[2]);
  } else if (count == 17) {
    subcount = 0;
  } else if (count >= 18 & count <= 23) {
    subcount = subcount + 1;
    dist$Region[subcount,subjectcount] = as.numeric(unlist(strsplit(distunparsed[overallcount],":",fixed=TRUE))[2]);
  } else if (count == 24) {
    count = 0;
    subcount = 0;
    subjectcount = subjectcount + 1;
  }
}

colnames(dist$Stage) = rownames(datTraits);
rownames(dist$Stage) = seq(1,15);
colnames(dist$Region) = rownames(datTraits);
rownames(dist$Region) = c("AMY","CBC","Cortex","HIP","STR","THAL");

datTraits1=datTraits[order(datTraits$Group),]
dist$Stage=dist$Stage[,rownames(datTraits1)]
dist$Region=dist$Region[,rownames(datTraits1)]
png("../StageMatrix.png", width=6,height=8,units="in",res=300);
Name = "Temporal Identity of Your Data";

labeledHeatmap(Matrix = dist$Stage,
               xLabels = datTraits1$Group,
               yLabels = rownames(dist$Stage),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               setStdMargins = TRUE,
               cex.lab.x = 0.55,
               cex.lab.y = 0.55,
               main = Name);
dev.off();

png("../RegionMatrix.png", width=6,height=4,units="in",res=300);

Name = "Regional Identity of Your Data";

labeledHeatmap(Matrix = dist$Region[c(3,4,1,6,5,2),],
               xLabels = datTraits1$Group,
               yLabels = rownames(dist$Region[c(3,4,1,6,5,2),]),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               setStdMargins = TRUE,
               cex.lab.x = 0.55,
               cex.lab.y = 0.55,
               main = Name);

dev.off();
