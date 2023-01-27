rm(list = ls());gc()
options(stringsAsFactors = F)

#load("F:/projects/organoid/rna1month/finallme2/wgcna/voom.forWGCNA.RData")
#load("F:/projects/ipsc/finallme/wgcna/sv1/voom.forWGCNA.RData")
#load("F:/projects/organoid/rna1month/forContext/1mOrganoid.forContext.allGenes.RData")
#load("F:/projects/ipsc/forContext/ipsc.forContext.allGenes.RData")
#load("F:/projects/organoid/rna3month/forContext/3mOrganoid.forContext.allGenes.RData")
#load("F:/projects/organoid/rna3month/forContext/3mOrganoid.forContext.exprGenes.RData")
#load("F:/projects/ipsc/forContext/ipsc.forContext.exprGenes.RData")
#load("F:/projects/organoid/rna1month/forContext/1mOrganoid.forContext.exprGenes.RData")
load("F:/projects/organoid/rna3month/forContext/3mOrganoid.forContext.exprGenes.fulldata.RData")
link=read.delim("F:/references/idconvert/entrez.ensemble.unambigous.txt")
table(rownames(datExpr) %in% link$ensembl_gene)
datExpr1=as.data.frame(datExpr)
datExpr1$entrezid=link$entrez_gene[match(rownames(datExpr1),link$ensembl_gene)]
datExpr2=datExpr1[complete.cases(datExpr1),]
rownames(datExpr2)=datExpr2$entrezid
datExpr2$entrezid=NULL
datExpr.HNP=datExpr2
rm(datExpr,datExpr1,datExpr2);gc()

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

cat("% Your Data\n",file=ARFFfnamefull);
cat("@RELATION 'YourData: -C 2'\n",file=ARFFfnamefull,append=TRUE);
cat("@ATTRIBUTE stage {",allstages,"}\n",file=ARFFfnamefull,append=TRUE,sep="");
cat("@ATTRIBUTE region {",allregions,"}\n",file=ARFFfnamefull,append=TRUE,sep="");
for (i in 1:length(annotnew$ENTREZ_GENE_ID)) {
  cat("@ATTRIBUTE",annotnew$ENTREZ_GENE_ID[i],"numeric\n",file=ARFFfnamefull,append=TRUE);
}
cat("\n@DATA\n",file=ARFFfnamefull,append=TRUE);
for (i in 1:ncol(datExpr.HNP)) {
  cat(separatebycommas(c(1,'AMY',t(datExpr.HNP[,i]))),"\n",file=ARFFfnamefull,append=TRUE,sep="");
}



separatebycommas <- function(x) {
  combo = "";
  for (i in 1:length(x)) {
    if (i < length(x)) {
      if (!is.na(x[i])) {
        combo = paste(combo,x[i],",",sep="");
      } else if (is.na(x[i])) {
        combo = paste(combo,"?",",",sep="");
      }
    } else {
      if (!is.na(x[i])) {
        combo = paste(combo,x[i],sep="");
      } else if (is.na(x[i])) {
        combo = paste(combo,"?",sep="");
      }
    }
  }
  return(combo);
}
