rm(list=ls());gc()
options(stringsAsFactors = F)

library(ggplot2)
library(SummarizedExperiment)


load("working_data/summarizedExperiment/fc_2pass/se_blocks_CPM_outlierRemoved.RData")
load("working_data/summarizedExperiment/fc_2pass/se_neuron_CPM_outlierRemoved.RData")


se_blocks$Diagnosis = factor(se_blocks$Diagnosis, levels=c("Control", "Autism"))
se_neuron$Diagnosis = factor(se_neuron$Diagnosis, levels= c("Control", "Autism"))


# blocks
datMeta=as.data.frame(colData(se_blocks))
X = model.matrix(~Diagnosis + Sex + Age + RIN + X260.280 + seqPC1 + seqPC2 + seqPC3,colData(se_blocks))
Y = assays(se_blocks)$log2CPM
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
#datExpr = Y - t(as.matrix(X[,4]) %*% t(as.matrix(beta[4,])))
datExpr = Y - t(X[,c(5:ncol(X))] %*% beta[c(5:nrow(beta)),])
table(rownames(datMeta) == colnames(datExpr))
save(file = "./working_data/wgcna/voom.forWGCNA.input.blocks.fc.RData",datExpr,datMeta)

# neuron
datMeta=as.data.frame(colData(se_neuron))
X = model.matrix(~Diagnosis + Sex + Age + Type_RNAseqRunNumber + seqPC1 + seqPC2 + seqPC3,colData(se_neuron))
Y = assays(se_neuron)$log2CPM
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
#datExpr = Y - t(as.matrix(X[,5]) %*% t(as.matrix(beta[5,])))
#datExpr = Y - t(X[,c(3:ncol(X))] %*% beta[c(3:nrow(beta)),])
datExpr = assays(se_neuron)$log2CPM
table(rownames(datMeta) == colnames(datExpr))
save(file = "./working_data/context/input_neuron_fc2pass_raw.RData",datExpr,datMeta)


# match with brainspan training data
rm(list=ls());gc()
options(stringsAsFactors = F)

load("./working_data/context/input_neuron_fc2pass_raw.RData")
rownames(datExpr)=sapply(rownames(datExpr),function(x) str_split_fixed(x,"\\.",Inf)[1,1])
link=read.delim("D:/references/idconvert/entrez.ensemble.unambigous.txt")
table(rownames(datExpr) %in% link$ensembl_gene)
datExpr1=as.data.frame(datExpr)
datExpr1$entrezid=link$entrez_gene[match(rownames(datExpr1),link$ensembl_gene)]
datExpr2=datExpr1[complete.cases(datExpr1),]
rownames(datExpr2)=datExpr2$entrezid
datExpr2$entrezid=NULL
datExpr.HNP=datExpr2
rm(datExpr,datExpr1,datExpr2);gc()

load("G:/.shortcut-targets-by-id/1iitsIdiJjzFtQeuRC3ppGzrfRsBh4q2u/context/ToYining/MEKAfiles/annot_1e-05_0.001.Rdata")
annotnew=annotnew[annotnew$ENTREZ_GENE_ID %in% rownames(datExpr.HNP),]
ind = match(annotnew$ENTREZ_GENE_ID,rownames(datExpr.HNP))
datExpr.HNP = datExpr.HNP[ind,]

datExpr = read.csv("G:/.shortcut-targets-by-id/1iitsIdiJjzFtQeuRC3ppGzrfRsBh4q2u/context/ToYining/InVivoData/Kang17565genes1340samples_fromGEO.txt")
annot = read.csv("G:/.shortcut-targets-by-id/1iitsIdiJjzFtQeuRC3ppGzrfRsBh4q2u/context/ToYining/InVivoData/annot.csv")
ind = match(rownames(datExpr),rownames(annot))
annot = annot[ind,]

ind = match(annotnew$ENTREZ_GENE_ID,annot$ENTREZ_GENE_ID);
datExprnew = datExpr[ind,];
datExpr.HNP = datExpr.HNP - median(apply(datExpr.HNP,1,median)) + median(apply(datExprnew,1,median))

ARFFfnamefull = "./working_data/context/input_neuron_raw.arff"



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

### run MEKA classifier in Linux bash
# java -cp '/u/project/gandalm/pampas/context_meka/MEKAfiles/meka-1.3-edit/lib/*' weka.classifiers.multitarget.CR -t /u/project/gandalm/pampas/context_meka/MEKAfiles/Sestan_1e-05_0.001.arff -W weka.classifiers.functions.SMO -T input_neuron.arff > neuron.output
### then plot the output

# plot output
library(WGCNA)
rm(list = ls());gc()
options(stringsAsFactors = F)
load("./working_data/context/input_neuron_fc2pass.RData")
datTraits = datMeta;
#rownames(datTraits)=datTraits$X
nsamples = nrow(datTraits);

## Read and parse distribution file
dist = list(Stage=matrix(NA,nrow=15,ncol=nsamples),Region=matrix(NA,nrow=6,ncol=nsamples));
distunparsed = read.table("working_data/context/neuron_raw.output",sep="}")$V1;

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

datTraits1=datTraits[order(datTraits$Diagnosis),]
dist$Stage=dist$Stage[,rownames(datTraits1)]
dist$Region=dist$Region[,rownames(datTraits1)]

cat('Plotting Stage...\n');
png("working_data/context/Neuron_StageMatrix.png", width=6,height=8,units="in",res=300);
#pdf(paste(base.dir,outputdir,"StageMatrix.pdf",sep="/"), width=6,height=8);
## Create output matrices for Stage and Region
Name = "Temporal Identity of Your Data";

##par(oma=c(1,3,1,1));
labeledHeatmap(Matrix = dist$Stage,
               xLabels = datTraits1$Diagnosis, #with(datTraits2,paste(Group,Lab)),
               yLabels = rownames(dist$Stage),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               setStdMargins = TRUE,
               cex.lab.x = 0.55,
               cex.lab.y = 0.55,  
               main = Name);    
dev.off();

cat('Plotting Region...\n');
png("../RegionMatrix.png", width=6,height=4,units="in",res=300);
#pdf(paste(base.dir,outputdir,"RegionMatrix.pdf",sep="/"), width=6,height=4);
## Create output matrices for Stage and Region
Name = "Regional Identity of Your Data";

##par(oma=c(1,3,1,1));
labeledHeatmap(Matrix = dist$Region[c(3,4,1,6,5,2),],
               xLabels = datTraits1$Diagnosis,#with(datTraits2,paste(Group,Lab)),
               yLabels = rownames(dist$Region[c(3,4,1,6,5,2),]),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               setStdMargins = TRUE,
               cex.lab.x = 0.55,
               cex.lab.y = 0.55,  
               main = Name);    

dev.off();
