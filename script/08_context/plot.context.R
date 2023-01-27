rm(list = ls());gc()
options(stringsAsFactors = F)
library(WGCNA);
datTraits = read.table("../1mOrganoid//metadata.forContext.txt");
nsamples = nrow(datTraits);

## Read and parse distribution file
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
cat('Plotting Stage...\n');
png("../StageMatrix.png", width=6,height=8,units="in",res=300);
#pdf(paste(base.dir,outputdir,"StageMatrix.pdf",sep="/"), width=6,height=8);
## Create output matrices for Stage and Region
Name = "Temporal Identity of Your Data";

##par(oma=c(1,3,1,1));
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

cat('Plotting Region...\n');
png("../RegionMatrix.png", width=6,height=4,units="in",res=300);
#pdf(paste(base.dir,outputdir,"RegionMatrix.pdf",sep="/"), width=6,height=4);
## Create output matrices for Stage and Region
Name = "Regional Identity of Your Data";

##par(oma=c(1,3,1,1));
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
