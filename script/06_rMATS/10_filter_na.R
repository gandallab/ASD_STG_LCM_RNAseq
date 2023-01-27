rm(list = ls());gc()
options(stringsAsFactors = F)

library(stringr)

juncs=list.files("./working_data/rmats/block.out/","*MATS.JC.txt")
alldat=data.frame()
for (i in juncs){
  event=str_split_fixed(i,"\\.",Inf)[1,1]
  if (event == "MXE"){
    startcol=23
  }else{
    startcol=21
  }
  dat=read.delim(paste("working_data/rmats/block.out/",i,sep = ""))
  ind=apply(dat,1,function(x){
    a=str_split_fixed(x[startcol],",",Inf)
    b=str_split_fixed(x[startcol+1],",",Inf)
    y=sum(is.na(as.numeric(a)))/length(a)
    z=sum(is.na(as.numeric(b)))/length(b)
    (y <= 0.2 & z <= 0.2)
  })
  currentdat=dat[ind,c("ID","GeneID","geneSymbol","chr","strand","IncLevel1","IncLevel2","IncLevelDifference","PValue")]
  currentdat$event=event
  alldat=rbind(alldat,currentdat)
}
table(alldat$event)
write.table(alldat,file = "./working_data/rmats/block.out/filter.combined.events.block.txt",
            quote = F,sep = "\t",row.names = F)

library(ggplot2)

alldat=read.delim("working_data/rmats/block.out/filter.combined.events.block.txt")
dat=alldat[alldat$PValue < 0.05,]
ggplot(dat,aes(event))+geom_bar()

plotdata=read.delim("working_data/rmats/event.by.datatype.filtered.txt")
ggplot(plotdata,aes(Event,Counts,fill=Datatype))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_y_continuous(trans = "log1p",breaks = c(10,50,100,500,1000,3000,5000,20000,40000))+
  theme(text = element_text(family = "Arial",size=15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))+
  ggtitle("rMATS Events Count")
