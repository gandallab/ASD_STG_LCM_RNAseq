library(tidyverse)
options(stringsAsFactors = F)

rm(list=ls());gc()

dat=read.delim("working_data/revision/IPA/Dopamine_RegulatorsToTargets_NeuronsDxAdjP0.1_11-14-22bss.txt")
dat1=dat[,1:3]

dat2 = dat[,-c(2,3)] %>%
  gather("regulator","effect",2:18) %>%
  filter(effect != " ")

dat21=dat2 %>%
  group_by(Target) %>%
  summarise(n=n()) %>%
  filter(n > 2)

dat10=dat1[dat1$Target %in% dat21$Target,]
dat11=data.frame(Target = colnames(dat)[4:20],Expr.Log.Ratio="", Molecule.Type="regulator")
dat12=rbind(dat10,dat11)

dat22=dat2[dat2$Target %in% dat21$Target,]

dat3=read.delim("working_data/revision/IPA/Dopamine_Neurons-DxAdjP0.1_MechNetworp0.05p0.005_AllRelationships_11-14-22bss.txt")
dat3$From.Molecule.s.=make.names(dat3$From.Molecule.s.)
dat3$To.Molecule.s.=make.names(dat3$To.Molecule.s.)
dat31=dat3 %>%
  group_by(To.Molecule.s. , From.Molecule.s.) %>%
  summarise(relationship=paste(Relationship.Type,collapse = ";"))
dat32=dat31[!dat31$To.Molecule.s. %in% setdiff(dat1$Target,dat10$Target),]
dat4=merge(dat22,dat32,by.x=c("Target","regulator"),by.y=c("To.Molecule.s.","From.Molecule.s."),all.y=TRUE)

#relationship=as.data.frame(table(dat3$Relationship.Type))

dat4$relationship2=sapply(dat4$relationship,function(x){
  if (length(intersect(x,c("chemical-protein interactions","modification","molecular cleavage","phosphorylation","protein-DNA interactions","ubiquitination")))!=0){
    a="direct"
  }else{
    a="indirect"
  }
  return(a)
})

dat4$effect=apply(dat4,1,function(x){
  if (is.na(x[3]) & grepl("activation",x[4])){
    a="Activated"
  }else{
    a=x[3]
  }
  return(a)
})

write.table(dat12,file = "working_data/revision/IPA/processed_node_table.txt",quote = F,sep = "\t",row.names = F)
write.table(dat4,file = "working_data/revision/IPA/processed_edge_table.txt",quote = F,sep = "\t",row.names = F)
