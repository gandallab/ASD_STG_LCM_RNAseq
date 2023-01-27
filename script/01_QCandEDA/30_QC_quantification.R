rm(list=ls());gc()
options(stringsAsFactors = F)

salmon.self=read.delim("working_data/salmon/BSS_H06_quant.sf")[,c(1,5)]
salmon.mike=read.delim("working_data/salmon_mike/BSS_H06_mike_quant.sf")[,c(1,5)]
rsem=read.delim("working_data/star_rsem/BSS_H06.rsem.isoforms.results")[,c(1,5)]

library(tidyverse)
salmon.mike$id=sapply(salmon.mike$Name,function(x){str_split_fixed(x,"\\.",Inf)[1,1]})
salmon.self$id=sapply(salmon.self$Name,function(x){str_split_fixed(x,"\\.",Inf)[1,1]})
rsem$id=sapply(rsem$transcript_id,function(x){str_split_fixed(x,"\\.",Inf)[1,1]})

salmon=merge(salmon.mike,salmon.self,by="id")
salmon=salmon[salmon$NumReads.x > 0 | salmon$NumReads.y > 0,]
pick=salmon[salmon$NumReads.x > 0 & salmon$NumReads.x < 1,]
ggplot(salmon,aes(log2(NumReads.x),log2(NumReads.y)))+
  geom_point()

rsem.mike=merge(salmon.mike,rsem,by="id")
rsem.mike=rsem.mike[rsem.mike$NumReads > 0 | rsem.mike$expected_count >0,]

ggplot(rsem.mike,aes(log2(NumReads),log2(expected_count)))+
  geom_point()

rsem.self=merge(salmon.self,rsem,by="id")
rsem.self=rsem.self[rsem.self$NumReads > 0 | rsem.self$expected_count >0,]

ggplot(rsem.self,aes(log2(NumReads),log2(expected_count)))+
  geom_point()

