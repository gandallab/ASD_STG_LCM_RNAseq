rm(list = ls());gc()
options(stringsAsFactors = F)
library(tidyverse)

dat=read.delim("working_data/revision/donor_info_toPlot.txt")
dat$Diagnosis=factor(dat$Diagnosis,levels = c("Control","Autism"))
dat$Brain_Bank_Source=str_trim(dat$Brain_Bank_Source)
dat$PMI=as.numeric(dat$PMI)

ggplot(dat)+
  geom_bar(aes(Diagnosis,fill=Diagnosis),width = 0.5)+
  xlab("")+
  theme_bw()+
  theme(text = element_text(family = "Arial",size=19,face = "bold"),legend.position = "none")

ggplot(dat)+
  geom_boxplot(aes(Diagnosis,Age,fill=Diagnosis),width = 0.5)+
  xlab("")+
  theme_bw()+
  theme(text = element_text(family = "Arial",size=19,face = "bold"),legend.position = "none")
  

dat %>% 
  group_by(Diagnosis,Gender) %>%
  summarise(n=n()) %>%
  ggplot()+
  geom_bar(aes(Diagnosis,n,fill=Gender),width = 0.5,stat = "identity",position = "fill")+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(text = element_text(family = "Arial",size=19,face = "bold"),legend.position = "none")

dat %>%
  group_by(Brain_Bank_Source,Diagnosis) %>%
  summarise(n=n()) %>%
  ggplot()+
  geom_bar(aes(Brain_Bank_Source,n,fill=Diagnosis),width = 0.5,stat = "identity",position = "stack")+
  xlab("")+
  ylab("count")+
  theme_bw()+
  theme(text = element_text(family = "Arial",size=19,face = "bold"),legend.position = "none",axis.text.x = element_text(size=14,angle = 90,vjust = 0.5, hjust=1))

ggplot(dat)+
  geom_boxplot(aes(Diagnosis,PMI,fill=Diagnosis),width = 0.5)+
  xlab("")+
  ylab("Hour")+
  theme_bw()+
  theme(text = element_text(family = "Arial",size=19,face = "bold"),legend.position = "none")

ggplot(dat)+
  geom_boxplot(aes(Diagnosis,RIN,fill=Diagnosis),width = 0.5)+
  xlab("")+
  scale_y_continuous(limits = c(2,8))+
  theme_bw()+
  theme(text = element_text(family = "Arial",size=19,face = "bold"),legend.position = "none")
plotdetail=ggplot_build(p)
