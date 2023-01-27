options(stringsAsFactors = F)
library(xlsx)
library(utils)
library(stringr)
rm(list = ls());gc()

#block

goes=list.files("D:\\projects\\ASD_STG_LCM_RNAseq\\working_data\\wgcna\\blocks.fc\\go\\",glob2rx("go*txt"))
outwb=createWorkbook()

for (go in goes){
  current=read.delim(paste("D:\\projects\\ASD_STG_LCM_RNAseq\\working_data\\wgcna\\blocks.fc\\go\\",go,sep=""))
  sheet=createSheet(outwb,sheetName = str_sub(go,4,-5))
  addDataFrame(current,sheet,row.names = F)
}
saveWorkbook(outwb,file = "./working_data/final/tables_v1/Supplementary Table 3_WGCNA_block_tissue_GO.xlsx")

# neuron

goes=list.files("D:/projects/ASD_STG_LCM_RNAseq/working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/go/",glob2rx("go*txt"))
outwb=createWorkbook()

for (go in goes){
  current=read.delim(paste("D:/projects/ASD_STG_LCM_RNAseq/working_data/wgcna/neuron.fc2pass/neuron_fc2pass_ds2mod100/go/",go,sep=""))
  sheet=createSheet(outwb,sheetName = str_sub(go,4,-5))
  addDataFrame(current,sheet,row.names = F)
}
saveWorkbook(outwb,file = "./working_data/final/tables_v1/Supplementary Table 8_WGCNA_LCM_neuron_GO.xlsx")