# import dataset 2
setwd("~/Box Sync/UCSF/Fall Classes 2015/BMI_206/Project/DiseaseDiseaseBMI206/data")
library(stringr)
DataS2_disease_genes <- read.delim("~/Box Sync/UCSF/Fall Classes 2015/BMI_206/Project/DiseaseDiseaseBMI206/data/DataS2_disease_genes.tsv", comment.char="#")
col_5_disease_genes=sapply(DataS2_disease_genes[,5, drop=FALSE],as.character)
col_5_split=str_split(col_5_disease_genes, ';')
col_6_disease_genes=sapply(DataS2_disease_genes[,6, drop=FALSE],as.character)
col_6_split=str_split(col_5_disease_genes, ';')
View(DataS2_disease_genes)
#View(col_5_split)
