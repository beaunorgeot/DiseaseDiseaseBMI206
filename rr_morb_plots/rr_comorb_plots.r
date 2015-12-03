# rr co morb plots

# add icd9 to library
library(icd9)
new_s_output <- read.delim("~/Box Sync/UCSF/Fall Classes 2015/BMI_206/Project/DiseaseDiseaseBMI206/kathleen_scripts/new_s_output.csv", comment.char="#")

hudinet <- read.delim("~/Box Sync/UCSF/Fall Classes 2015/BMI_206/Project/DiseaseDiseaseBMI206/data/hudinet_3.txt", header=FALSE, row.names=NULL, colClasses = c("character"))
colnames(hudinet) <- c('ICD9_dis_1', 'ICD9_dis_2', 'Prevalence disease 1', 'Prevalence disease 2', 'Co-occurence between diseases 1 and 2', 'Relative risk', 'Relative Risk 99% Conf. Interval (left)', 'Relative Risk 99% Conf. Interval (right)', 'Phi-correlation', 't-test value')
hudinet<-hudinet[!(hudinet$ICD9_dis_1=="000"),]
hudinet<-hudinet[!(hudinet$ICD9_dis_2=="000"),]
hudinet$`Prevalence disease 1`=NULL
hudinet$`Prevalence disease 2`=NULL
hudinet$`Co-occurence between diseases 1 and 2`=NULL
hudinet$`Relative Risk 99% Conf. Interval (left)`=NULL
hudinet$`Relative Risk 99% Conf. Interval (right)`=NULL
hudinet$`t-test value`=NULL
hudinet$`Phi-correlation`=NULL
hudinet$ICD9_dis_1 <- paste0("'", hudinet$ICD9_dis_1,"'")
hudinet$ICD9_dis_2 <- paste0("'", hudinet$ICD9_dis_2,"'")
hudinet$dis_1 <- icd9Explain(hudinet$ICD9_dis_1)
#hudinet$dis_2 <- lapply(hudinet$ICD9_dis_2, icd9Explain)
#write.csv(hudinet, file = 'hudinet_rev.csv')