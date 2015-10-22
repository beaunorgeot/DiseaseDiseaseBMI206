setwd("/Users/beaunorgeot/BMI206/DiseaseDiseaseBMI206/data")

indata <- read.table("DataS1_interactome.tsv", sep="\t", header = T)

names(indata) <- c("gene_ID_1","gene_ID_2","data_source")

dataSummary <- summary(indata)

save(dataSummary,file = "dataSummary.Rdata")