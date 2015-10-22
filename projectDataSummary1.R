setwd("/Users/beaunorgeot/BMI206/DiseaseDiseaseBMI206/data")

#load data
indata1 <- read.table("DataS1_interactome.tsv", sep="\t", header = T)
#indata2 <- read.delim("DataS2_disease_genes.tsv", sep="\t",header = T)
# data2 has funky format, cols 5&6 are comma seperated. Fix: load cols 1-4 as tsv, load 5&6 as csv, col bind the 2
indata3 <- read.table("DataS3_localization.tsv", sep="\t", header = T)
indata4 <- read.table("DataS4_disease_pairs.tsv", sep="\t", header = T)

#assign names
names(indata1) <- c("gene_ID_1","gene_ID_2","data_source")
#names(indata2) <- c("disease","number_of_all_genes","number_of_OMIM_genes","number_of_GWAS_genes","OMIM_genes","GWAS_genes")
names(indata3) <- c("disease","number_of_all_genes","LCC_size","LCC_z_score","d_s","d_s_Glass_delta", "d_s_p_value")
names(indata4) <- c("disease_A","disease_B", "s_AB(observed)" ,"d_AB(observed)" ,"z(full rand)"
,"p(full rand)","q(full rand)","z(MeSH rand)","p(MeSH rand)","q(MeSH rand)" )

dataSummary1 <- summary(indata1)
#dataSummary2 <- summary(indata2)
dataSummary3 <- summary(indata3)
dataSummary4 <- summary(indata4)

save(dataSummary1,file = "dataSummary1.Rdata")
#save(dataSummary2,file = "dataSummary2.Rdata")
save(dataSummary3,file = "dataSummary3.Rdata")
save(dataSummary4,file = "dataSummary4.Rdata")

