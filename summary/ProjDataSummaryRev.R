# revised stuff for a new summary because I like full points woo
library(stringr)

# data set 1
DataS1_interactome <- read.table("~/Box Sync/UCSF/Fall Classes 2015/BMI_206/Project/DiseaseDiseaseBMI206/data/DataS1_interactome.tsv", header=TRUE, quote="\"")

#nrow(DataS1_interactome)

col_3_interactome=sapply(DataS1_interactome[,3, drop=FALSE],as.character) # select third column (types of interactions)
#head(col_3_interactome,5)

# count # of ; per entry to determine # of interaction types


#col_3_split=strsplit(col_3_interactome,';')

col_3_split=str_split_fixed(col_3_interactome, ';',7)
#head(col_3_split,5)

# ignoring data set 2 again (wtf is up with that format)

# data set 3
DataS3_localization <- read.delim("~/Box Sync/UCSF/Fall Classes 2015/BMI_206/Project/DiseaseDiseaseBMI206/data/DataS3_localization.tsv", comment.char="#")
dataSummary3 <- summary(DataS3_localization)
table_sum=table(DataS3_localization$disease,DataS3_localization$number_of_all_genes)
head(table_sum,10)
mean(DataS3_localization$number_of_all_genes)

# data set 4

DataS4_disease_pairs <- read.delim("~/Box Sync/UCSF/Fall Classes 2015/BMI_206/Project/DiseaseDiseaseBMI206/data/DataS4_disease_pairs.tsv", comment.char="#")
summ_4=summary(DataS4_disease_pairs)
table(DataS4_disease_pairs$disease_A,DataS4_disease_pairs$disease_B)
nrow(DataS4_disease_pairs)