setwd("/Users/beaunorgeot/bmi206/BMI206class/DiseaseDiseaseBMI206")
netdf = read.csv("interactomeNodeTable.csv")
summary(netdf)
# NumberOfUndirectedEdges is an empty column
# PartnerOfMultiEdgedNodePairs is an empty column
# selected is false for all observations

#can consider starting with just 2 clusters, assumming that genes are either targets or not targets
# and further assuming that the features that explain "targetness" are described by networkNode statisitcs

#library(NbClust)
set.seed(1234)

library(dplyr)
#remove empty cols
netdf = netdf %>% select(-c(NumberOfUndirectedEdges,PartnerOfMultiEdgedNodePairs,selected,SUID))
# convert isSingleNode from boolean to binary
netdfClean = netdf %>% mutate(IsSingleNode = ifelse(IsSingleNode == "false",0,1))

#nc <- NbClust(netdfClean, min.nc=2, max.nc=15, method="kmeans")
#fails: The TSS matrix is indefinite. There must be too many missing values. The index cannot be calculated.
#summary(netdfClean)
#no NAs & no negatives
#check if any infinite values
#any( sapply(netdfClean, is.infinite) ) #nothing is infinite. 
#netdfClean = as.matrix(netdfClean)
#nc <- NbClust(netdfClean, min.nc=2, max.nc=15, method="kmeans") #also plans
#is the matrix simply too big?
#netSmall = netdfClean[1:20,]
#ncSmall <- NbClust(netSmall, min.nc=2, max.nc=15, method="kmeans") #size isn't the issue
#this sucks, trying a different package

library(stats)
#kmeans(x, centers, iter.max = 10, nstart = 1, algorithm = c("Hartigan-Wong", "Lloyd", "Forgy","MacQueen"), trace=FALSE)
netdfClean = as.matrix(netdfClean)
# Brief note of concern: Nothing in my data has been centered/scaled
km.out <- kmeans(netdfClean,centers=2,iter.max=5,nstart=25)
#look at cluster assignments
km.out$cluster
summary(km.out)
#plot
plot(netdfClean,col=(km.out$cluster),main = "K-means 2 clusters", pch=20,cex=2)
# this gives me Betweeness centrality vs averageShortestPathLength. How would I view highly dimensional data?

#split the clusters
#group1 = netdfClean[km.out$cluster == 1]
#group2 = netdfClean[km.out$cluster == 2]
#length(group1) #945
#length(group2) #200995
#this is cool, huge imbalance.
dim(netdfClean) #13460 x 15
# WAIT, there were only 13460 observations, where are the big numbers for the groups coming from?
#you can get the size of each cluster in an easier manner
km.out$size # 13397 63 
#ratio's of size of groups and group1/group2 are identical. What's happening with the groups?

#add cluster labels as a column
combdf <- cbind(netdfClean, clusterNum = km.out$cluster)

#Daniel's data: https://raw.githubusercontent.com/dhimmel/integrate/master/compile/CbG-binding.tsv
#library(RCurl)
#data <- getURL("https://raw.githubusercontent.com/dhimmel/integrate/master/compile/CbG-binding.tsv",
 #              ssl.verifypeer=0L, followlocation=1L)
#kdTargets <- read.delim(data, sep = "\t")
# sadly, my attempt to use getURL failed. Will reinvestiage later. Just curled instead.
kdTarget = read.delim("CbG-binding.tsv", header = T)
#now just extract the entrez_gene_id column and compare that column to each group
names(kdTarget)
#Next: create dummy var "in target db". ifelse(0,1)
combdf = as.data.frame(combdf)
combdf = combdf %>% mutate(isTarget = ifelse(name %in% kdTarget$entrez_gene_id,1,0))


#compare clusters to known drug targets using a table
#assuming known drug targets = kdTargets and there is an outcome of 0 or 1 for target or not
tab = table(combdf$isTarget, combdf$clusterNum) #abcd
#oddsratio = ad/bc

#There are 2353 known target genes in cluster 1 and 16 in cluster 2
# there are 11044 non-targets in cluster 1 and 47 in cluster 2.
#13397 63 
#ratios of known targets per cluster: probably nothing interesting here
kt1 = 2353/13397 #0.1756363
kt2 = 16/63 #0.2539683 

#create a vector for each cluster with labels for whether each gene is a known target
group1 = combdf %>% filter(clusterNum == 1) %>% select(isTarget)
group2 = combdf %>% filter(clusterNum == 2) %>% select(isTarget)
#perform a chi-square on the vectors to determine if there is a
chisq.test(tab) #p-value = 0.1435
#perform oddsratio
fisher.test(tab) #odds ratio 1.597748
#oddsratio = ad/bc for col1: odds of being a 0 (known target) are a/c, for col2 odds of being a 0 (known target) are b/d
# odds of being in b/d is 1.6 times greater
# for row 0, odds of being in group1 is a/b
# for row1, odds of bieng in group1 is c/d


#todo
#1.split into 2 or 3 groups (good targets, bad targets, average)
#2.for earch group, do boolean for whether a gene is a known target or not (make sure using same gene
# id system)
#3. do a t-test comparing the 2 groups. do we keep the nullH that membership to a specific cluster
# does not make the gene more likely to be a target? Or do we reject the null (alt hyp is that membership matters)
# 4. Deal w/fact that data values might be on different scales..