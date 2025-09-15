

##Load in the MCF7 data 
##load the library to bring extra pipe things (must say with a French accent apparently)
library(magrittr)

##load in the MCF7 CE and E2 data. 
fcData = read.table("mcf7_counts.txt", sep="\t", header=TRUE)

##view the top part of the data and then see dimensions and the names of the cols. 
fcData %>% head()
dim(fcData)
names(fcData)

##rename the cols from 7-14 
names(fcData)[7:14] = c("SRS491582E2", "SRS491584E2", "SRS491585E2", "SRS491586E2", "SRS491590CE", "SRS491591CE", "SRS491592CE", "SRS491593CE")
fcData %>% head()

##Remove the annotation cols and rename the rows with the geneIDs
counts = fcData[, 7:14] 
rownames(counts) = fcData$Geneid 
counts %>% head()

##read counts per sample
colSums(counts)

##Change the margins to make it more readable
par(mar= c(7, 4, 4, 2) + 0.1)
##make a barplot. use the counts on the y axis and and cols on the x.
colSums(counts) %>% barplot(., las=3, ylab="Reads mapped per sample")


library(DESeq2)
##Specify the groupings
conds = c("E2","E2","E2", "E2", "CE","CE","CE", "CE")

##Create object of class CountDataSet derived from eSet class
dds = DESeqDataSetFromMatrix(countData = as.matrix(counts), 
                             colData = data.frame(conds=factor(conds)),
                             design = formula(~conds))

##See dimensions of the model - genes by groupings/samples
dim(dds)

##fit the model
dds = DESeq(dds)
res = DESeq2::results(dds)
####DESeq2 package DE analysis
#extract the results from that
knitr::kable(res[1:20,])

## Remove rows with NAs
res = na.omit(res)
## Get the rows of "res" with significant adjusted p-values (default is FDR)
resPadj = res[res$padj <= 0.01 , ]
## Get dimensions, the row number is the sum of significantly differentially expressed genes. 
dim(resPadj)

head (resPadj)
##sort the results in order (such as top table) by fold change. 
resultsSorted <- resPadj[order(abs(resPadj$log2FoldChange), decreasing = TRUE), ]
head(resultsSorted, 20)

