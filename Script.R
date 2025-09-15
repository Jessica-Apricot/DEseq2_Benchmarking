
######code from stats
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

<<<<<<< HEAD





=======
>>>>>>> 55825e64ae8bebbce5d5a6e79f8e2449871f9ba3







###randomise the gene names
##select out the treatment ones. Only randomise these
shuffle_samples <- c("SRS491582E2", "SRS491584E2", 
                     "SRS491585E2", "SRS491586E2")

##extract them
subset_counts <- counts[, shuffle_samples]

##randomise their gene names (idk why this code works but it does)
set.seed(123)
shuffled_subset <- subset_counts
rownames(shuffled_subset) <- sample(rownames(subset_counts))
##sanity check. make sure they are different.
head(rownames(counts))          #original order
head(rownames(shuffled_subset)) #after shuffling

##now reinsert into the full matrix
full_shuffled <- counts #for new matrix with og control data too
full_shuffled[, shuffle_samples] <- shuffled_subset[rownames(full_shuffled), ]

library(DESeq2)

# colData from before (Treatment vs Control)
samples <- colnames(full_shuffled)
condition <- ifelse(grepl("E2$", samples), "Treatment", "Control")
colData <- data.frame(row.names = samples,
                      condition = condition)

dds <- DESeqDataSetFromMatrix(countData = full_shuffled,
                              colData   = colData,
                              design    = ~ condition)

dds <- DESeq(dds)
res_shuffled <- results(dds)
res_shuffled = na.omit(res_shuffled)
resPadj_shuffled = res_shuffled[res_shuffled$padj <= 0.01 , ]
head(resultsSorted, 20)
dim(resPadj_shuffled)


##set data for plot
res_df <- as.data.frame(res)
resShuffle_df <- as.data.frame(res_shuffled)

##for significnace
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "yes", "no")
resShuffle_df$significant <- ifelse(resShuffle_df$padj < 0.05 & abs(resShuffle_df$log2FoldChange) > 1, "yes", "no")

##volcano plot of normal run
library("ggplot2")
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 p-value") +
  ggtitle("Volcano Plot")

##volcano plot of shuffled run (for paul)
library("ggplot2")
ggplot(resShuffle_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 p-value") +
  ggtitle("Volcano Plot")



#install.packages("openxlsx")  
library(openxlsx)
library(DESeq2)


SULT2_Data <- read.xlsx("GSE129746_counts.xlsx", sheet = 1)


#Rename
rownames(SULT2_Data) <- SULT2_Data[[1]]
SULT2_Data <- SULT2_Data[, -1]
SULT2_Data <- SULT2_Data[,1:8]

##Specify the groupings
conds2 = c("M","M","M", "M", "Mock","Mock","Mock", "Mock")

##Create object of class CountDataSet derived from eSet class
dds2 <- DESeqDataSetFromMatrix(
  countData = as.matrix(SULT2_Data),
  colData = data.frame(conds = factor(conds2)),
  design = ~conds
)

##fit the model
dds2 = DESeq(dds2)
res2 = DESeq2::results(dds2)

res2 = na.omit(res2)
## Get the rows of "res" with significant adjusted p-values (default is FDR)
res2Padj = res2[res2$padj <= 0.01 , ]

results2Sorted <- res2Padj[order(abs(res2Padj$log2FoldChange), decreasing = TRUE), ]
head(results2Sorted, 20)



