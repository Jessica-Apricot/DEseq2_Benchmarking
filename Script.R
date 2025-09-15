

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

#AUCOTT STUFFs

#HLD Brain HLD_Br_conds <- c("DMSO", "Low", "DMSO", "DMSO", "Low", "Low", "High", "High", "High", "High", "DMSO") 

HLD_metadata <- data.frame(row.names = colnames(H_L_D_Brain), conds = factor(HLD_Br_conds)) HLD_Br_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(H_L_D_Brain), colData = HLD_metadata, design = ~ conds) 

keep <- rowSums(cpm(counts(HLD_Br_conds_dds)) >= 1) >= 6 HLD_Br_conds_dds <- HLD_Br_conds_dds[keep,] 

HLD_Br_conds_dds <- DESeq(HLD_Br_conds_dds) 

#Pairwise 

#High vs Low 

Br_res_HL <- results(HLD_Br_conds_dds, contrast = c("conds", "High", "Low")) 

#High vs DMSO 

Br_res_HD <- results(HLD_Br_conds_dds, contrast = c("conds", "High", "DMSO")) 

#Low vs DMSO 

Br_res_LD <- results(HLD_Br_conds_dds, contrast = c("conds", "Low", "DMSO")) 

#remove na values  

Br_res_HL <- na.omit(Br_res_HL) Br_res_HD <- na.omit(Br_res_HD) Br_res_LD <- na.omit(Br_res_LD) #significant p-values Br_res_s_HL <- Br_res_HL[ Br_res_HL$padj <= 0.05 & (Br_res_HL$log2FoldChange > 1 | Br_res_HL$log2FoldChange < -1), ] dim(Br_res_s_HL)  

Br_res_s_HD <- Br_res_HD[ Br_res_HD$padj <= 0.05 & (Br_res_HD$log2FoldChange > 1 | Br_res_HD$log2FoldChange < -1), ] dim(Br_res_s_HD)  

Br_res_s_HD Br_res_s_LD <- Br_res_LD[ Br_res_LD$padj <= 0.05 & (Br_res_LD$log2FoldChange > 1 | Br_res_LD$log2FoldChange < -1), ] dim(Br_res_s_LD)  
