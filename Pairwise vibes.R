#AUCOTT STUFF



#HLD Brain 
mcf7_conds <- c("E", "E", "E", "E", "Control", "Control", "Control", "Control", "Random", "Random", "Random", "Random") 

mcf7_metadata <- data.frame(row.names = colnames(H_L_D_Brain), conds = factor(HLD_Br_conds)) HLD_Br_conds_dds <- DESeqDataSetFromMatrix(countData = as.matrix(H_L_D_Brain), colData = HLD_metadata, design = ~ conds) 

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