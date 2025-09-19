### Limma Analysis ####
library(limma)
#BiocManager::install("edgeR")
library(edgeR)

mcf7.design = model.matrix(~conds)
head(mcf7.design)

Lm_mcf7 = DGEList(counts=counts)
Lm_mcf7 = calcNormFactors(Lm_mcf7)
keep_mcf7 <- filterByExpr(Lm_mcf7, mcf7.design)
Lm_mcf7 <- Lm_mcf7[keep_mcf7, , keep.lib.sizes=FALSE]

#voom transformation
v.mcf7 <- voom(Lm_mcf7, mcf7.design, plot = TRUE) 
fit.mcf7 <- lmFit(v.mcf7, mcf7.design)
fit.mcf7 <- eBayes(fit.mcf7)

#Extracting results
results_mcf7 <- topTable(fit.mcf7, coef = ncol(mcf7.design), number = Inf, adjust.method = "BH")

#Significant p-values
results_mcf7 <- results_mcf7[results_mcf7$adj.P.Val <= 0.05, ]
dim(results_mcf7) #0 significant genes


##volcano plot of normal run
#plotting df
volc_mcf7_df <- data.frame(
  logFC = results_mcf7$logFC,
  negLogPval = -log10(results_mcf7$adj.P.Val),
  adj.P.Val = results_mcf7$adj.P.Val,
  ID = rownames(results_mcf7))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_mcf7_df$category <- ifelse(results_mcf7$adj.P.Val <= pval_thres,
                               ifelse(volc_mcf7_df$logFC >= fc_thres, "Upregulated", 
                                      ifelse(volc_mcf7_df$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                               "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_mcf7_df, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Limma Data set mcf7 Normal Comparison",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

### Randomised  mcf7 ### 

mcf7_R.design = model.matrix(~conds)
Lm_mcf7_R = DGEList(counts=full_shuffled)
Lm_mcf7_R = calcNormFactors(Lm_mcf7_R)
keep_mcf7_R <- filterByExpr(Lm_mcf7_R, mcf7_R.design)
Lm_mcf7_R <- Lm_mcf7_R[keep_mcf7_R, , keep.lib.sizes=FALSE]

#voom transformation
v.mcf7_R <- voom(Lm_mcf7_R, mcf7_R.design, plot = TRUE) 
fit.mcf7_R <- lmFit(v.mcf7_R, mcf7_R.design)
fit.mcf7_R <- eBayes(fit.mcf7_R)

#Extracting results
results_mcf7_R <- topTable(fit.mcf7_R, coef = ncol(mcf7_R.design), number = Inf, adjust.method = "BH")

#Significant p-values
results_mcf7_R <- results_mcf7_R[results_mcf7_R$adj.P.Val <= 0.05, ]
dim(results_mcf7_R) #0 significant genes


##volcano plot of normal run
#plotting df
volc_mcf7_df_R <- data.frame(
  logFC = results_mcf7_R$logFC,
  negLogPval = -log10(results_mcf7_R$adj.P.Val),
  adj.P.Val = results_mcf7_R$adj.P.Val,
  ID = rownames(results_mcf7_R))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_mcf7_df_R$category <- ifelse(results_mcf7_R$adj.P.Val <= pval_thres,
                                ifelse(volc_mcf7_df_R$logFC >= fc_thres, "Upregulated", 
                                       ifelse(volc_mcf7_df_R$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_mcf7_df_R, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Limma Data set mcf7 Random Comparison",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

## Sult2 Normal

Sult2.design = model.matrix(~conds2)
Lm_Sult2 = DGEList(counts=SULT2_Data)
Lm_Sult2 = calcNormFactors(Lm_Sult2)
keep_sult2 <- filterByExpr(Lm_Sult2, Sult2.design)
Lm_Sult2 <- Lm_Sult2[keep_sult2, , keep.lib.sizes=FALSE]

#voom transformation
v.sult2 <- voom(Lm_Sult2, Sult2.design, plot = TRUE)
fit.sult2 <- lmFit(v.sult2, Sult2.design)
fit.sult2 <- eBayes(fit.sult2)

#Extracting results
results_sult2 <- topTable(fit.sult2, coef = ncol(Sult2.design), number = Inf, adjust.method = "BH")

#Significant p-values
results_sult2 <- results_sult2[results_sult2$adj.P.Val <= 0.05, ]
dim(results_sult2) #0 significant genes


##volcano plot of normal run
#plotting df
volc_sult2_df <- data.frame(
  logFC = results_sult2$logFC,
  negLogPval = -log10(results_sult2$adj.P.Val),
  adj.P.Val = results_sult2$adj.P.Val,
  ID = rownames(results_sult2))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_sult2_df$category <- ifelse(results_sult2$adj.P.Val <= pval_thres,
                                ifelse(volc_sult2_df$logFC >= fc_thres, "Upregulated", 
                                       ifelse(volc_sult2_df$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_sult2_df, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Limma Data set Sult2 Normal Comparison",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

## Sult2 Random ##


sult2_R.design = model.matrix(~conds2)
Lm_sult2_R = DGEList(counts=full_shuffled_SULT2)
Lm_sult2_R = calcNormFactors(Lm_sult2_R)
keep_sult2_R <- filterByExpr(Lm_sult2_R, sult2_R.design)
Lm_sult2_R <- Lm_sult2_R[keep_sult2_R, , keep.lib.sizes=FALSE]

#voom transformation
v.sult2_R <- voom(Lm_sult2_R, sult2_R.design, plot = TRUE) 
fit.sult2_R <- lmFit(v.sult2_R, sult2_R.design)
fit.sult2_R <- eBayes(fit.sult2_R)

#Extracting results
results_sult2_R <- topTable(fit.sult2_R, coef = ncol(sult2_R.design), number = Inf, adjust.method = "BH")

#Significant p-values
results_sult2_R <- results_sult2_R[results_sult2_R$adj.P.Val <= 0.05, ]
dim(results_sult2_R) #0 significant genes


##volcano plot of normal run
#plotting df
volc_sult2_df_R <- data.frame(
  logFC = results_sult2_R$logFC,
  negLogPval = -log10(results_sult2_R$adj.P.Val),
  adj.P.Val = results_sult2_R$adj.P.Val,
  ID = rownames(results_sult2_R))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_sult2_df_R$category <- ifelse(results_sult2_R$adj.P.Val <= pval_thres,
                                  ifelse(volc_sult2_df_R$logFC >= fc_thres, "Upregulated", 
                                         ifelse(volc_sult2_df_R$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                  "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_sult2_df_R, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "Limma Data set Sult2 Random Comparison",
    subtitle = paste("Thresholds: |log2FC| >", fc_thres, "and adjusted p-value <", pval_thres),
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Differential Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )
