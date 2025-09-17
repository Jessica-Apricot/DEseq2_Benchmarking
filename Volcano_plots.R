#Volcano plots
#install.packages(c("systemfonts", "ragg","ggrepel", "plotly", "ggiraph"))
library(ggplot2)
library(ggiraph)
library(plotly)
library(ggrepel)
library(systemfonts)
library(ragg)

#install.packages(c("ggplot2", "grDevices", "Cairo"))

#Define significance thresholds
pval_thres <- 0.05
fc_thres <- 1 #log2 fold

##volcano plot of normal run
#plotting df
volc_res_df <- data.frame(
  logFC = res_df$log2FoldChange,
  negLogPval = -log10(res_df$padj),
  adj.P.Val = res_df$padj,
  ID = rownames(res_df))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_res_df$category <- ifelse(res_df$padj <= pval_thres,
                              ifelse(volc_res_df$logFC >= fc_thres, "Upregulated", 
                                     ifelse(volc_res_df$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                              "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_res_df, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "DESeq2 Data set mcf7 Normal Comparison",
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

#Mcf7 Randomised
#plotting df

volc_resShuffle_df <- data.frame(
  logFC = resShuffle_df$log2FoldChange,
  negLogPval = -log10(resShuffle_df$padj),
  adj.P.Val = resShuffle_df$padj,
  ID = rownames(resShuffle_df))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated

volc_resShuffle_df$category <- ifelse(resShuffle_df$padj <= pval_thres,
                               ifelse(volc_resShuffle_df$logFC >= fc_thres, "Upregulated", 
                                      ifelse(volc_resShuffle_df$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                               "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_resShuffle_df, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "DESeq Data mcf7 Randomised Comparison",
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

#SULT2 Normal
#plotting df
volc_res2_df <- data.frame(
  logFC = res2_df$log2FoldChange,
  negLogPval = -log10(res2_df$padj),
  adj.P.Val = res2_df$padj,
  ID = rownames(res2_df))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_res2_df$category <- ifelse(res2_df$padj <= pval_thres,
                                      ifelse(volc_res2_df$logFC >= fc_thres, "Upregulated", 
                                             ifelse(volc_res2_df$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                      "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_res2_df, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "DESeq2 Data set SULT2 dataset Normal Comparison",
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


#SULT2 Random
#plotting df
volc_res_shuffled_SS_df <- data.frame(
  logFC = res_shuffled_SS_df$log2FoldChange,
  negLogPval = -log10(res_shuffled_SS_df$padj),
  adj.P.Val = res_shuffled_SS_df$padj,
  ID = rownames(res_shuffled_SS_df))

#Column to catagorise genes, if above fc threshold labelled upregulated (greater than 1), if below threshold (-1), downregulated
volc_res_shuffled_SS_df$category <- ifelse(res_shuffled_SS_df$padj <= pval_thres,
                                ifelse(volc_res_shuffled_SS_df$logFC >= fc_thres, "Upregulated", 
                                       ifelse(volc_res_shuffled_SS_df$logFC <= -fc_thres, "Downregulated", "Not Significant")),
                                "Not Significant")

# Create the volcano plot using ggplot2
ggplot(volc_res_shuffled_SS_df, aes(x = logFC, y = negLogPval, color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-fc_thres, fc_thres), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_thres), linetype = "dashed") +
  labs(
    title = "DESeq2 Data set SULT2 Randomised Comparison",
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


