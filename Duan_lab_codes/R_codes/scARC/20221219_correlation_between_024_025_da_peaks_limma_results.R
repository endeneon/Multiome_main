# Chuxuan Li 12/19/2022
# Script adapted from 14Dec2022_correlation_between_limma_024_025.R,
#Correlate the p-value and logFC of 025 and 024 limma results - ATACseq

# init ####
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)

paths_025 <- list.files(path = "./dapeaks_limma_doover/unfiltered_by_padj", 
                        pattern = "[1|6]v0", full.names = T)
res_025 <- vector("list", length(paths_025))
for (i in 1:length(paths_025)) {
  res_025[[i]] <- as.data.frame(read_csv(paths_025[i]))[,2:8]
}
names(res_025) <- str_replace(str_extract(paths_025, "[A-Za-z]+_res_[1|6]v0"), "_res_", " ")

paths_024 <- list.files(path = "../Duan_project_024_ATAC/da_peaks_limma/unfiltered_by_padj", 
                        pattern = "[1|6]v0", full.names = T)
res_024 <- vector("list", length(paths_024))
for (i in 1:length(paths_024)) {
  res_024[[i]] <- as.data.frame(read_csv(paths_024[i]))[,2:8]
}
names(res_024) <- str_replace(str_extract(paths_024, "[A-Za-z]+_res_[1|6]v0"), "_res_", " ")

# find intersect gene list ####
peaks <- vector("list", length(res_024))
gene_counts_table <- matrix(nrow = 3, ncol = 3, 
                            dimnames = list(c("024", "025", "shared"), 
                                            str_remove(names(res_024)[c(1,3,5)], " 1v0")))
for (i in 1:length(res_024)) {
  peaks[[i]] <- intersect(res_024[[i]]$gene, res_025[[i]]$gene)
  if (i %in% c(1,3,5)) {
    gene_counts_table[1, i %/% 2 + 1] <- length(res_024[[i]]$gene)
    gene_counts_table[2, i %/% 2 + 1] <- length(res_025[[i]]$gene)
    gene_counts_table[3, i %/% 2 + 1] <- length(peaks[[i]])
  }
  res_024[[i]] <- res_024[[i]][res_024[[i]]$gene %in% peaks[[i]], ]
  res_025[[i]] <- res_025[[i]][res_025[[i]]$gene %in% peaks[[i]], ]
  
  res_024[[i]]$neglogp <- (-1) * log10(res_024[[i]]$P.Value)
  res_025[[i]]$neglogp <- (-1) * log10(res_025[[i]]$P.Value)
}
write.table(gene_counts_table, 
            file = "./Analysis_part2_GRCh38/024_025_limma_correlation_plots/shared_high_expression_gene_counts.csv",
            quote = F, sep = ",", row.names = T, col.names = T)
# test = res_024[[i]] %>% left_join(res_025[[i]], by = "gene")

# all peaks shared between 024 and 025 ####
for (i in 1:length(res_024)) {
  data <- res_024[[i]] %>% left_join(res_025[[i]], by = "gene")
  # correlation between p-values
  png(paste0("./Analysis_part2_GRCh38/024_025_limma_correlation_plots/all_shared_peaks/", 
             str_replace(names(res_024)[i], " ", "_"), "_all_gene_pval_corr.png"))
  p1 <- ggscatter(data = data,
                  x = "P.Value.x", y = "P.Value.y", size = 0.01, add = "reg.line", conf.int = T, 
                  cor.coef = T, cor.method = "pearson", xlab = "024 limma p-value",
                  ylab = "025 limma p-value", xlim = c(0, 1), ylim = c(0, 1), 
                  title = paste0("Correlation between 024 and 025 ", names(res_024)[i], "\nlimma differential gene expression results - p value"))
  print(p1)
  dev.off()
  # correlation between FDR
  png(paste0("./Analysis_part2_GRCh38/024_025_limma_correlation_plots/all_shared_peaks/", 
             str_replace(names(res_024)[i], " ", "_"), "_all_gene_fdr_corr.png"))
  p2 <- ggscatter(data = data,
                  x = "adj.P.Val.x", y = "adj.P.Val.y", size = 0.01, add = "reg.line", conf.int = T, 
                  cor.coef = T, cor.method = "pearson", xlab = "024 limma FDR",
                  ylab = "025 limma FDR", xlim = c(0, 1), ylim = c(0, 1), 
                  title = paste0("Correlation between 024 and 025 ", names(res_024)[i], "\nlimma differential gene expression results - FDR"))
  print(p2)
  dev.off()
  # correlation between logFC
  png(paste0("./Analysis_part2_GRCh38/024_025_limma_correlation_plots/all_shared_peaks/", 
             str_replace(names(res_024)[i], " ", "_"), "_all_gene_logfc_corr.png"))
  p3 <- ggscatter(data = data,
                  x = "logFC.x", y = "logFC.y", size = 0.01, add = "reg.line", conf.int = T, 
                  cor.coef = T, cor.method = "pearson", xlab = "024 limma logFC",
                  ylab = "025 limma logFC", #xlim = c(0, 1), ylim = c(0, 1), 
                  title = paste0("Correlation between 024 and 025 ", names(res_024)[i], "\nlimma differential gene expression results - logFC"))
  print(p3)
  dev.off()
  # correlation between -log10(p)
  png(paste0("./Analysis_part2_GRCh38/024_025_limma_correlation_plots/all_shared_peaks/", 
             str_replace(names(res_024)[i], " ", "_"), "_all_gene_neglogp_corr.png"))
  p4 <- ggscatter(data = data,
                  x = "neglogp.x", y = "neglogp.y", size = 0.01, add = "reg.line", conf.int = T, 
                  cor.coef = T, cor.method = "pearson", xlab = "024 limma -log10(p)",
                  ylab = "025 limma -log10(p)",
                  title = paste0("Correlation between 024 and 025 ", names(res_024)[i], "\nlimma differential gene expression results - -log10p"))
  print(p4)
  dev.off()
}

# peaks significant in either 024 or 025 ####
peaks <- vector("list", length(res_024))
res_024 <- res_024
res_025 <- res_025
sig_counts_table <- matrix(nrow = 3, ncol = 6, 
                           dimnames = list(c("024", "025", "union"), names(res_024)))

for (i in 1:length(res_024)) {
  intersect_peaks <- intersect(res_024[[i]]$gene, res_025[[i]]$gene)
  peaks[[i]] <- union(res_024[[i]]$gene[res_024[[i]]$gene %in% intersect_peaks & 
                                          res_024[[i]]$adj.P.Val < 0.05], 
                      res_025[[i]]$gene[res_025[[i]]$gene %in% intersect_peaks &
                                          res_025[[i]]$adj.P.Val < 0.05])
  sig_counts_table[1, i] <- length(res_024[[i]]$gene[res_024[[i]]$gene %in% intersect_peaks & 
                                                       res_024[[i]]$adj.P.Val < 0.05])
  sig_counts_table[2, i] <- length(res_025[[i]]$gene[res_025[[i]]$gene %in% intersect_peaks & 
                                                       res_025[[i]]$adj.P.Val < 0.05])
  sig_counts_table[3, i] <- length(peaks[[i]])
  
  res_024[[i]] <- res_024[[i]][res_024[[i]]$gene %in% peaks[[i]], ]
  res_025[[i]] <- res_025[[i]][res_025[[i]]$gene %in% peaks[[i]], ]
  
  res_024[[i]]$neglogp <- (-1) * log10(res_024[[i]]$P.Value)
  res_025[[i]]$neglogp <- (-1) * log10(res_025[[i]]$P.Value)
}
write.table(sig_counts_table, file = "./Analysis_part2_GRCh38/024_025_limma_correlation_plots/shared_significant_gene_counts.csv",
            sep = ",", quote = F, row.names = T, col.names = T)

for (i in 1:length(res_024)) {
  data <- res_024[[i]] %>% left_join(res_025[[i]], by = "gene")
  # correlation between p-values
  png(paste0("./Analysis_part2_GRCh38/024_025_limma_correlation_plots/significant_peaks_only/", 
             str_replace(names(res_024)[i], " ", "_"), "_sig_gene_pval_corr.png"))
  p1 <- ggscatter(data = data,
                  x = "P.Value.x", y = "P.Value.y", size = 0.01, add = "reg.line", conf.int = T, 
                  cor.coef = T, cor.method = "pearson", xlab = "024 limma p-value",
                  ylab = "025 limma p-value", xlim = c(0, 1), ylim = c(0, 1), 
                  title = paste0("Correlation between 024 and 025 ", names(res_024)[i], "\nlimma FDR<0.05 peaks - p value"))
  print(p1)
  dev.off()
  # correlation between FDR
  png(paste0("./Analysis_part2_GRCh38/024_025_limma_correlation_plots/significant_peaks_only/", 
             str_replace(names(res_024)[i], " ", "_"), "_sig_gene_fdr_corr.png"))
  p2 <- ggscatter(data = data,
                  x = "adj.P.Val.x", y = "adj.P.Val.y", size = 0.01, add = "reg.line", conf.int = T, 
                  cor.coef = T, cor.method = "pearson", xlab = "024 limma FDR",
                  ylab = "025 limma FDR", xlim = c(0, 1), ylim = c(0, 1), 
                  title = paste0("Correlation between 024 and 025 ", names(res_024)[i], "\nlimma FDR<0.05 peaks - FDR"))
  print(p2)
  dev.off()
  # correlation between logFC
  png(paste0("./Analysis_part2_GRCh38/024_025_limma_correlation_plots/significant_peaks_only/", 
             str_replace(names(res_024)[i], " ", "_"), "_sig_gene_logfc_corr.png"))
  p3 <- ggscatter(data = data,
                  x = "logFC.x", y = "logFC.y", size = 0.01, add = "reg.line", conf.int = T, 
                  cor.coef = T, cor.method = "pearson", xlab = "024 limma logFC",
                  ylab = "025 limma logFC", #xlim = c(0, 1), ylim = c(0, 1), 
                  title = paste0("Correlation between 024 and 025 ", names(res_024)[i], "\nlimma FDR<0.05 peaks - logFC"))
  print(p3)
  dev.off()
  # correlation between -log10(p)
  png(paste0("./Analysis_part2_GRCh38/024_025_limma_correlation_plots/significant_peaks_only/", 
             str_replace(names(res_024)[i], " ", "_"), "_sig_gene_neglogp_corr.png"))
  p4 <- ggscatter(data = data,
                  x = "neglogp.x", y = "neglogp.y", size = 0.01, add = "reg.line", conf.int = T, 
                  cor.coef = T, cor.method = "pearson", xlab = "024 limma -log10(p)",
                  ylab = "025 limma -log10(p)",
                  title = paste0("Correlation between 024 and 025 ", names(res_024)[i], "\nlimma FDR<0.05 peaks - -log10p"))
  print(p4)
  dev.off()
  
}
