# Chuxuan li 03/18/2022
# For 5-line data, plot the log2FC of genes from certain disease-related gene lists using the 
#pseudobulk results

# init ####
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(rowr)

SZ <- read_excel("SZ_BP_ASD_MD_PTSD_genes.xlsx", 
                 sheet = "prioritized SZ single genes")
colnames(SZ) <- "gene"
ASD <- read_excel("SZ_BP_ASD_MD_PTSD_genes.xlsx", 
                  sheet = "ASD102_Cell_exomSeq")
colnames(ASD) <- "gene"
BP <- read_excel("SZ_BP_ASD_MD_PTSD_genes.xlsx", 
                 sheet = "BP_Nat Gent 2021_MAGMA genes")
colnames(BP) <- "gene"
MD <- read_excel("SZ_BP_ASD_MD_PTSD_genes.xlsx", 
                 sheet = "MDD_NN2019_MAGMA", skip = 1)
MD <- as.data.frame(MD$`Gene Name`)
colnames(MD) <- "gene"
PTSD <- read_excel("SZ_BP_ASD_MD_PTSD_genes.xlsx", 
                   sheet = "PTSD_NG2021", col_names = FALSE, 
                   skip = 2)
PTSD <- as.data.frame(PTSD)
colnames(PTSD) <- "gene"

plotlist <- vector(mode = "list", length = 4L)
title_list <- c("prioritized SZ single genes",
                "ASD102_Cell_exomSeq",
                "BP_Nat Gent 2021_MAGMA genes",
                "MDD_NN2019_MAGMA",
                "PTSD_NG2021")

res_0v1_list <- vector(mode = "list", length = 4L)
res_0v6_list <- vector(mode = "list", length = 4L)
path_1 <- list.files(path = "/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_pseudobulk/5line_DE_results/filtered_by_basemean",
           pattern = "*1v0*", full.names = T)
path_6 <- list.files(path = "/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_pseudobulk/5line_DE_results/filtered_by_basemean",
                           pattern = "*6v0*", full.names = T)
for (i in 1:length(res_0v1_list)){
  res_0v1_list[[i]] <- read.csv(path_1[i], header = F, skip = 1)
  colnames(res_0v1_list[[i]]) <- c("gene", "baseMean",
                                   "log2FoldChange", "lfcSE",	"stat",
                                   "pvalue",	"padj",	"q_value")
  res_0v6_list[[i]] <- read.csv(path_6[i], header = F, skip = 1)
  colnames(res_0v6_list[[i]]) <- c("gene", "baseMean",
                                   "log2FoldChange", "lfcSE",	"stat",
                                   "pvalue",	"padj",	"q_value")
}

# plot ####
k = 0
for (file in list(SZ, ASD, BP, MD, PTSD)){
  k = k + 1
  gene_list <- res_0v1_list[[1]][res_0v1_list[[1]]$gene %in% file$gene, ]$gene
  print(head(gene_list))
  
  for (i in 1:length(res_0v1_list)){
    if (i == 1){
      df_to_plot <- res_0v1_list[[i]][res_0v1_list[[i]]$gene %in% gene_list, c(1, 3)]
      df_0v6 <- res_0v6_list[[i]][res_0v6_list[[i]]$gene %in% gene_list, c(1, 3)]
      df_to_plot <- merge(x = df_to_plot, y = df_0v6, by = "gene")
    } else {
      df_to_append_0v1 <- res_0v1_list[[i]][res_0v1_list[[i]]$gene %in% gene_list, c(1, 3)]
      df_to_append_0v6 <- res_0v6_list[[i]][res_0v6_list[[i]]$gene %in% gene_list, c(1, 3)]
      df_to_plot <- merge(x = df_to_plot, y = df_to_append_0v1, by = "gene")
      df_to_plot <- merge(x = df_to_plot, y = df_to_append_0v6, by = "gene")
    }
  }
  colnames(df_to_plot) <- c("gene", "GABA_1v0", "GABA_6v0", "nmglut_1v0", "nmglut_6v0",
                            "NPC_1v0", "NPC_6v0",
                            "npglut_1v0", "npglut_6v0")
  rownames(df_to_plot) <- df_to_plot$gene
  df_to_plot <- df_to_plot[, 2:9]
  df_clean <- df_to_plot[apply(df_to_plot, 1, var) != 0, ]
  gene_list <- rownames(df_clean)
  print(head(gene_list))

  for (i in 1:length(res_0v1_list)){
    if (i == 1){
      df_qval_sig <- res_0v1_list[[i]][res_0v1_list[[i]]$gene %in% gene_list, c(1, 8)]
      df_qval_num <- df_qval_sig
      df_qval_num[df_qval_num$q_value >= 0.05, 2] <- 0
      df_0v6 <- res_0v6_list[[i]][res_0v6_list[[i]]$gene %in% gene_list, c(1, 8)]
      df_0v6_num <- df_0v6
      df_0v6_num[df_0v6$q_value >= 0.05, 2] <- 0
      df_qval_num <- merge(x = df_qval_num, y = df_0v6_num, by = "gene")
    } else {
      df_to_append_0v1 <- res_0v1_list[[i]][res_0v1_list[[i]]$gene %in% gene_list, c(1, 8)]
      df_to_append_0v6 <- res_0v6_list[[i]][res_0v6_list[[i]]$gene %in% gene_list, c(1, 8)]
      df_to_append_0v1[df_to_append_0v1$q_value >= 0.05, 2] <- 0
      df_to_append_0v6[df_to_append_0v6$q_value >= 0.05, 2] <- 0
      
      df_qval_num <- merge(x = df_qval_num, y = df_to_append_0v1, by = "gene")
      df_qval_num <- merge(x = df_qval_num, y = df_to_append_0v6, by = "gene")
    }
  }
  colnames(df_qval_num) <- c("gene", "GABA_1v0", "GABA_6v0", "nmglut_1v0", "nmglut_6v0",
                         "NPC_1v0", "NPC_6v0",
                         "npglut_1v0", "npglut_6v0")
  df_qval_num$gene <- NULL
  
  df_clean <- df_clean[rowSums(df_qval_num) > 0, ]
  gene_list <- rownames(df_clean)
  print(head(gene_list))

  for (i in 1:length(res_0v1_list)){
    if (i == 1){
      df_qval <- res_0v1_list[[i]][res_0v1_list[[i]]$gene %in% gene_list, c(1, 8)]
      df_qval[df_qval$q_value < 0.05, 2] <- "*"
      df_qval[df_qval$q_value >= 0.05, 2] <- ""
      df_0v6 <- res_0v6_list[[i]][res_0v6_list[[i]]$gene %in% gene_list, c(1, 8)]
      df_0v6[df_0v6$q_value < 0.05, 2] <- "*"
      df_0v6[df_0v6$q_value >= 0.05, 2] <- ""
      df_qval <- merge(x = df_qval, y = df_0v6, by = "gene")
    } else {
      df_to_append_0v1 <- res_0v1_list[[i]][res_0v1_list[[i]]$gene %in% gene_list, c(1, 8)]
      df_to_append_0v6 <- res_0v6_list[[i]][res_0v6_list[[i]]$gene %in% gene_list, c(1, 8)]
      df_to_append_0v1[df_to_append_0v1$q_value < 0.05, 2] <- "*"
      df_to_append_0v1[df_to_append_0v1$q_value >= 0.05, 2] <- ""
      df_to_append_0v6[df_to_append_0v6$q_value < 0.05, 2] <- "*"
      df_to_append_0v6[df_to_append_0v6$q_value >= 0.05, 2] <- ""
      df_qval <- merge(x = df_qval, y = df_to_append_0v1, by = "gene")
      df_qval <- merge(x = df_qval, y = df_to_append_0v6, by = "gene")
    }
  }
  colnames(df_qval) <- c("gene", "GABA_1v0", "GABA_6v0", "nmglut_1v0", "nmglut_6v0",
                            "NPC_1v0", "NPC_6v0",
                            "npglut_1v0", "npglut_6v0")
  rownames(df_qval) <- df_qval$gene
  df_qval$gene <- NULL
  print(head(rownames(df_clean)))
  print(head(rownames(df_qval)))
  print(sum(is.na(df_clean)))
  df_clean <- df_clean[order(df_clean$GABA_1v0, decreasing = T), ]
  label <- rownames(df_clean)
  p <- pheatmap(df_clean, cluster_rows = F, cluster_cols = F, scale = "none",
           display_numbers = df_qval,
           clustering_distance_rows = "euclidean",
           treeheight_row = 0, labels_row = label, angle_col = 45,
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(45),
           breaks = seq(-2, 2, 0.1),
           border_color = NA,
           fontsize = 8, fontsize_row = 6,
           main = title_list[k])
  plotlist[[k]] <- p
}


print(plotlist[[1]])

