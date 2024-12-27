# Chuxuan li 03/18/2022
# For 18-line data, plot the log2FC of genes from disease-related gene lists using the 
#pseudobulk_DE_by_limma_batch_in_design_and_voom.R results

# init ####
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(readxl)
library(rowr)

SZ <- read_excel("./SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                 sheet = "prioritized SZ single genes")
colnames(SZ) <- "gene"
ASD <- read_excel("SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                  sheet = "ASD102_Cell_exomSeq")
colnames(ASD) <- "gene"
BP <- read_excel("SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                 sheet = "BP_Nat Gent 2021_MAGMA genes")
colnames(BP) <- "gene"
MDD <- read_excel("SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                 sheet = "MDD_NN2019_MAGMA", skip = 1)
MDD <- as.data.frame(MD$`Gene Name`)
colnames(MDD) <- "gene"
PTSD <- read_excel("SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                   sheet = "PTSD_NG2021", col_names = FALSE, 
                   skip = 2)
PTSD <- as.data.frame(PTSD)
colnames(PTSD) <- "gene"

SZ_SCHEMA <- read_excel("SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                 sheet = "SZ_SCHEMA")
SZ_SCHEMA <- as.data.frame(SZ_SCHEMA$...17)
colnames(SZ_SCHEMA) <- "gene"

title_list <- c("prioritized SZ single genes",
                "ASD102_Cell_exomSeq",
                "BP_Nat Gent 2021_MAGMA genes",
                "MDD_NN2019_MAGMA",
                "PTSD_NG2021",
                "SZ_SCHEMA")
plotlist <- vector(mode = "list", length = length(title_list))

res_0v1_list <- vector(mode = "list", length = 3L)
res_0v6_list <- vector(mode = "list", length = 3L)
path_1 <- sort(list.files(path = "./pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm_9in18_res/",
                     pattern = "1v0_all", full.names = T, recursive = T))
path_6 <- sort(list.files(path = "./pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm_9in18_res/",
                     pattern = "6v0_all", full.names = T, recursive = T))
for (i in 1:length(res_0v1_list)){
  res_0v1_list[[i]] <- read_csv(path_1[i])
  res_0v6_list[[i]] <- read_csv(path_6[i])
}

ttlabel <- rep(str_remove(str_extract(path_1, "[A-Za-z]+_res_1v0"), "_res"), each = 2)
ttlabel[c(2, 4, 6)] <- str_replace(ttlabel[c(2, 4, 6)], "1v0", "6v0")
ttlabel

# plot ####
k = 0
for (file in list(SZ, ASD, BP, MDD, PTSD, SZ_SCHEMA)){
  k = k + 1
  gene_list <- res_0v1_list[[1]][res_0v1_list[[1]]$genes %in% file$gene, ]$genes
  # select rows from full list based on genelist
  for (i in 1:length(res_0v1_list)){
    if (i == 1){
      df_to_plot <- res_0v1_list[[i]][res_0v1_list[[i]]$genes %in% gene_list, 1:2]
      df_0v6 <- res_0v6_list[[i]][res_0v6_list[[i]]$genes %in% gene_list, 1:2]
      df_to_plot <- merge(x = df_to_plot, y = df_0v6, by = "genes")
    } else {
      df_to_append_0v1 <- res_0v1_list[[i]][res_0v1_list[[i]]$genes %in% gene_list, 1:2]
      df_to_append_0v6 <- res_0v6_list[[i]][res_0v6_list[[i]]$genes %in% gene_list, 1:2]
      df_to_plot <- merge(x = df_to_plot, y = df_to_append_0v1, by = "genes")
      df_to_plot <- merge(x = df_to_plot, y = df_to_append_0v6, by = "genes")
    }
  }
  colnames(df_to_plot) <- c("genes", ttlabel) # genes + cell type x time
  rownames(df_to_plot) <- df_to_plot$genes
  df_to_plot$genes <- NULL # remove the gene column
  df_clean <- df_to_plot[apply(df_to_plot, 1, var) != 0, ] # remove rows with zero variance (for heatmap)
  gene_list <- rownames(df_clean)
  print(head(gene_list))
  
  # get q values for each gene, store in a df
  for (i in 1:length(res_0v1_list)){
    if (i == 1){
      df_qval_sig <- res_0v1_list[[i]][res_0v1_list[[i]]$genes %in% gene_list, c(1, 6)]
      df_qval_num <- df_qval_sig
      df_qval_num[df_qval_num$adj.P.Val >= 0.05, 2] <- 0 # all nonsig ones assigned 0
      df_0v6 <- res_0v6_list[[i]][res_0v6_list[[i]]$genes %in% gene_list, c(1, 6)]
      df_0v6_num <- df_0v6
      df_0v6_num[df_0v6$adj.P.Val >= 0.05, 2] <- 0
      df_qval_num <- merge(x = df_qval_num, y = df_0v6_num, by = "genes")
    } else {
      df_to_append_0v1 <- res_0v1_list[[i]][res_0v1_list[[i]]$genes %in% gene_list, c(1, 6)]
      df_to_append_0v6 <- res_0v6_list[[i]][res_0v6_list[[i]]$genes %in% gene_list, c(1, 6)]
      df_to_append_0v1[df_to_append_0v1$adj.P.Val >= 0.05, 2] <- 0
      df_to_append_0v6[df_to_append_0v6$adj.P.Val >= 0.05, 2] <- 0
      
      df_qval_num <- merge(x = df_qval_num, y = df_to_append_0v1, by = "genes")
      df_qval_num <- merge(x = df_qval_num, y = df_to_append_0v6, by = "genes")
    }
  }
  colnames(df_qval_num) <- c("genes", ttlabel)
  rownames(df_qval_num) <- df_qval_num$genes
  df_qval_num$genes <- NULL
  df_clean <- df_clean[rowSums(df_qval_num) > 0, ] # use this to remove genes with all nonsig qval from the fc df
  df_clean <- df_clean[order(df_clean$GABA_1v0, decreasing = T), ] # order genes from most positive fc to most negative
  gene_list <- rownames(df_clean) # get the genes from cleaned q value df
  label <- rownames(df_clean)
  df_qval_num <- df_qval_num[gene_list, ]
  df_qval <- df_qval_num
  df_qval[df_qval_num != 0] <- "*" # significant ones has a *
  df_qval[df_qval_num == 0] <- "" # nonsignificant qvalues has no labeling
  
  print(head(rownames(df_clean)))
  print(head(rownames(df_qval)))
  print(sum(is.na(df_clean)))

  p <- pheatmap(df_clean, cluster_rows = F, cluster_cols = F, # don't reorder rows or cols or
                scale = "none", # don't calculate Z score
                display_numbers = df_qval, # show * for significant q values
                clustering_distance_rows = "euclidean",
                treeheight_row = 0, labels_row = label, angle_col = 45,
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(45), # positive = red, negative = blue
                breaks = seq(-2, 2, 0.1), # limit the range of values (any value outside of (-2,2) will cohere to (-2,2))
                border_color = NA,
                fontsize = 8, fontsize_row = 6, fontsize_number = 5,
                cellwidth = 28, cellheight = 10,
                main = title_list[k])
  plotlist[[k]] <- p
}

print(plotlist[[1]])

