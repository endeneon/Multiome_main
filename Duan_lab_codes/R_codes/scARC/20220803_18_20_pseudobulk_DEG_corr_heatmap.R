# Chuxuan Li 08/03/2022
# plot correlation between the log2FC of timextype specific DEGs in 
#18- and 20-line data

# init ####
library(ggplot2)
library(gplots)
library(Hmisc)
library(corrplot)
library(stringr)

# 1. all passed filter genes ####
res18 <- list.files(path = "/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm_9in18_res/", 
                      pattern = "0_all_DEGs.csv", full.names = T, recursive = T, include.dirs = F)
res20 <- list.files(path = "/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_pseudobulk/20line_limma/",
                    pattern = "0_all_DEGs.csv", full.names = T, recursive = T, include.dirs = F)

deg_lst_18 <- vector("list", length(res18))
deg_lst_20 <- vector("list", length(res20))

names <- rep_len("", 6*2)
for (i in 1:length(deg_lst_18)) {
  deg_lst_18[[i]] <- read_csv(res18[i], show_col_types = F)
  deg_lst_20[[i]] <- read_csv(res20[i], show_col_types = F)
  names[2*i-1] <- paste("18", 
                        str_extract(res18[i], "[A-Za-z]+_[A-Za-z]+_[1|6]v[0|1|6]"),
                        sep = "_")
  names[2*i] <- paste("20", 
                      str_extract(res20[i], "[A-Za-z]+_[A-Za-z]+_[1|6]v[0|1|6]"),
                      sep = "_")
}

for (i in 1:length(deg_lst_18)) {
  if (i == 1) {
    common_genes <- intersect(deg_lst_18[[i]]$genes, deg_lst_20[[i]]$genes)
  } else {
    temp <- intersect(deg_lst_18[[i]]$genes, deg_lst_20[[i]]$genes)
    common_genes <- intersect(common_genes, temp)
  }
}

logFC_df <- array(dim = c(length(common_genes), 6*2))
for (i in 1:length(deg_lst_18)) {
  fc <- deg_lst_18[[i]][deg_lst_18[[i]]$genes %in% common_genes, 1:2]
  fc <- fc[order(fc$genes), ]
  logFC_df[, 2*i-1] <- fc$logFC
  fc <- deg_lst_20[[i]][deg_lst_20[[i]]$genes %in% common_genes, 1:2]
  fc <- fc[order(fc$genes), ]
  logFC_df[, 2*i] <- fc$logFC
}

rownames(logFC_df) <- common_genes
colnames(logFC_df) <- str_replace(str_replace(str_replace_all(str_remove(names, "res"), "_", " "), "18", "18line"), "20", "20line")

# correlation matrix
corr_mat <- cor(logFC_df)
#corr_mat <- corr_mat[, str_detect(colnames(corr_mat), "20")]
#corr_mat <- corr_mat[str_detect(rownames(corr_mat), "18"), ]
corrplot(corr_mat)

palette = colorRampPalette(c("lightsteelblue1", "steelblue3", "black"))(50)
pheatmap::pheatmap(mat = corr_mat, color = palette, 
                   cluster_cols = T, cluster_rows = T, 
                   border_color = NA, angle_col = 315)

# 2. significant genes ####
