# Chuxuan Li 03/14/2023
# Regress against PRS within 3 time points separately

# init ####
library(Seurat)
library(readr)
library(readxl)
library(stringr)
library(limma)
library(edgeR)

load("covar_table_only_lines_w_PRS.RData")
load("sep_by_type_time_col_by_line_adj_mat_lst_w_only_lines_w_PRS.RData")

# Limma ####
createDGE <- function(count_matrix){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix))
  A <- rowSums(y$counts)
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}
filterByCpm <- function(y, cutoff, nsample) {
  cpm <- cpm(y)
  passfilter <- (rowSums(cpm >= cutoff) >= nsample)
  return(passfilter)
}
cnfV <- function(y, design) {
  y <- calcNormFactors(y)
  # voom (no logcpm)
  v <- voom(y, design, plot = T)
  
  return(v)  
}

# run limma by cell type and time point ####
covar_table_0hr <- covar_table_PRS[covar_table_PRS$time == "0hr", ]
covar_table_1hr <- covar_table_PRS[covar_table_PRS$time == "1hr", ]
covar_table_6hr <- covar_table_PRS[covar_table_PRS$time == "6hr", ]
covar_table_lst <- list(covar_table_0hr, covar_table_1hr, covar_table_6hr,
                        covar_table_0hr, covar_table_1hr, covar_table_6hr,
                        covar_table_0hr, covar_table_1hr, covar_table_6hr)
names(covar_table_lst) <- names(adj_mat_lst)

fit_all <- vector("list", length(adj_mat_lst))

for (i in 1:length(adj_mat_lst)) {
  type <- str_extract(names(adj_mat_lst)[i], "^[A-Za-z]+")
  time <- str_extract(names(adj_mat_lst)[i], "[0|1|6]hr")
  if (type == "GABA") {
    design <- model.matrix(~PRS + batch + age + sex + GABA_fraction, 
                                data = covar_table_lst[[i]])
  } else if (type == "nmglut") {
    design <- model.matrix(~PRS + batch + age + sex + nmglut_fraction, 
                           data = covar_table_lst[[i]])
  } else {
    design <- model.matrix(~PRS + batch + age + sex + npglut_fraction, 
                           data = covar_table_lst[[i]])
  }
  y <- createDGE(adj_mat_lst[[i]]) # plug in adjusted or raw matrix
  ind.keep <- filterByCpm(y, 1, 24)
  print(sum(ind.keep)) #17775\
  v <- cnfV(y[ind.keep, ], design)
  fit <- lmFit(v, design)
  fit_all[[i]] <- eBayes(fit)
}


setwd("./sep_by_celltype_time_DE_results/")
types <- names(adj_mat_lst)
for (i in 1:length(fit_all)){
  #res <- topTable(fit_all[[i]], coef = "PRS", p.value = 0.05, sort.by = "P", number = Inf)
  res <- topTable(fit_all[[i]], coef = "PRS", p.value = Inf, sort.by = "P", number = Inf)
  res$gene <- rownames(res)
  
  filename <- paste0("./upregulated_significant/", types[i], "_PRS_upregulated_DEGs.csv")
  print(filename)
  write.table(res[res$logFC > 0 & res$P.Value < 0.05, ], file = filename, quote = F, sep = ",", 
              row.names = F, col.names = T)
  filename <- paste0("./downregulated_significant/", types[i], "_PRS_downregulated_DEGs.csv")
  print(filename)
  write.table(res[res$logFC < 0 & res$P.Value < 0.05, ], file = filename, quote = F, sep = ",", 
              row.names = F, col.names = T)  
  
  res <- topTable(fit_all[[i]], coef = "PRS", p.value = Inf, sort.by = "P", number = Inf)
  res$gene <- rownames(res)
  
  filename <- paste0("./unfiltered_by_padj/", types[i], "_PRS_upregulated_DEGs.csv")
  print(filename)
  write.table(res[res$logFC > 0, ], file = filename, 
              quote = F, sep = ",", row.names = F, col.names = T)
  
  filename <- paste0("./unfiltered_by_padj/", types[i], "_PRS_downregulated_DEGs.csv")
  print(filename)
  write.table(res[res$logFC < 0, ], file = filename, 
              quote = F, sep = ",", row.names = F, col.names = T)
  
  filename <- paste0("./unfiltered_by_padj/", types[i], "_PRS_all_DEGs.csv")
  print(filename)
  write.table(res, file = filename, 
              quote = F, sep = ",", row.names = F, col.names = T)
}

# summary table
cnames <- c("beta > 0", "beta < 0", "beta > 0 (p < 0.05)", "beta < 0 (p < 0.05)",
            "total passed filter")
deg_counts <- array(dim = c(9, 5), dimnames = list(types, cnames)) # 3 cell types are rows, 1/6 up/down are cols
for (i in 1:length(fit_all)){                                                               
  res <- topTable(fit_all[[i]], coef = "PRS", p.value = Inf, sort.by = "P", number = Inf)
  deg_counts[i, 1] <- sum(res$logFC > 0 & res$adj.P.Val < 0.05)
  deg_counts[i, 2] <- sum(res$logFC < 0 & res$adj.P.Val < 0.05)
  deg_counts[i, 3] <- sum(res$logFC > 0 & res$P.Value < 0.05)
  deg_counts[i, 4] <- sum(res$logFC < 0 & res$P.Value < 0.05)
  deg_counts[i, 5] <- nrow(res)
}
write.table(deg_counts, file = "./PRS_assoc_gene_count_summary.csv", quote = F, sep = ",", 
            row.names = T, col.names = T)

# count number of genes overlapping between time points ####
overlap_counts <- matrix(nrow = length(fit_all), ncol = 10, 
                         dimnames = list(names(adj_mat_lst), c("beta > 0 overlap with 0hr",
                                                               "beta > 0 overlap with 1hr",
                                                               "beta > 0 overlap with 6hr",
                                                               "beta > 0 3-way overlap",
                                                               "total beta > 0",
                                                               "beta < 0 overlap with 0hr",
                                                               "beta < 0 overlap with 1hr",
                                                               "beta < 0 overlap with 6hr",
                                                               "beta < 0 3-way overlap",
                                                               "total beta < 0")))
posb_lst <- vector("list", length(fit_all))
negb_lst <- vector("list", length(fit_all))
for (i in 1:length(fit_all)){                                                               
  res <- topTable(fit_all[[i]], coef = "PRS", p.value = Inf, sort.by = "P", number = Inf)
  posb_lst[[i]] <- res[res$logFC > 0 & res$P.Value < 0.05, ]
  negb_lst[[i]] <- res[res$logFC < 0 & res$P.Value < 0.05, ]
  overlap_counts[i, 5] <- sum(res$logFC > 0 & res$P.Value < 0.05)
  overlap_counts[i, 10] <- sum(res$logFC < 0 & res$P.Value < 0.05)
}
names(posb_lst) <- names(adj_mat_lst)
names(negb_lst) <- names(adj_mat_lst)
overlap_counts[c(1,4,7), c(1, 6)] <- "/"
overlap_counts[c(2,5,8), c(2, 7)] <- "/"
overlap_counts[c(3,6,9), c(3, 8)] <- "/"
for (i in c(1, 4, 7)) {
  # positive beta
  overlap_0n1 <- intersect(posb_lst[[i]]$genes, posb_lst[[i + 1]]$genes)
  write.table(overlap_0n1, 
              file = paste0("./overlapped_genes_pval_significant/", 
                            names(posb_lst)[i], "_overlap_w_1hr_genes_positive_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_0n6 <- intersect(posb_lst[[i]]$genes, posb_lst[[i + 2]]$genes)
  write.table(overlap_0n6, 
              file = paste0("overlapped_genes_pval_significant/", 
                            names(posb_lst)[i], "_overlap_w_6hr_genes_positive_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_counts[i, 2] <- length(overlap_0n1)
  overlap_counts[i + 1, 1] <- length(overlap_0n1)
  overlap_counts[i, 3] <- length(overlap_0n6)
  overlap_counts[i + 2, 1] <- length(overlap_0n6)
  # negative beta
  overlap_0n1 <- intersect(negb_lst[[i]]$genes, negb_lst[[i + 1]]$genes)
  write.table(overlap_0n1, 
              file = paste0("overlapped_genes_pval_significant/", 
                            names(negb_lst)[i], "_overlap_w_1hr_genes_negative_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_0n6 <- intersect(negb_lst[[i]]$genes, negb_lst[[i + 2]]$genes)
  write.table(overlap_0n6, 
              file = paste0("overlapped_genes_pval_significant/", 
                            names(negb_lst)[i], "_overlap_w_6hr_genes_negative_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_counts[i, 7] <- length(overlap_0n1)
  overlap_counts[i + 1, 6] <- length(overlap_0n1)
  overlap_counts[i, 8] <- length(overlap_0n6)
  overlap_counts[i + 2, 6] <- length(overlap_0n6)
  
  overlap_3w <- intersect(intersect(posb_lst[[i]]$genes, posb_lst[[i + 1]]$genes), posb_lst[[i + 2]]$genes)
  write.table(overlap_3w, file = paste0("overlapped_genes_pval_significant/", 
                                        names(posb_lst)[i], "_3way_overlap_positive_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_counts[i:(i+2), 4] <- length(overlap_3w)
  only_in_0hr <- posb_lst[[i]]$genes[!posb_lst[[i]]$genes %in% overlap_3w]
  write.table(overlap_3w, file = paste0("unique_to_time_point_genes_pval_significant/", 
                                        names(posb_lst)[i], "_only_in_0hr_positive_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  only_in_1hr <- posb_lst[[i + 1]]$genes[!posb_lst[[i + 1]]$genes %in% overlap_3w]
  write.table(overlap_3w, file = paste0("unique_to_time_point_genes_pval_significant/", 
                                        names(posb_lst)[i + 1], "_only_in_1hr_positive_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  only_in_6hr <- posb_lst[[i + 2]]$genes[!posb_lst[[i + 2]]$genes %in% overlap_3w]
  write.table(overlap_3w, file = paste0("unique_to_time_point_genes_pval_significant/", 
                                        names(posb_lst)[i + 2], "_only_in_6hr_positive_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  
  overlap_3w <- intersect(intersect(negb_lst[[i]]$genes, negb_lst[[i + 1]]$genes), negb_lst[[i + 2]]$genes)
  write.table(overlap_3w, file = paste0("overlapped_genes_pval_significant/", 
                                        names(negb_lst)[i], "_3way_overlap_negative_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_counts[i:(i+2), 9] <- length(overlap_3w)
  only_in_0hr <- negb_lst[[i]]$genes[!negb_lst[[i]]$genes %in% overlap_3w]
  write.table(overlap_3w, file = paste0("unique_to_time_point_genes_pval_significant/", 
                                        names(negb_lst)[i], "_only_in_0hr_negative_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  only_in_1hr <- negb_lst[[i + 1]]$genes[!negb_lst[[i + 1]]$genes %in% overlap_3w]
  write.table(overlap_3w, file = paste0("unique_to_time_point_genes_pval_significant/", 
                                        names(negb_lst)[i + 1], "_only_in_1hr_negative_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  only_in_6hr <- negb_lst[[i + 2]]$genes[!negb_lst[[i + 2]]$genes %in% overlap_3w]
  write.table(overlap_3w, file = paste0("unique_to_time_point_genes_pval_significant/", 
                                        names(negb_lst)[i + 2], "_only_in_6hr_negative_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
}
for (i in c(2, 5, 8)) {
  # positive beta
  overlap_1n6 <- intersect(posb_lst[[i]]$genes, posb_lst[[i + 1]]$genes)
  write.table(overlap_1n6, 
              file = paste0("overlapped_genes_pval_significant/", 
                            names(posb_lst)[i], "_overlap_w_6hr_genes_positive_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_counts[i, 3] <- length(overlap_1n6)
  overlap_counts[i + 1, 2] <- length(overlap_1n6)
  # negative beta
  overlap_1n6 <- intersect(negb_lst[[i]]$genes, negb_lst[[i + 1]]$genes)
  write.table(overlap_1n6, 
              file = paste0("overlapped_genes_pval_significant/", 
                            names(negb_lst)[i], "_overlap_w_6hr_genes_negative_beta.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_counts[i, 8] <- length(overlap_1n6)
  overlap_counts[i + 1, 7] <- length(overlap_1n6)
}
write.table(overlap_counts, file = "overlapping_gene_count_summary.csv", 
            quote = F, sep = ",", row.names = T, col.names = T)

