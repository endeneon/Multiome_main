# Chuxuan Li 03/14/2023
# Regression with response to clozapine as the main variable of interest

# init ####
library(Seurat)
library(readr)
library(readxl)
library(stringr)
library(limma)
library(edgeR)

load("covar_table_only_lines_w_clozapine_response.RData")
load("sep_by_type_time_col_by_line_adj_mat_lst_w_only_lines_w_clozapine.RData")

# auxiliary functions ####
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
covar_table_0hr <- covar_table_clozapine[covar_table_clozapine$time == "0hr", ]
covar_table_1hr <- covar_table_clozapine[covar_table_clozapine$time == "1hr", ]
covar_table_6hr <- covar_table_clozapine[covar_table_clozapine$time == "6hr", ]
covar_table_lst <- list(covar_table_0hr, covar_table_1hr, covar_table_6hr,
                        covar_table_0hr, covar_table_1hr, covar_table_6hr,
                        covar_table_0hr, covar_table_1hr, covar_table_6hr)
names(covar_table_lst) <- names(adj_mat_lst)

fit_all <- vector("list", length(adj_mat_lst))

for (i in 1:length(adj_mat_lst)) {
  type <- str_extract(names(adj_mat_lst)[i], "^[A-Za-z]+")
  time <- str_extract(names(adj_mat_lst)[i], "[0|1|6]hr")
  if (type == "GABA") {
    design <- model.matrix(~0 + clozapine + batch + age + sex + GABA_fraction, 
                           data = covar_table_lst[[i]])
  } else if (type == "nmglut") {
    design <- model.matrix(~0 + clozapine + batch + age + sex + nmglut_fraction, 
                           data = covar_table_lst[[i]])
  } else {
    design <- model.matrix(~0 + clozapine + batch + age + sex + npglut_fraction, 
                           data = covar_table_lst[[i]])
  }
  y <- createDGE(adj_mat_lst[[i]]) # plug in adjusted or raw matrix
  ind.keep <- filterByCpm(y, 1, 13)
  print(sum(ind.keep)) #17775\
  v <- cnfV(y[ind.keep, ], design)
  fit <- lmFit(v, design)
  contr.matrix <- makeContrasts(
      cloyesvsno = clozapineyes-clozapineno, levels = colnames(design))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit_all[[i]] <- eBayes(fit, trend = T)
  print(summary(decideTests(fit)))
}

# output results ####
setwd("./sep_by_celltype_time_DE_results/")
types <- names(adj_mat_lst)
for (i in 1:length(fit_all)){
  #res <- topTable(fit_all[[i]], coef = "cloyesvsno", p.value = 0.05, sort.by = "P", number = Inf)
  res <- topTable(fit_all[[i]], coef = "cloyesvsno", p.value = Inf, sort.by = "P", number = Inf)
  res$gene <- rownames(res)
  
  filename <- paste0("./upregulated_significant/", types[i], "_clozapine_upregulated_DEGs.csv")
  print(filename)
  write.table(res[res$logFC > 0 & res$P.Value < 0.05, ], file = filename, quote = F, sep = ",", 
              row.names = F, col.names = T)
  filename <- paste0("./downregulated_significant/", types[i], "_clozapine_downregulated_DEGs.csv")
  print(filename)
  write.table(res[res$logFC < 0 & res$P.Value < 0.05, ], file = filename, quote = F, sep = ",", 
              row.names = F, col.names = T)  
  
  res <- topTable(fit_all[[i]], coef = "cloyesvsno", p.value = Inf, sort.by = "P", number = Inf)
  res$gene <- rownames(res)
  
  filename <- paste0("./unfiltered_by_padj/", types[i], "_clozapine_upregulated_DEGs.csv")
  print(filename)
  write.table(res[res$logFC > 0, ], file = filename, 
              quote = F, sep = ",", row.names = F, col.names = T)
  
  filename <- paste0("./unfiltered_by_padj/", types[i], "_clozapine_downregulated_DEGs.csv")
  print(filename)
  write.table(res[res$logFC < 0, ], file = filename, 
              quote = F, sep = ",", row.names = F, col.names = T)
  
  filename <- paste0("./unfiltered_by_padj/", types[i], "_clozapine_all_DEGs.csv")
  print(filename)
  write.table(res, file = filename, 
              quote = F, sep = ",", row.names = F, col.names = T)
}

# summary table
cnames <- c("upregulated in clozapine responsive", "downregulated in clozapine responsive",
            "upregulated (p < 0.05)", "downregulated (p < 0.05)",
            "total passed filter")
deg_counts <- array(dim = c(9, 5), dimnames = list(types, cnames)) # 3 cell types are rows, 1/6 up/down are cols
for (i in 1:length(fit_all)){                                                               
  res <- topTable(fit_all[[i]], coef = "cloyesvsno", p.value = Inf, sort.by = "P", number = Inf)
  deg_counts[i, 1] <- sum(res$logFC > 0 & res$adj.P.Val < 0.05)
  deg_counts[i, 2] <- sum(res$logFC < 0 & res$adj.P.Val < 0.05)
  deg_counts[i, 3] <- sum(res$logFC > 0 & res$P.Value < 0.05)
  deg_counts[i, 4] <- sum(res$logFC < 0 & res$P.Value < 0.05)
  deg_counts[i, 5] <- nrow(res)
}
write.table(deg_counts, file = "./clozapine_assoc_gene_count_summary.csv", quote = F, sep = ",", 
            row.names = T, col.names = T)

# count number of genes overlapping between time points ####
overlap_counts <- matrix(nrow = length(fit_all), ncol = 10, 
                         dimnames = list(names(adj_mat_lst), c("upregulated overlap with 0hr",
                                                           "upregulated overlap with 1hr",
                                                           "upregulated overlap with 6hr",
                                                           "upregulated 3-way overlap",
                                                           "total upregulated",
                                                           "downregulated overlap with 0hr",
                                                           "downregulated overlap with 1hr",
                                                           "downregulated overlap with 6hr",
                                                           "downregulated 3-way overlap",
                                                           "total downregulated")))
posb_lst <- vector("list", length(fit_all))
negb_lst <- vector("list", length(fit_all))
for (i in 1:length(fit_all)){                                                               
  res <- topTable(fit_all[[i]], coef = "cloyesvsno", p.value = Inf, sort.by = "P", number = Inf)
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
                            names(posb_lst)[i], "_overlap_w_1hr_genes_upregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_0n6 <- intersect(posb_lst[[i]]$genes, posb_lst[[i + 2]]$genes)
  write.table(overlap_0n6, 
              file = paste0("overlapped_genes_pval_significant/", 
                            names(posb_lst)[i], "_overlap_w_6hr_genes_upregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_counts[i, 2] <- length(overlap_0n1)
  overlap_counts[i + 1, 1] <- length(overlap_0n1)
  overlap_counts[i, 3] <- length(overlap_0n6)
  overlap_counts[i + 2, 1] <- length(overlap_0n6)
  # negative beta
  overlap_0n1 <- intersect(negb_lst[[i]]$genes, negb_lst[[i + 1]]$genes)
  write.table(overlap_0n1, 
              file = paste0("overlapped_genes_pval_significant/", 
                            names(negb_lst)[i], "_overlap_w_1hr_genes_downregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_0n6 <- intersect(negb_lst[[i]]$genes, negb_lst[[i + 2]]$genes)
  write.table(overlap_0n6, 
              file = paste0("overlapped_genes_pval_significant/", 
                            names(negb_lst)[i], "_overlap_w_6hr_genes_downregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_counts[i, 7] <- length(overlap_0n1)
  overlap_counts[i + 1, 6] <- length(overlap_0n1)
  overlap_counts[i, 8] <- length(overlap_0n6)
  overlap_counts[i + 2, 6] <- length(overlap_0n6)
  
  overlap_3w <- intersect(intersect(posb_lst[[i]]$genes, posb_lst[[i + 1]]$genes), posb_lst[[i + 2]]$genes)
  write.table(overlap_3w, file = paste0("overlapped_genes_pval_significant/", 
                                        names(posb_lst)[i], "_3way_overlap_upregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_counts[i:(i+2), 4] <- length(overlap_3w)
  only_in_0hr <- posb_lst[[i]]$genes[!posb_lst[[i]]$genes %in% overlap_3w]
  write.table(overlap_3w, file = paste0("unique_to_time_point_genes_pval_significant/", 
                                        names(posb_lst)[i], "_only_in_0hr_upregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  only_in_1hr <- posb_lst[[i + 1]]$genes[!posb_lst[[i + 1]]$genes %in% overlap_3w]
  write.table(overlap_3w, file = paste0("unique_to_time_point_genes_pval_significant/", 
                                        names(posb_lst)[i + 1], "_only_in_1hr_upregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  only_in_6hr <- posb_lst[[i + 2]]$genes[!posb_lst[[i + 2]]$genes %in% overlap_3w]
  write.table(overlap_3w, file = paste0("unique_to_time_point_genes_pval_significant/", 
                                        names(posb_lst)[i + 2], "_only_in_6hr_upregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  
  overlap_3w <- intersect(intersect(negb_lst[[i]]$genes, negb_lst[[i + 1]]$genes), negb_lst[[i + 2]]$genes)
  write.table(overlap_3w, file = paste0("overlapped_genes_pval_significant/", 
                                        names(negb_lst)[i], "_3way_overlap_downregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_counts[i:(i+2), 9] <- length(overlap_3w)
  only_in_0hr <- negb_lst[[i]]$genes[!negb_lst[[i]]$genes %in% overlap_3w]
  write.table(overlap_3w, file = paste0("unique_to_time_point_genes_pval_significant/", 
                                        names(negb_lst)[i], "_only_in_0hr_downregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  only_in_1hr <- negb_lst[[i + 1]]$genes[!negb_lst[[i + 1]]$genes %in% overlap_3w]
  write.table(overlap_3w, file = paste0("unique_to_time_point_genes_pval_significant/", 
                                        names(negb_lst)[i + 1], "_only_in_1hr_downregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  only_in_6hr <- negb_lst[[i + 2]]$genes[!negb_lst[[i + 2]]$genes %in% overlap_3w]
  write.table(overlap_3w, file = paste0("unique_to_time_point_genes_pval_significant/", 
                                        names(negb_lst)[i + 2], "_only_in_6hr_downregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
}
for (i in c(2, 5, 8)) {
  # positive beta
  overlap_1n6 <- intersect(posb_lst[[i]]$genes, posb_lst[[i + 1]]$genes)
  write.table(overlap_1n6, 
              file = paste0("overlapped_genes_pval_significant/", 
                            names(posb_lst)[i], "_overlap_w_6hr_genes_upregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_counts[i, 3] <- length(overlap_1n6)
  overlap_counts[i + 1, 2] <- length(overlap_1n6)
  # negative beta
  overlap_1n6 <- intersect(negb_lst[[i]]$genes, negb_lst[[i + 1]]$genes)
  write.table(overlap_1n6, 
              file = paste0("overlapped_genes_pval_significant/", 
                            names(negb_lst)[i], "_overlap_w_6hr_genes_downregulated.txt"),
              sep = "\t", quote = F, row.names = F, col.names = F)
  overlap_counts[i, 8] <- length(overlap_1n6)
  overlap_counts[i + 1, 7] <- length(overlap_1n6)
}
write.table(overlap_counts, file = "overlapping_gene_count_summary.csv", 
            quote = F, sep = ",", row.names = T, col.names = T)
