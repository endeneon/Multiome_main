# Chuxuan Li 03/13/2023
# Regression with PRS as the main variable of interest, regress out time, aff, etc.

# init ####
library(Seurat)
library(readr)
library(readxl)
library(stringr)
library(limma)
library(edgeR)

load("covar_table_only_lines_w_PRS.RData")
load("sep_by_type_col_by_linextime_adj_mat_lst_w_only_lines_w_PRS.RData")

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
  colInd0 <- str_detect(colnames(y), "0hr")
  colInd1 <- str_detect(colnames(y), "1hr")
  colInd6 <- str_detect(colnames(y), "6hr")
  
  cpm_0hr <- cpm[, colInd0]
  cpm_1hr <- cpm[, colInd1]
  cpm_6hr <- cpm[, colInd6]
  
  passfilter <- (rowSums(cpm_0hr >= cutoff) >= nsample |
                   rowSums(cpm_1hr >= cutoff) >= nsample |
                   rowSums(cpm_6hr >= cutoff) >= nsample)
  return(passfilter)
}
cnfV <- function(y, design) {
  y <- calcNormFactors(y)
  # voom (no logcpm)
  v <- voom(y, design, plot = T)
  
  return(v)  
}

# run limma by cell type ####
design_GABA <- model.matrix(~PRS + time + batch + age + sex + GABA_fraction, 
                            data = covar_table_PRS)
y_GABA <- createDGE(adj_mat_lst$GABA) # plug in adjusted or raw matrix
ind.keep <- filterByCpm(y_GABA, 1, 24)
sum(ind.keep) #17775
v_GABA <- cnfV(y_GABA[ind.keep, ], design_GABA)
fit_GABA <- lmFit(v_GABA, design_GABA)
fit_GABA <- eBayes(fit_GABA)
topTable(fit_GABA, coef = ncol(design_GABA))

design_nmglut <- model.matrix(~PRS + time + batch + age + sex + nmglut_fraction, 
                            data = covar_table_PRS)
y_nmglut <- createDGE(adj_mat_lst$nmglut) # plug in adjusted or raw matrix
ind.keep <- filterByCpm(y_nmglut, 1, 24)
sum(ind.keep) #17272
v_nmglut <- cnfV(y_nmglut[ind.keep, ], design_nmglut)
fit_nmglut <- lmFit(v_nmglut, design_nmglut)
fit_nmglut <- eBayes(fit_nmglut)
topTable(fit_nmglut, coef = ncol(fit_nmglut))

design_npglut <- model.matrix(~PRS + time + batch + age + sex + npglut_fraction, 
                              data = covar_table_PRS)
y_npglut <- createDGE(adj_mat_lst$npglut) # plug in adjusted or raw matrix
ind.keep <- filterByCpm(y_npglut, 1, 24)
sum(ind.keep) #17378
v_npglut <- cnfV(y_npglut[ind.keep, ], design_npglut)
fit_npglut <- lmFit(v_npglut, design_npglut)
fit_npglut <- eBayes(fit_npglut)
topTable(fit_npglut, coef = ncol(fit_npglut))

fit_all <- list(fit_GABA, fit_nmglut, fit_npglut)
setwd("./sep_by_celltype_DE_results/")
types <- names(adj_mat_lst)
for (i in 1:3){
  res <- topTable(fit_all[[i]], coef = "PRS", p.value = 0.05, sort.by = "P", number = Inf)
  res$gene <- rownames(res)
  
  filename <- paste0("./upregulated_significant/", types[i], "_PRS_upregulated_DEGs.csv")
  print(filename)
  write.table(res[res$logFC > 0, ], file = filename, quote = F, sep = ",", 
              row.names = F, col.names = T)
  filename <- paste0("./downregulated_significant/", types[i], "_PRS_downregulated_DEGs.csv")
  print(filename)
  write.table(res[res$logFC < 0, ], file = filename, quote = F, sep = ",", 
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
cnames <- c("up", "down", "total passed filter")
deg_counts <- array(dim = c(3, 3), dimnames = list(types, cnames)) # 3 cell types are rows, 1/6 up/down are cols
for (i in 1:3){                                                               
  res <- topTable(fit_all[[i]], coef = "PRS", p.value = Inf, sort.by = "P", number = Inf)
  res_6v0 <- topTable(fit_all[[i]], coef = "PRS", p.value = Inf, sort.by = "P", number = Inf)
  deg_counts[i, 1] <- sum(res$logFC > 0 & res$adj.P.Val < 0.05)
  deg_counts[i, 2] <- sum(res$logFC < 0 & res$adj.P.Val < 0.05)
  deg_counts[i, 3] <- nrow(res)
}
write.table(deg_counts, file = "./PRS_assoc_gene_count_summary.csv", quote = F, sep = ",", 
            row.names = T, col.names = T)
