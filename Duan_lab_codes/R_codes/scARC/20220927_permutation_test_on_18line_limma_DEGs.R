# Chuxuan Li 09/27/2022
# permutation test on the limma-generated pseudobulk DE results 

# init ####
library(Seurat)
library(limma)
library(edgeR)
library(sva)

library(readr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

library(doParallel)
library(foreach)

load("../../../covariates_pooled_together_for_all_linextime_5+18+20.RData")
load("./pseudobulk_DE/res_use_combatseq_mat/by_line_pseudobulk_mat_adjusted.RData")
celllines <- str_remove(colnames(GABA_mat_adj), "_[0|1|6]hr$")

covar_18line <- covar_table[covar_table$cell_line %in% lines, ]
covar_18line <- covar_18line[order(covar_18line$cell_line), ]
rm(covar_table)

npglut_fraction <- covar_18line$npglut_fraction

createDGE <- function(count_matrix, lines){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = lines)
  A <- rowSums(y$counts)
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}
filterByCpm <- function(y, cutoff, nsample) {
  cpm <- cpm(y)
  colInd0 <- seq(1, ncol(y), by = 3)
  colInd1 <- seq(2, ncol(y), by = 3)
  colInd6 <- seq(3, ncol(y), by = 3)
  
  cpm_0hr <- cpm[, colInd0]
  cpm_1hr <- cpm[, colInd1]
  cpm_6hr <- cpm[, colInd6]
  #passfilter <- rowSums(zerohr_cpm >= cutoff) >= nsample
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
contrastFit <- function(fit, design) {
  contr.matrix <- makeContrasts(
    onevszero = time1hr-time0hr, 
    sixvszero = time6hr-time0hr, 
    sixvsone = time6hr-time1hr, levels = colnames(design))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit <- eBayes(fit, trend = T)
  print(summary(decideTests(fit)))
  return(fit)
}


# Limma with voom####
design <- model.matrix(~0 + time + group + age + sex + disease + npglut_fraction, data = covar_18line) #29400 genes
permuteLimmaDE <- function(mat, colnumber) {
  colorder <- sample(colnames(mat), ncol(mat), F)
  reordered_mat <- mat[, colorder]
  colnames(reordered_mat) <- colnames(mat)
  reordered_lines <- str_extract(colnames(reordered_mat), "^CD_[0-9][0-9]")
  y <- createDGE(mat, reordered_lines)
  ind.keep <- filterByCpm(y, 1, 9)
  v <- cnfV(y[ind.keep, ], design)
  fit <- lmFit(v, design)
  fit_contr <- contrastFit(fit, design)
  pval_lst <- as.vector(fit_contr$p.value[,colnumber])
  return(pval_lst)
}
registerDoParallel(cores = 4)
pval_lst_1v0_limmaDE <- foreach(i=1:100, .combine = c) %dopar% permuteLimmaDE(npglut_mat_adj, 1)
save(pval_lst_1v0_limmaDE, file = "permuted_pval_lst_from_Limma_DE_analysis.RData")

# plot results ####
hist(pval_lst_1v0, breaks = 500, col = "grey50", border = "transparent", xlab = "p-value",
     main = "Distribution of p-values from 100 permutations\nof Limma differential expression tests")

npglut_res_1v0 <- read_csv("pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm1in9_res/unfiltered_by_padj/npglut_res_1v0_all_DEGs.csv")
hist(npglut_res_1v0$P.Value, breaks = 300, 
     col = "grey50", border = "transparent", xlab = "p-value",
     main = "Distribution of p-values from\noriginal Limma differential expression results (not permuted)")

pval_lst_sample <- -log10(sample(pval_lst_1v0, 10000, F))
unif <- -log10(runif(length(pval_lst_sample)))
qqplot(unif, pval_lst_sample, xlab = "expected -log10(p)", 
       ylab = "observed -log10(p)", cex = 0.1, 
       main = "pseudobulk permuted p-value against uniform distribution")
abline(0, 1, col = "red3")

res_sample <- -log10(sample(npglut_res_1v0$P.Value, 10000, F))
unif <- -log10(sort(runif(length(res_sample))))
qqplot(unif, res_sample, xlab = "expected -log10(p)", 
       ylab = "observed -log10(p)", cex = 0.1, 
       main = "pseudobulk unpermuted p-value against uniform distribution")
abline(0, 1, col = "red3")

# Check with only cell line as covar####
design <- model.matrix(~0 + time + cell_line, data = covar_18line) #29400 genes
permuteLimmaDE <- function(mat, colnumber) {
  colorder <- sample(colnames(mat), ncol(mat), F)
  reordered_mat <- mat[, colorder]
  colnames(reordered_mat) <- colnames(mat)
  reordered_lines <- str_extract(colnames(reordered_mat), "^CD_[0-9][0-9]")
  y <- createDGE(mat, reordered_lines)
  ind.keep <- filterByCpm(y, 1, 9)
  v <- cnfV(y[ind.keep, ], design)
  fit <- lmFit(v, design)
  fit_contr <- contrastFit(fit, design)
  pval_lst <- as.vector(fit_contr$p.value[,colnumber])
  return(pval_lst)
}
registerDoParallel(cores = 4)
pval_lst_1v0_limmaDE <- foreach(i=1:100, .combine = c) %dopar% permuteLimmaDE(npglut_mat_adj, 1)

hist(pval_lst_1v0, breaks = 500, col = "grey50", border = "transparent", xlab = "p-value",
     main = "Distribution of p-values from 100 permutations\nof Limma differential expression tests")

npglut_res_1v0 <- read_csv("pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm1in9_res/unfiltered_by_padj/npglut_res_1v0_all_DEGs.csv")
hist(npglut_res_1v0$P.Value, breaks = 300, 
     col = "grey50", border = "transparent", xlab = "p-value",
     main = "Distribution of p-values from\noriginal Limma differential expression results (not permuted)")

pval_lst_sample <- -log10(sample(pval_lst_1v0_limmaDE, 10000, F))
unif <- -log10(runif(length(pval_lst_sample)))
qqplot(unif, pval_lst_sample, xlab = "expected -log10(p)", 
       ylab = "observed -log10(p)", cex = 0.1, 
       main = "pseudobulk permuted p-value against uniform distribution\nusing cell line as covar")
abline(0, 1, col = "red3")

res_sample <- -log10(sample(npglut_res_1v0$P.Value, 10000, F))
unif <- -log10(sort(runif(length(res_sample))))
qqplot(unif, res_sample, xlab = "expected -log10(p)", 
       ylab = "observed -log10(p)", cex = 0.1, 
       main = "pseudobulk unpermuted p-value against uniform distribution")
abline(0, 1, col = "red3")

# test with different filtering criteria ####
design <- model.matrix(~0 + time + group + age + sex + disease + npglut_fraction, data = covar_18line) #29400 genes
permuteLimmaDE <- function(mat, colnumber) {
  colorder <- sample(colnames(mat), ncol(mat), F)
  reordered_mat <- mat[, colorder]
  colnames(reordered_mat) <- colnames(mat)
  reordered_lines <- str_extract(colnames(reordered_mat), "^CD_[0-9][0-9]")
  y <- createDGE(mat, reordered_lines)
  ind.keep <- filterByCpm(y, 10, 9)
  v <- cnfV(y, design)
  fit <- lmFit(v, design)
  fit_contr <- contrastFit(fit, design)
  pval_lst <- as.vector(fit_contr$p.value[,colnumber])
  return(pval_lst)
}
registerDoParallel(cores = 4)
pval_lst_1v0_limmaDE <- foreach(i=1:100, .combine = c) %dopar% permuteLimmaDE(npglut_mat_adj, 1)
hist(pval_lst_1v0_limmaDE, breaks = 500, col = "grey50", border = "transparent", xlab = "p-value",
     main = "Distribution of p-values from 100 permutations\nof Limma differential expression tests")

npglut_res_1v0 <- read_csv("pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm1in9_res/unfiltered_by_padj/npglut_res_1v0_all_DEGs.csv")
hist(npglut_res_1v0$P.Value, breaks = 300, 
     col = "grey50", border = "transparent", xlab = "p-value",
     main = "Distribution of p-values from\noriginal Limma differential expression results (not permuted)")

pval_lst_sample <- -log10(sample(pval_lst_1v0_limmaDE, 10000, F))
unif <- -log10(runif(length(pval_lst_sample)))
qqplot(unif, pval_lst_sample, xlab = "expected -log10(p)", 
       ylab = "observed -log10(p)", cex = 0.1, 
       main = "pseudobulk permuted p-value against uniform distribution\nno filtering")
abline(0, 1, col = "red3")

res_sample <- -log10(sample(npglut_res_1v0$P.Value, 10000, F))
unif <- -log10(sort(runif(length(res_sample))))
qqplot(unif, res_sample, xlab = "expected -log10(p)", 
       ylab = "observed -log10(p)", cex = 0.1, 
       main = "pseudobulk unpermuted p-value against uniform distribution")
abline(0, 1, col = "red3")

