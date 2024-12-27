# Chuxuan Li 09/23/2022
# Do permutation test on 08Aug2022_combat_limma_18line_da_peaks results

# init ####
library(limma)
library(edgeR)
library(sva)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(readr)
library(Signac)
library(Seurat)
library(doParallel)
registerDoParallel(cores = 4)

load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/da_peaks_limma/post-combat_count_matrices.RData")
load("../../covariates_pooled_together_for_all_linextime_5+18+20.RData")

# define covariates from covar table
covar_18line <- covar_table[covar_table$cell_line %in% lines, ]
covar_18line <- covar_18line[order(covar_18line$cell_line), ]
rm(covar_table)
npglut_fraction <- covar_18line$npglut_fraction

# use npglut as an example ####
batch = rep(c("12", "11", "12", "12", "8", "8", "9", "9", "10", "8", "8", "9", 
              "10", "10", "11", "11", "11", "10"), each = 3)
batch
length(batch) #54
ncol(npglut_mat) #54

# auxillary functions ####
createDGE <- function(count_matrix, lines){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = lines)
  A <- rowSums(y$counts)
  #isexpr <- A > 10
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
library(foreach)
design <- model.matrix(~0 + time + group + age + sex + disease + npglut_fraction, data = covar_18line) #29400 genes

permuteLimmaDAPeaks <- function(mat, colnumber) {
  colorder <- sample(colnames(mat), ncol(mat), F)
  reordered_mat <- mat[, colorder]
  colnames(reordered_mat) <- colnames(mat)
  reordered_lines <- str_extract(colnames(reordered_mat), "^CD_[0-9][0-9]")
  y <- createDGE(reordered_mat, lines = reordered_lines)
  ind.keep <- filterByCpm(y, 1, 9)
  v <- cnfV(y[ind.keep, ], design)
  fit <- lmFit(v, design)
  fit_contr <- contrastFit(fit, design)
  pval_lst <- as.vector(fit_contr$p.value[,colnumber])
  return(pval_lst)
}

pval_lst_1v0 <- foreach(i=1:100, .combine = c) %dopar% permuteLimmaDAPeaks(npglut_mat_adj, 1)
save(pval_lst_1v0, file = "pseudobulk_da_peaks_permutation_test_pvals.RData")

# plot p-value distribution ####
hist(pval_lst_1v0, breaks = 1000, col = "grey50", border = "transparent", 
     xlab = "p-value", 
     main = "p-value distribution of pseudobulk da peak results\nfrom 1000 permutations on NEFM+ glut")
sum(pval_lst_1v0 < 0.05)/length(pval_lst_1v0)


npglut_res_1v0 <- read_csv("da_peaks_limma/unfiltered_by_padj/npglut_res_1v0_all_peaks.csv")
npglut_res_1v0$gene <- NULL
colnames(npglut_res_1v0) <- c("peak", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
hist(npglut_res_1v0$P.Value, breaks = 500, col = "grey50", border = "transparent", 
     xlab = "p-value", 
     main = "p-value distribution of pseudobulk da peaks\noriginal results (not permuted)")

# make a uniform distribution
pval_lst_sample <- -log10(sample(pval_lst_1v0, 10000, F))
unif <- -log10(runif(length(pval_lst_sample)))
qqplot(unif, pval_lst_sample, xlab = "expected -log10(p)", 
       ylab = "observed -log10(p)", cex = 0.1, 
       main = "pseudobulk permuted p-value against uniform distribution")
abline(0, 1, col = "red3")

res_sample <- -log10(sample(npglut_res_1v0$P.Value, 10000, F))
unif <- -log10(sort(runif(length(res_sample))))
qqplot(unif, res_sample, xlab = "uniform", ylab = "original p", cex = 0.1, 
       main = "pseudobulk unpermuted p-value against uniform distribution")
abline(0, 1, col = "red3")

# calculate empirical p-value ####
npglut_res_1v0$empirical.pval <- ""
calcEmpiricalP <- function(row) {
  return(sum(pval_lst_1v0 < row[5]) / length(pval_lst_1v0))
}
ep_lst <- apply(npglut_res_1v0, 1, calcEmpiricalP)
