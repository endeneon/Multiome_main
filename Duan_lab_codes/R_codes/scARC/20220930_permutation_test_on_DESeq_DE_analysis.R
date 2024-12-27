# Chuxuan Li 09/30/2022
# permutation test on DE analysis using DESeq2 (after combat)

# init ####

library(DESeq2)
library(Seurat)
library(sva)

library(stringr)
library(future)

library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(ggrepel)
library(reshape2)

load("pseudobulk_DE_08jul2022_cellline_as_covar_dds.RData")
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/pseudobulk_DE/res_use_combatseq_mat/by_line_pseudobulk_mat_adjusted.RData")

# permute ####
permuteDESeq2DEanalysis <- function(mat) {
  colorder <- sample(colnames(mat), ncol(mat), F)
  reordered_mat <- mat[, colorder]
  colnames(reordered_mat) <- colnames(mat)
  coldata <- data.frame(condition = str_extract(string = colnames(reordered_mat),
                                                       pattern = "[0|1|6]"),
                               type = str_extract(string = colnames(reordered_mat),
                                                  pattern = "^CD_[0-9][0-9]"))
  rownames(coldata) <- colnames(reordered_mat)
  dds <- DESeqDataSetFromMatrix(countData = reordered_mat,
                                       colData = coldata,
                                       design = ~ condition)
  keep <- rowSums(counts(dds)) >= 200
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res_0v1 <- results(dds, contrast =c ("condition", "1", "0"))
  return(res_0v1$pvalue)
}

registerDoParallel(cores = 4)
pval_lst_1v0 <- foreach(i=1:100, .combine = c) %dopar% permuteDESeq2DEanalysis(npglut_mat_adj)

# plot results ####
hist(pval_lst_1v0, breaks = 500, col = "grey50", border = "transparent", xlab = "p-value",
     main = "Distribution of p-values from 100 permutations\nof DESeq2 differential expression tests")

res_npglut_0v1 <- results(dds_npglut, contrast =c ("condition", "1", "0"))
hist(res_npglut_0v1$pvalue, breaks = 300, 
     col = "grey50", border = "transparent", xlab = "p-value",
     main = "Distribution of p-values from\noriginal Limma differential expression results (not permuted)")

pval_lst_sample <- -log10(sample(pval_lst_1v0, 10000, F))
unif <- -log10(runif(length(pval_lst_sample)))
qqplot(unif, pval_lst_sample, xlab = "expected -log10(p)", 
       ylab = "observed -log10(p)", cex = 0.1, 
       main = "pseudobulk permuted p-value against uniform distribution")
abline(0, 1, col = "red3")

res_sample <- -log10(sample(res_npglut_0v1$pvalue, 10000, F))
unif <- -log10(sort(runif(length(res_sample))))
qqplot(unif, res_sample, xlab = "expected -log10(p)", 
       ylab = "observed -log10(p)", cex = 0.1, 
       main = "pseudobulk unpermuted p-value against uniform distribution")
abline(0, 1, col = "red3")
