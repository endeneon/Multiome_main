# Chuxuan Li 09/26/2022
# Do permutation test on da peak results from 26Apr2022_call_DA_peaks_on_ATAC_new_obj
#by scrambling the time tags on cells

# init ####
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)

library(future)
library(readr)
library(stringr)
library(foreach)
library(doParallel)
registerDoParallel(cores = 4)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)

set.seed(1886)
# set threads and parallelization
plan("multisession", workers = 1)
options(expressions = 20000)
options(future.globals.maxSize = 11474836480)
load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/ATAC_new_obj_with_MACS2_new_peaks_and_RNA_labels.RData")

# list idents ####
unique(ATAC_new$timextype.ident)
idents_list <- c("GABA_0hr", "GABA_1hr", "GABA_6hr",
                 "NEFM_neg_glut_0hr", "NEFM_neg_glut_1hr", "NEFM_neg_glut_6hr",
                 "NEFM_pos_glut_0hr", "NEFM_pos_glut_1hr", "NEFM_pos_glut_6hr")

ATAC <- subset(ATAC_new, timextype.ident %in% idents_list)
DefaultAssay(ATAC) <- "peaks"

# permute ####
unique(Idents(ATAC))
chrs <- unique(str_extract(rownames(ATAC), "chr[0-9]+-"))
chrs <- chrs[!is.na(chrs)]
peaks <- sample(rownames(ATAC), 10000, F)
peaks.use <- rownames(ATAC)[rownames(ATAC) %in% peaks]
obj <- ATAC[rownames(ATAC) %in% peaks.use, ]
rm(ATAC_new)
Idents(obj) <- "timextype.ident"
permuteScDAPeaks <- function(object) {
  tt <- as.vector(object$timextype.ident)
  tt_reordered <- sample(tt, length(tt), F)
  object$timextype.ident <- tt_reordered
  Idents(object) <- "timextype.ident"
  res <- FindMarkers(
    object = object,
    ident.1 = "NEFM_pos_glut_1hr",
    ident.2 = "NEFM_pos_glut_0hr",
    test.use = 'MAST',
    logfc.threshold = 0, 
    min.pct = 0.01, 
    random.seed = 100)
  return(res$p_val)
}
system.time({pval_lst_1v0 <- foreach(i=1:5, .combine = c) %dopar% permuteScDAPeaks(obj)})
save(pval_lst_1v0, file = "single_cell_da_peak_permutation_pval_lst_seg1.RData")
save(pval_lst_1v0, file = "single_cell_da_peak_permutation_pval_lst_seg2.RData")

load("single_cell_da_peak_permutation_pval_lst_seg1.RData")
pval_lst_1v0_seg1 <- pval_lst_1v0
load("single_cell_da_peak_permutation_pval_lst_seg2.RData")
pval_lst_1v0_seg2 <- pval_lst_1v0
pval_lst_1v0 <- c(pval_lst_1v0_seg1, pval_lst_1v0_seg2)

# plot distribution ####

hist(pval_lst_1v0, breaks = 1000, xlab = "p-value",
     main = "p-value distribution of single-cell da peak results\nfrom 10 permutations on NEFM+ glut")
load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/chr11_peaks_1v0_0.01.RData")
hist(peaklist1[[3]]$p_val, breaks = 500, xlab = "p-value",
     main = "p-value distribution of single-cell original da peak results\n on NEFM+ glut (not permuted)")
hist(peaklist1[[3]]$p_val_adj, 
     main = "p-value distribution of single-cell original da peak results\n on NEFM+ glut (not permuted)")


pval_lst_sample <- -log10(pval_lst_1v0)
unif <- -log10(runif(length(pval_lst_sample)))
qqplot(unif, pval_lst_sample, 
       xlab = "expected -log10(p)", ylab = "observed -log10(p)", cex = 0.1, 
       main = "single cell permuted p-values against uniform distribution")
abline(0, 1, col = "blue")

res_sample <- -log10(sample(peaklist1[[3]]$p_val, 10000, F))
unif <- -log10(sort(runif(length(res_sample))))
qqplot(unif, res_sample, xlab = "expected -log10(p)", ylab = "observed -log10(p)", cex = 0.1, 
       main = "single cell original p-values (not permuted)\nagainst uniform distribution")
abline(0, 1, col = "blue")

# test with fewer iterations, see effect
pval_lst_test <- foreach(i=1:5, .combine = c) %dopar% permuteScDAPeaks(obj)
hist(pval_lst_test, breaks = 500) #not different in terms of how tall the freq at 0 is
