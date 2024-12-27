# Chuxuan Li 10/7/2022
# permutation test on single cell signac-MAST DE results 

# init ####
library(Seurat)
library(SeuratWrappers)

library(dplyr)
library(tidyr)
library(stringr)

library(ggplot2)
library(RColorBrewer)

library(future)
library(doParallel)
library(foreach)
plan("multisession", workers = 2)
options(future.globals.maxSize = 207374182400)

load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")

unique(integrated_labeled$cell.line.ident)
obj <- subset(integrated_labeled, cell.line.ident != "unmatched")
rm(integrated_labeled)

unique(obj$cell.type)
npglut <- subset(obj, cell.type == "NEFM_pos_glut")
npglut <- ScaleData(npglut)
rm(obj)

# compute deg ####
permuteSignacMASTDE <- function(object) {
  t <- as.vector(object$time.ident)
  t_reordered <- sample(t, length(t), F)
  object$time.ident <- t_reordered
  DefaultAssay(test) <- "RNA"
  object <- ScaleData(object)
  Idents(object) <- "time.ident"
  res <- FindMarkers(
    object = object, slot = "scale.data",
    group.by = 'time.ident', 
    ident.1 = "1hr", ident.2 = "0hr",
    test.use = 'MAST', latent.vars = "orig.ident", 
    logfc.threshold = 0, min.pct = 0.01, 
    random.seed = 100)
  res$gene_symbol <- rownames(res)
  return(res$p_val)
}

registerDoParallel(cores = 4)
pval_lst_mast <- foreach(i=1:20, .combine = c) %dopar% permuteSignacMASTDE(npglut)
save(pval_lst_mast, file = "permuted_pval_lst_from_Signac_MAST_DE_analysis.RData")

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
