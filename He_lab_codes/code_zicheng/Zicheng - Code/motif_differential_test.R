library(tidyverse)
library(lme4)
library(parallel)
library(multcomp)
library(RhpcBLASctl)

blas_set_num_threads(1)
ncores <- 20

union_degs <- readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi_union.rds")

deg_tfs <- readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/combined_trajMM.rds") %>%
  rownames() %>%
  .[. %in% union_degs]

motif_diff_test <- function(cell_type) {
  cat("Started motif differential test for", cell_type, "\n")
  
  load(paste0("/project/xinhe/zicheng/neuron_stimulation/data/scMatrix/combined/scMatrix-", cell_type, ".rda"))
  rm(atac_agg, rna_agg, Peak2GeneList)
  mm_agg <- mm_agg[deg_tfs, ]
  gc()
  
  wilcox_df <- mclapply(rownames(mm_agg), function(tf) {
    test_1hrVS0hr <- wilcox.test(mm_agg[tf, rownames(subset_cellMeta)[subset_cellMeta$Time == "1hr"]],
                                 mm_agg[tf, rownames(subset_cellMeta)[subset_cellMeta$Time == "0hr"]],
                                 alternative = "greater")
    test_6hrVS0hr <- wilcox.test(mm_agg[tf, rownames(subset_cellMeta)[subset_cellMeta$Time == "6hr"]],
                                 mm_agg[tf, rownames(subset_cellMeta)[subset_cellMeta$Time == "0hr"]],
                                 alternative = "greater")
    test_6hrVS1hr <- wilcox.test(mm_agg[tf, rownames(subset_cellMeta)[subset_cellMeta$Time == "6hr"]],
                                 mm_agg[tf, rownames(subset_cellMeta)[subset_cellMeta$Time == "1hr"]],
                                 alternative = "greater")
    
    return(c(test_1hrVS0hr = test_1hrVS0hr$p.value,
             test_6hrVS0hr = test_6hrVS0hr$p.value,
             test_6hrVS1hr = test_6hrVS1hr$p.value))
  }, mc.cores = ncores) %>%
    bind_rows() %>%
    as.data.frame()
  
  rownames(wilcox_df) <- rownames(mm_agg)
  
  return(list(wilcox.p = wilcox_df))
}

motif_diff_pval_list <- lapply(c("GABA", "nmglut", "npglut"), motif_diff_test)
names(motif_diff_pval_list) <- c("GABA", "nmglut", "npglut")

motif_diff_fdr_list <- lapply(motif_diff_pval_list, function(x) {
  lapply(x, function(y) {
    apply(y, 2, function(z) p.adjust(z, method = "fdr")) %>%
      as.data.frame()
  })
})

save(motif_diff_pval_list, motif_diff_fdr_list, file = "../output/Diff_Motif/motif_diff_pval_fdr.RData")

fdr_threshold <- 0.01

# Classify early and late response TF types based on motif changes
early_response_list <- lapply(motif_diff_fdr_list, function(x) {
  lapply(x, function(y) {
    rownames(y)[y[, 1] <= fdr_threshold]
  })
})

late_response_list <- lapply(motif_diff_fdr_list, function(x) {
  lapply(x, function(y) {
    rownames(y)[y[, 2] <= fdr_threshold & y[, 3] <= fdr_threshold]
  })
})

save(early_response_list, late_response_list, file = "../output/Diff_Motif/response_list.RData")

# Define early and late response TFs based on motif changes and expression changes
library(vroom)

setwd("/project/xinhe/zicheng/neuron_stimulation/script/")

deg_res_files <- list.files("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi", full.names = TRUE)

res_list <- lapply(deg_res_files, vroom)
names(res_list) <- gsub("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi/", "", deg_res_files) %>%
  gsub("res_", "", .) %>%
  gsub("_all_DEGs.csv", "", .)

load("../output/Diff_Motif/response_list.RData")

sig_list <- lapply(res_list, function(x) x[x$adj.P.Val <= 0.01, ])

early_response_tfs <- lapply(c("GABA", "nmglut", "npglut"), function(x) {
  res_1vs0 <- sig_list[[paste0(x, "_1v0")]]
  up_1vs0 <- res_1vs0$genes[res_1vs0$logFC > 0]
  
  return(intersect(up_1vs0, early_response_list[[x]]$wilcox.p) %>% sort())
})
names(early_response_tfs) <- c("GABA", "nmglut", "npglut")

late_response_tfs <- lapply(c("GABA", "nmglut", "npglut"), function(x) {
  res_6vs0 <- sig_list[[paste0(x, "_6v0")]]
  up_6vs0 <- res_6vs0$genes[res_6vs0$logFC > 0]
  
  res_6vs1 <- sig_list[[paste0(x, "_6v1")]]
  up_6vs1 <- res_6vs1$genes[res_6vs1$logFC > 0]
  
  return(Reduce(intersect, list(up_6vs0, up_6vs1, late_response_list[[x]]$wilcox.p)) %>% sort())
})
names(late_response_tfs) <- c("GABA", "nmglut", "npglut")

save(early_response_tfs, late_response_tfs, file = "../output/Diff_Motif/response_list_w_expr_filter.RData")

# UpSet plot
library(ComplexHeatmap)

load("../output/Diff_Motif/response_list_w_expr_filter.RData")
m1 <- make_comb_mat(early_response_tfs)
UpSet(m1, set_order = 1:3)

m2 <- make_comb_mat(late_response_tfs)
UpSet(m2, set_order = 1:3)

# Calculate average activity matrix for 9 contexts for plotting
library(edgeR)

load("/project/xinhe/zicheng/neuron_stimulation/data/DEG/Lexi_matrix_4_Zicheng/by_line_pseudobulk_mat_adjusted.RData")

combat_counts_list <- list(GABA_mat_adj, nmglut_mat_adj, npglut_mat_adj)
names(combat_counts_list) <- c("GABA", "nmglut", "npglut")

norm_counts_list <- lapply(combat_counts_list, function(x) {
  dge <- DGEList(x)
  dge <- calcNormFactors(dge)
  return(cpm(dge, log = TRUE))
})
saveRDS(norm_counts_list, "/project/xinhe/zicheng/neuron_stimulation/data/diff_motif/Log2CPM_Expr_Bulk_18lines.rds")

norm_counts_list <- readRDS("/project/xinhe/zicheng/neuron_stimulation/data/diff_motif/Log2CPM_Expr_Bulk_18lines.rds")
avg_expr_contexts <- lapply(names(norm_counts_list), function(x) {
  avg_df <- apply(norm_counts_list[[x]], 1, function(y) {
    avg_0hr <- mean(y[grep("0hr", names(y), fixed = TRUE)])
    avg_1hr <- mean(y[grep("1hr", names(y), fixed = TRUE)])
    avg_6hr <- mean(y[grep("6hr", names(y), fixed = TRUE)])
    
    return(c(avg_0hr, avg_1hr, avg_6hr))
  }) %>%
    t()
  colnames(avg_df) <- paste0(c("0hr", "1hr", "6hr"), "_", x)
  return(avg_df)
}) %>%
  do.call(cbind, .)

avg_expr_contexts["FOS", 2] - avg_expr_contexts["FOS", 1]
avg_expr_contexts["EGR1", 2] - avg_expr_contexts["EGR1", 1]

saveRDS(avg_expr_contexts, "/project/xinhe/zicheng/neuron_stimulation/data/diff_motif/Log2CPM_Avg_Expr_9contexts.rds")

# Average motif activity for 9 contexts
library(Matrix)

avg_motif_contexts <- lapply(c("GABA", "nmglut", "npglut"), function(x) {
  load(paste0("/project/xinhe/zicheng/neuron_stimulation/data/scMatrix/combined/scMatrix-", x, ".rda"))
  rm(atac_agg, rna_agg, Peak2GeneList)
  gc()
  
  avg_df <- apply(mm_agg, 1, function(y) {
    avg_0hr <- mean(y[rownames(subset_cellMeta)[subset_cellMeta$Time == "0hr"]])
    avg_1hr <- mean(y[rownames(subset_cellMeta)[subset_cellMeta$Time == "1hr"]])
    avg_6hr <- mean(y[rownames(subset_cellMeta)[subset_cellMeta$Time == "6hr"]])
    
    return(c(avg_0hr, avg_1hr, avg_6hr))
  }) %>%
    t()
  colnames(avg_df) <- paste0(c("0hr", "1hr", "6hr"), "_", x)
  return(avg_df)
}) %>%
  do.call(cbind, .)

saveRDS(avg_motif_contexts, "/project/xinhe/zicheng/neuron_stimulation/data/diff_motif/Raw_Avg_Motif_9contexts.rds")

