#################### Dynamic QTL input preparation
library(tidyverse)
library(vroom)
library(Matrix)
library(RNOmni)

geno_PCs <- vroom("/project/xinhe/zicheng/neuron_stimulation/data/genotypes/plink_bed/CIRM_v1.2_MGS_merged_100cells_hg38_rsq0.7_maf0.05.eigenvec") %>%
  as.data.frame()
rownames(geno_PCs) <- geno_PCs$`#IID`
geno_PCs <- geno_PCs %>% select(PC1, PC2, PC3)
colnames(geno_PCs) <- paste0("Geno_", colnames(geno_PCs))

logTPM_RankNorm <- function(time, cell_type) {
  count_matrix <- readRDS(paste0("/project/xinhe/zicheng/neuron_stimulation/caQTL/data/count_matrix/archr/Raw/", cell_type, "__", time, ".rds"))
  norm_df <- Seurat::NormalizeData(count_matrix)
  colnames(norm_df) <- paste0(colnames(norm_df), ".", time)
  return(norm_df)
}

gen_dynamic_qtl_input <- function(cell_type) {
  atac_mtx_agg_list <- lapply(c("0hr", "1hr", "6hr"), function(time) logTPM_RankNorm(time, cell_type))
  stopifnot(all(rownames(atac_mtx_agg_list[[1]]) == rownames(atac_mtx_agg_list[[2]])))
  stopifnot(all(rownames(atac_mtx_agg_list[[3]]) == rownames(atac_mtx_agg_list[[2]])))
  
  atac_mtx_agg <- bind_cols(atac_mtx_agg_list) %>% apply(., 1, RankNorm) %>% t()
  rownames(atac_mtx_agg) <- rownames(atac_mtx_agg_list[[1]])
  
  agg_pca <- prcomp(atac_mtx_agg, scale. = TRUE, center = TRUE)$rotation[, 1:10]
  sample_id <- gsub('\\..*', "", rownames(agg_pca))
  combined_covar <- cbind(geno_PCs[sample_id, ], agg_pca) %>% scale() %>% t()
  colnames(combined_covar) <- rownames(agg_pca)
  
  save(atac_mtx_agg, combined_covar, sample_id, file = paste0("/project/xinhe/zicheng/neuron_stimulation/caQTL/data/LMM/", cell_type, "_Pheno_Covar.RData"))
  print(paste("Finished", cell_type))
  return("Done!")
}

lapply(c("GABA", "nmglut", "npglut"), gen_dynamic_qtl_input)

########## Generate plotting counts
#################### Dynamic QTL input prep
######## Generate plot object
library(tidyverse)
library(vroom)
library(Matrix)

time_vector <- rep(c("0hr", "1hr", "6hr"), 3)
cell_type_vector <- c(rep("GABA", 3), rep("nmglut", 3), rep("npglut", 3))
clusters <- paste0(time_vector, "__", cell_type_vector)
clusters_rev <- paste0(cell_type_vector, "__", time_vector)

logTPM_counts <- function(cluster) {
  count_matrix <- readRDS(paste0("/project/xinhe/zicheng/neuron_stimulation/caQTL/data/count_matrix/archr/Raw/", cluster, ".rds"))
  count_matrix <- Matrix(as.matrix(count_matrix), sparse = TRUE)
  norm_df <- Seurat::NormalizeData(count_matrix)
  return(norm_df)
}

logTPM_plot <- lapply(clusters_rev, logTPM_counts)
names(logTPM_plot) <- clusters

saveRDS(logTPM_plot, "/project/xinhe/zicheng/neuron_stimulation/caQTL/data/count_matrix/TMM_plot/logTPM_List_for_Plot.rds")

################ Dynamic QTL tests with permutations
library(MASS)
library(RhpcBLASctl)
library(tidyverse)
library(vroom)
library(data.table)
library(parallel)
library(qvalue)

blas_set_num_threads(1)

load("/project/xinhe/zicheng/neuron_stimulation/caQTL/output/caQTL_mapping/res_analysis/cPeak_list_25kb.RData")

topQTL_list <- lapply(c("GABA", "nmglut", "npglut"), function(cell_type) {
  cPeak_list[grep(cell_type, names(cPeak_list))] %>%
    bind_rows() %>%
    group_by(phenotype_id) %>%
    filter(pval_true_df == min(pval_true_df)) %>%
    ungroup() %>%
    select(phenotype_id, variant_id) %>%
    as.data.frame()
})

names(topQTL_list) <- c("GABA", "nmglut", "npglut")

unique_variant <- lapply(topQTL_list, function(x) x$variant_id) %>% Reduce(union, .)

snp_pos <- vroom("/project2/xinhe/lifanl/neuron_stim/geno_info_100lines.txt", col_names = c("chr", "pos", "a1", "a2", "rsid")) %>%
  as.data.frame()
rownames(snp_pos) <- snp_pos$rsid
snp_pos <- snp_pos[unique_variant, ]

snp <- data.frame(fread("/project/xinhe/lifan/neuron_stim/mateqtl_input_100lines/snps.txt"), row.names = 1) %>%
  .[unique_variant, ] %>% as.matrix()

permutation_LM <- function(iter, QTL_list, atac_mtx, snp_mtx, time_col, covar_mtx) {
  curr_snp <- QTL_list[iter, "variant_id"]
  curr_pheno <- QTL_list[iter, "phenotype_id"]
  
  reg_df <- data.frame(atac = atac_mtx[curr_pheno, ], snp = snp_mtx[curr_snp, ], time = time_col, covar_mtx)
  
  mod <- lm(atac ~ snp + time + snp:time + Geno_PC1*time + Geno_PC2*time + Geno_PC3*time +
              PC1*time + PC2*time + PC3*time + PC4*time + PC5*time, data = reg_df)
  reg_coef <- coef(summary(mod))[c("snp", "snp:time1hr", "snp:time6hr"), ]
  
  mod0 <- lm(atac ~ snp + time + Geno_PC1*time + Geno_PC2*time + Geno_PC3*time +
               PC1*time + PC2*time + PC3*time + PC4*time + PC5*time, data = reg_df)
  anova_res <- anova(mod0, mod)
  
  perm_pvals <- lapply(1:1000, function(x) {
    perm_df <- reg_df %>% group_by(sample_id) %>% mutate(perm_time = sample(time, size = n(), replace = FALSE)) %>% ungroup() %>%
      as.data.frame()
    
    perm_mod <- lm(atac ~ snp + time + snp:perm_time + Geno_PC1*time + Geno_PC2*time + Geno_PC3*time +
                     PC1*time + PC2*time + PC3*time + PC4*time + PC5*time, data = perm_df)
    perm_coef <- coef(summary(perm_mod))[c("snp:perm_time1hr", "snp:perm_time6hr"), 4]
    perm_anova_res <- anova(mod0, perm_mod)
    return(list(lm_perm = perm_coef, anova_perm = perm_anova_res$`Pr(>F)`[2]))
  })
  
  perm_1hr <- lapply(perm_pvals, function(x) x$lm_perm["snp:perm_time1hr"]) %>% unlist()
  perm_6hr <- lapply(perm_pvals, function(x) x$lm_perm["snp:perm_time6hr"]) %>% unlist()
  perm_anova <- lapply(perm_pvals, function(x) x$anova_perm) %>% unlist()
  
  emp_1hr <- est_perm_pval(reg_coef["snp:time1hr", 4], perm_1hr)
  emp_6hr <- est_perm_pval(reg_coef["snp:time6hr", 4], perm_6hr)
  emp_anova <- est_perm_pval(anova_res$`Pr(>F)`[2], perm_anova)
  emp_pvals <- c(emp_1hr, emp_6hr, emp_anova)
  names(emp_pvals) <- c("emp_1hr_pval", "emp_6hr_pval", "emp_anova_pval")
  
  return(list(coef = reg_coef, anova_res = anova_res, emp_pvals = emp_pvals))
}

est_perm_pval <- function(obs_pval, perm_pvals) {
  n_perm <- length(perm_pvals)
  emp_pval <- (sum(perm_pvals <= obs_pval) + 1) / (n_perm + 1)
  return(c(emp_pval = emp_pval))
}

# Run permutation tests
run_dyn_test <- function(cell_type) {
  load(paste0("/project/xinhe/zicheng/neuron_stimulation/caQTL/data/LMM/", cell_type, "_Pheno_Covar.RData"))
  
  atac_mtx_agg <- atac_mtx_agg[unique(topQTL_list[[cell_type]]$phenotype_id), ]
  stopifnot(all(gsub('\\..*', "", colnames(atac_mtx_agg)) == gsub('\\..*', "", colnames(combined_covar))))
  
  subset_snp <- snp[unique(topQTL_list[[cell_type]]$variant_id), sample_id]
  colnames(subset_snp) <- colnames(atac_mtx_agg)
  
  time_points <- gsub(".*\\.", "", colnames(atac_mtx_agg))
  
  combined_covar <- t(combined_covar) %>% as.data.frame()
  combined_covar$sample_id <- sub("\\..*hr", "", rownames(combined_covar))
  
  gc()
  
  dyn_res <- mclapply(1:nrow(topQTL_list[[cell_type]]), function(iter) permutation_LM(iter, topQTL_list[[cell_type]], atac_mtx_agg, subset_snp, time_points, combined_covar), mc.cores = 38)
  return(dyn_res)
}

dyn_test_res_list <- lapply(c("GABA", "nmglut", "npglut"), run_dyn_test)
names(dyn_test_res_list) <- c("GABA", "nmglut", "npglut")

save(topQTL_list, dyn_test_res_list, file = "/project/xinhe/zicheng/neuron_stimulation/caQTL/output/dynamic_QTL/final/LM_Perm_output.RData")

### Analyze results
library(tidyverse)
library(qqman)
library(qvalue)

load("/project/xinhe/zicheng/neuron_stimulation/caQTL/output/dynamic_QTL/final/LM_Perm_output.RData")

dynamic_qtl_perm <- lapply(names(dyn_test_res_list), function(x) {
  lapply(dyn_test_res_list[[x]], function(y) y$emp_pvals) %>% bind_rows() %>%
    bind_cols(topQTL_list[[x]], .)
})
names(dynamic_qtl_perm) <- names(dyn_test_res_list)

dynamic_qtl_fdr <- lapply(dynamic_qtl_perm, function(x) {
  fdr_df <- apply(x[, 3:5], 2, function(y) qvalue(y)$qvalues)
  colnames(fdr_df) <- c("emp_1hr_qval", "emp_6hr_qval", "emp_anova_qval")
  cbind(x[, 1:2], fdr_df)
})

dynamic_cPeaks <- lapply(dynamic_qtl_fdr, function(x) {
  x %>% filter(emp_anova_qval <= 0.05) %>% .$phenotype_id
})

dyn_agg_df <- lapply(c("GABA", "nmglut", "npglut"), function(cell_type) {
  lapply(dyn_test_res_list[[cell_type]], function(x) x$coef[c("snp", "snp:time1hr", "snp:time6hr"), 4]) %>% bind_rows() %>%
    bind_cols(topQTL_list[[cell_type]], .) %>%
    mutate(qval_0hr = qvalue(snp)$qvalues,
           qval_1hr = qvalue(`snp:time1hr`)$qvalues,
           qval_6hr = qvalue(`snp:time6hr`)$qvalues)
})
names(dyn_agg_df) <- c("GABA", "nmglut", "npglut")

anova_agg_df <- lapply(c("GABA", "nmglut", "npglut"), function(cell_type) {
  lapply(dyn_test_res_list[[cell_type]], function(x) x$anova_res$`Pr(>F)`[2]) %>% unlist() %>%
    bind_cols(topQTL_list[[cell_type]], anova_p = .) %>%
    mutate(anova_qval = qvalue(anova_p)$qvalues)
})

qq(anova_agg_df[[1]]$anova_p)

save(dyn_agg_df, anova_agg_df, file = "/project/xinhe/zicheng/neuron_stimulation/caQTL/output/LMM/final/LM_dynamicQTL_All.RData")
