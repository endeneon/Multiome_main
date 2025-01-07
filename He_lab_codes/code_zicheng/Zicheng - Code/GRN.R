# Load necessary libraries
library(tidyverse)
library(Matrix)

# Set working directory
setwd("/project/xinhe/zicheng/neuron_stimulation/script/")

# Load data
union_degs <- readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi_union.rds")
combined_trajGEX <- readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/combined_trajGEX.rds")
combined_trajMM <- readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/combined_trajMM.rds")

# Preprocess data
combined_trajGEX <- combined_trajGEX[rowSums(combined_trajGEX) > 0, ]
combined_trajGEX <- t(apply(combined_trajGEX, 1, scale))
rownames(combined_trajGEX) <- gsub(".*:", "", rownames(combined_trajGEX))

combined_trajMM <- t(apply(combined_trajMM, 1, scale))

# Find intersected transcription factors (TFs)
intersected_tfs <- intersect(rownames(combined_trajGEX), rownames(combined_trajMM))
intersected_tfs <- intersect(intersected_tfs, union_degs)

# Calculate self-correlated TFs
tf_self_corr <- lapply(c("GABA", "nmglut", "npglut"), function(cell_type) {
  sapply(intersected_tfs, function(tf) {
    sel_col <- switch(cell_type,
                      "GABA" = 1:100,
                      "nmglut" = 101:200,
                      "npglut" = 201:300)
    cor(combined_trajGEX[tf, sel_col], combined_trajMM[tf, sel_col], method = "spearman")
  })
})

# Filter TFs with correlation > 0.3
tf_regulator_list <- lapply(tf_self_corr, function(x) names(x)[x > 0.3])
names(tf_regulator_list) <- c("GABA", "nmglut", "npglut")

# Save results
saveRDS(tf_regulator_list, "../output/self_cor_tfs_cellType.rds")

# Load additional libraries
library(SummarizedExperiment)
library(ArchR)

# Load TF regulator list
tf_regulator_list <- readRDS("../output/self_cor_tfs_cellType.rds")

# Set working directory and load ArchR project
setwd("/project/xinhe/zicheng/neuron_stimulation/script/")
addArchRThreads(12)
ArchR_subset <- loadArchRProject("/scratch/midway3/zichengwang/ArchR_subset/")
peakset <- paste0(ArchR_subset@peakSet@seqnames, ":", ArchR_subset@peakSet@ranges)

# Load peak-to-gene results
p2g_res_list <- readRDS("../output/scmixedGLM/log_Mito0.01_p2g_res_list.rds")

# Filter significant peak-gene pairs
sig_peak_gene_pairs <- lapply(p2g_res_list, function(x) {
  x %>% filter(fdr <= 0.1, beta > 0, Gene %in% union_degs) %>% select(Gene, Peak)
})

sig_peak_gene_pairs_list <- lapply(sig_peak_gene_pairs, function(x) split(x$Peak, x$Gene))

# Load motif matrix
motif_mtx <- readRDS("/scratch/midway3/zichengwang/ArchR_subset/Annotations/Motif-Matches-In-Peaks.rds") %>% assay()
colnames(motif_mtx) <- gsub("_.*", "", colnames(motif_mtx))
rownames(motif_mtx) <- peakset
motif_mtx <- motif_mtx[unique(unlist(sig_peak_gene_pairs_list)), unique(unlist(tf_regulator_list))]

# Load parallel processing libraries
library(parallel)
library(RhpcBLASctl)
blas_set_num_threads(1)

# Identify TF-gene test pairs
tf_gene_test_pair_list <- lapply(c("GABA", "nmglut", "npglut"), function(cell_type) {
  tf_gene_test_pairs <- mclapply(sig_peak_gene_pairs_list[[cell_type]], function(peaks) {
    subset_motif_mtx <- motif_mtx
    colnames(subset_motif_mtx)[colSums(subset_motif_mtx[peaks, , drop = FALSE]) > 0]
  }, mc.cores = 20)
  tf_gene_test_pairs <- tf_gene_test_pairs[lapply(tf_gene_test_pairs, length) > 0]
  tf_gene_test_pairs
})

names(tf_gene_test_pair_list) <- c("GABA", "nmglut", "npglut")

# Calculate correlations
corr_res_list <- lapply(c("GABA", "nmglut", "npglut"), function(cell_type) {
  sel_col <- switch(cell_type,
                    "GABA" = 1:100,
                    "nmglut" = 101:200,
                    "npglut" = 201:300)
  corr_res_cellType <- mclapply(names(tf_gene_test_pair_list[[cell_type]]), function(gene) {
    corr_res <- lapply(tf_gene_test_pair_list[[cell_type]][[gene]], function(tf) {
      cor(combined_trajGEX[gene, sel_col], combined_trajMM[tf, sel_col], method = "spearman")
    }) %>% unlist()
    names(corr_res) <- tf_gene_test_pair_list[[cell_type]][[gene]]
    corr_res
  }, mc.cores = 20)
  names(corr_res_cellType) <- names(tf_gene_test_pair_list[[cell_type]])
  corr_res_cellType
})

names(corr_res_list) <- c("GABA", "nmglut", "npglut")

# Create correlation data frame
corr_df <- lapply(corr_res_list, function(x) {
  imap_dfr(x, ~ tibble(tf = names(.x), gene = .y, correlation = .x))
})

# Save correlation results
saveRDS(corr_df, "../output/scmixedGLM_GRN/log_Mito0.01_eGRN_corr_res_full.rds")

# Filter significant gene regulatory networks (GRNs)
sig_grn <- lapply(corr_df, function(x) x %>% filter(abs(correlation) > 0.5))
saveRDS(sig_grn, "../output/scmixedGLM_GRN/log_Mito0.01_cellType_eGRN_corr_res_tf0.2_corr0.5.rds")

# Save significant GRNs to TSV files
lapply(names(sig_grn), function(cell_type) {
  write_tsv(sig_grn[[cell_type]], paste0("../output/scmixedGLM_GRN/log_Mito0.01_", cell_type, "_corr_res_tf0.2_corr0.5.tsv"))
})

# Load additional libraries for ASD enrichment analysis
library(tidyverse)
library(parallel)

# Load significant GRNs and ASD genes
sig_grn <- readRDS("/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/log_Mito0.01_cellType_eGRN_corr_res_tf0.2_corr0.5.rds")
asd_genes <- scan("/project/xinhe/zicheng/neuron_stimulation/data/ASD_Fu_et_al_0.05.txt", character())
union_degs <- readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi_union.rds")

# Perform ASD enrichment analysis
asd_enriched_tfs <- mclapply(c("GABA", "nmglut", "npglut"), function(cell_type) {
  gene_targets <- sig_grn[[cell_type]]
  gene_targets <- split(gene_targets$gene, gene_targets$tf)
  gene_targets <- lapply(gene_targets, unique)
  
  asd_match <- union_degs %in% asd_genes
  
  asd_enrichment <- lapply(names(gene_targets), function(tf) {
    tf_targets_match <- union_degs %in% gene_targets[[tf]]
    res_fisher <- fisher.test(table(asd_match, tf_targets_match))
    res_vector <- c(res_fisher$estimate, res_fisher$p.value)
    names(res_vector) <- c("odd_ratio", "p_value")
    res_vector
  }) %>% bind_rows() %>% cbind(TF = names(gene_targets), .) %>% mutate(FDR = p.adjust(p_value, method = "fdr"))
  
  sig_asd_enrichment <- asd_enrichment %>% filter(odd_ratio > 1, FDR <= 0.05) %>% arrange(FDR)
  sig_asd_enrichment
}, mc.cores = 3)

names(asd_enriched_tfs) <- c("GABA", "nmglut", "npglut")

# Save ASD enrichment results to Excel file, Table S21
library(openxlsx)
write.xlsx(asd_enriched_tfs, "/project/xinhe/zicheng/neuron_stimulation/output/scmixedGLM_GRN/ASD_Enriched_TFs_by_CellType.xlsx")

