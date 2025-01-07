library(ArchR)
library(tidyverse)

setwd("/project/xinhe/zicheng/neuron_stimulation/script")
addArchRThreads(8)

# Load ArchR project
ArchR_subset <- loadArchRProject("/scratch/midway3/zichengwang/ArchR_subset/")

# Initialize list to store trajectory gene expression matrices
trajGEX_list <- list()
for (cellType in c("GABA", "nmglut", "npglut")) {
  traj <- paste0(cellType, "U")
  
  trajGEX <- getTrajectory(
    ArchRProj = ArchR_subset,
    name = traj,
    useMatrix = "GeneExpressionMatrix",
    log2Norm = TRUE
  )
  
  trajGEX_mat <- assay(trajGEX, "mat")
  colnames(trajGEX_mat) <- paste0(traj, ".", colnames(trajGEX_mat))
  trajGEX_list[[traj]] <- trajGEX_mat
}

# Check dimensions and row names consistency
lapply(trajGEX_list, dim)
all(rownames(trajGEX_list[[1]]) == rownames(trajGEX_list[[2]]))
all(rownames(trajGEX_list[[1]]) == rownames(trajGEX_list[[3]]))

# Combine trajectory gene expression matrices
combined_trajGEX <- do.call(cbind, trajGEX_list)
saveRDS(combined_trajGEX, "/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/combined_trajGEX.rds")

# Load DEG results
deg_res_files <- list.files("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi", full.names = TRUE)
res_list <- lapply(deg_res_files, vroom)
names(res_list) <- gsub("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi/", "", deg_res_files) %>%
  gsub("res_", "", .) %>%
  gsub("_all_DEGs.csv", "", .)

# Filter significant DEGs
sig_list <- lapply(res_list, function(x) x[x$adj.P.Val <= 0.05 & abs(x$logFC) >= 1, ])
union_degs <- lapply(sig_list, function(x) x$gene) %>% Reduce(union, .)
saveRDS(union_degs, "/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi_union_logFC1_Padj0.05.rds")

# Load combined trajectory gene expression matrix and union DEGs
combined_trajGEX <- readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/combined_trajGEX.rds")
rownames(combined_trajGEX) <- gsub(".*:", "", rownames(combined_trajGEX))

union_degs <- readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi_union_logFC1_Padj0.05.rds")
union_degs <- intersect(union_degs, rownames(combined_trajGEX))

top_trajGEX <- combined_trajGEX[union_degs, ]

# Scale gene expression data
scaled_top_trajGEX <- lapply(list(1:100, 101:200, 201:300), function(time_vec) {
  t(apply(top_trajGEX[, time_vec], 1, rescale))
}) %>% do.call(cbind, .)

# Perform k-means clustering
k <- 15
set.seed(42)
km_res <- kmeans(scaled_top_trajGEX, k, iter.max = 50, nstart = 20)
saveRDS(km_res, "/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/kmeans_k15_logFC1_unionDEGs.rds")

# Save clustering results
write_excel_csv(
  data.frame(Gene = names(km_res$cluster), Cluster = km_res$cluster) %>% arrange(Cluster),
  "/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/kmeans_k15_logFC1_unionDEGs.csv"
)

# Enrichment analysis
library(enrichR)
library(openxlsx)

km_res <- readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/kmeans_k15_logFC1_unionDEGs.rds")
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "DisGeNET", "PanglaoDB_Augmented_2021")
enriched_list <- lapply(1:15, function(x) enrichr(names(km_res$cluster)[km_res$cluster == x], dbs))
saveRDS(enriched_list, "/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/kmeans_k15_logFC1_enriched_list.rds")

enriched_list <- readRDS("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/kmeans_k15_logFC1_enriched_list.rds")
enrich_datasets <- names(enriched_list[[1]])

# Save enrichment results to Excel files
for (dat in enrich_datasets) {
  xl_lst <- list()
  for (i in 1:k) {
    xl_lst[[paste0("Cluster ", i)]] <- as.data.frame(enriched_list[[i]][[dat]])
  }
  write.xlsx(xl_lst, file = paste0("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/enrichment/", dat, "_enrichment.xlsx"))
}

# Filter and save significant enrichment results
enrich_res <- list()
for (dat in enrich_datasets) {
  dat_list <- list()
  for (i in 1:k) {
    dat_list[[paste0("Cluster ", i)]] <- as.data.frame(enriched_list[[i]][[dat]]) %>%
      filter(Adjusted.P.value <= 0.1) %>%
      mutate(Cluster = i)
  }
  enrich_res[[dat]] <- bind_rows(dat_list) %>% arrange(-Combined.Score)
  write_excel_csv(enrich_res[[dat]], paste0("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/enrichment/sig_res/", dat, "sig_Padj0.1_enrichment.csv"))
}


