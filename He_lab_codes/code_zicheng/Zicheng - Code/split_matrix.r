library(Matrix)
library(tidyverse)
library(vroom)

# Load gene to peak pairs
gene2peakPairs <- readRDS("/project/xinhe/zicheng/neuron_stimulation/data/Gene2Peak_500kb.rds")

# Load and filter cell metadata
cellMeta <- readRDS("/project/xinhe/zicheng/neuron_stimulation/data/cellMeta.rds") %>% 
  as.data.frame() %>% filter(grepl("Batch_024", batch_wo_time))

# Load DEG results and extract selected genes
deg_res_files <- list.files("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi", full.names = TRUE)
res_list <- lapply(deg_res_files, vroom)
names(res_list) <- gsub("/project/xinhe/zicheng/neuron_stimulation/output/gene_modules/DEG_Lexi/", "", deg_res_files) %>%
  gsub("res_", "", .) %>% gsub("_all_DEGs.csv", "", .)

selected_genes <- lapply(res_list, function(x) x$gene) %>% Reduce(union, .)
saveRDS(selected_genes, "/project/xinhe/zicheng/neuron_stimulation/data/scGLM_genes.rds")

# Create combined matrix by chromosome for gene peak mapping
library(Seurat)

# Load motif deviations
motif_deviations <- readRDS("/project/xinhe/zicheng/neuron_stimulation/data/motif_deviations_18lines.rds")

create_scGLM_mtx <- function(cell_type) {
  subset_cellMeta <- cellMeta[cellMeta$cell_type == cell_type, ]
  keep_cells <- rownames(subset_cellMeta)

  rna_list <- list()
  atac_list <- list()

  for (chr in 1:22) {
    rna <- readRDS(paste0("/project/xinhe/zicheng/neuron_stimulation/data/scMatrix/scRNA-", cell_type, "-chr", chr, ".rds"))
    rna <- rna[, keep_cells]

    atac <- readRDS(paste0("/project/xinhe/zicheng/neuron_stimulation/data/scMatrix/scATAC-", cell_type, "-chr", chr, ".rds"))
    atac <- atac[, keep_cells]

    rna_list[[chr]] <- rna
    atac_list[[chr]] <- atac
    rm(rna, atac)
    cat("Finished chr", chr, "\n")
  }

  rna_agg <- do.call(rbind, rna_list)
  atac_agg <- do.call(rbind, atac_list)
  rownames(atac_agg) <- gsub("^(.*?)-", "\\1:", rownames(atac_agg))

  Peak2GeneList <- gene2peakPairs %>% filter(Gene %in% rownames(rna_agg), Peak %in% rownames(atac_agg))
  mm_agg <- motif_deviations[, keep_cells]

  threads <- 20
  output_file <- paste0("/project/xinhe/zicheng/neuron_stimulation/data/scMatrix/combined/scMatrix-", cell_type, ".rda")
  con <- pipe(paste0("/scratch/midway2/zichengwang/software/pigz-2.8/pigz -p", threads, " > ", output_file), "wb")
  save(subset_cellMeta, Peak2GeneList, rna_agg, atac_agg, mm_agg, file = con)
  close(con)

  return(paste0("Finished ", cell_type))
}

lapply(c("GABA", "nmglut", "npglut"), create_scGLM_mtx)

# Create regression input
create_reg_input <- function(cell_type, chr, time = "all", max_row = 5000, threshold = 0.05) {
  if (time == "all") {
    subset_cellMeta <- cellMeta[cellMeta$cell_type == cell_type, ]
    keep_cells <- rownames(subset_cellMeta)
  } else {
    subset_cellMeta <- cellMeta[cellMeta$cell_cluster == paste0(cell_type, "__", time), ]
    keep_cells <- rownames(subset_cellMeta)
  }

  atac <- readRDS(paste0("/project/xinhe/zicheng/neuron_stimulation/data/scMatrix/scATAC-", cell_type, "-chr", chr, ".rds")) %>%
    .[, keep_cells]
  
  # Temporary fix for row names
  rownames(atac) <- gsub("^(.*?)-", "\\1:", rownames(atac))
  atac <- atac[rowMeans(atac > 0) >= threshold, ]

  depth_vec <- subset_cellMeta[colnames(atac), "ReadsInPeaks"]
  scale.factor <- median(depth_vec)
  atac <- t(t(atac) / depth_vec) * scale.factor

  rna <- readRDS(paste0("/project/xinhe/zicheng/neuron_stimulation/data/scMatrix/scRNA-", cell_type, "-chr", chr, ".rds")) %>%
    .[, keep_cells]
  rna <- rna[rownames(rna) %in% selected_genes, ]

  filtered_gene2peakPairs <- gene2peakPairs %>% filter(Gene %in% rownames(rna), Peak %in% rownames(atac))
  gc()

  num_split <- nrow(filtered_gene2peakPairs) %/% max_row + 1
  num_row <- ceiling(nrow(filtered_gene2peakPairs) / num_split)

  split_matrices <- split.data.frame(filtered_gene2peakPairs, f = rep(1:num_split, each = num_row, length.out = nrow(filtered_gene2peakPairs)))
  names(split_matrices) <- paste0("chr", chr, "_split", 1:num_split)

  mclapply(names(split_matrices), function(split) {
    subset_atac <- atac[unique(split_matrices[[split]]$Peak), ]
    subset_rna <- rna[unique(split_matrices[[split]]$Gene), ]
    assoc_obj <- list(Peak2GeneList = split_matrices[[split]], cellMeta = subset_cellMeta, atac_mtx = subset_atac, mrna_mtx = subset_rna)
    save(list = names(assoc_obj), file = paste0("/scratch/midway3/zichengwang/Temp_GenePeakMtx/", time, "_", cell_type, "_", split, ".RData"), envir = list2env(assoc_obj), compress = TRUE)
  }, mc.cores = 12)

  return(split_matrices)
}

reg_batch_list <- list()

for (cell_type in c("GABA", "nmglut", "npglut")) {
  reg_batch <- lapply(1:22, function(chr) create_reg_input(cell_type, chr))
  reg_batch_list[[cell_type]] <- unlist(reg_batch, recursive = FALSE)
  rm(reg_batch)
}

saveRDS(reg_batch_list, "/project/xinhe/zicheng/neuron_stimulation/data/all_Peak2Gene_Batch.rds")

# Generate input file list for snakemake
reg_batch_list <- readRDS("/project/xinhe/zicheng/neuron_stimulation/data/all_Peak2Gene_Batch.rds")
batch_names <- names(unlist(reg_batch_list, recursive = FALSE)) %>% gsub(".", "_", ., fixed = TRUE)
cat(batch_names, file = "/project/xinhe/zicheng/neuron_stimulation/data/all_Peak2Gene_Batch.txt", sep = "\n")
