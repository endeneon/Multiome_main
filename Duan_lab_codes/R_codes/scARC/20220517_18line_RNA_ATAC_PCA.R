# Chuxuan Li 05/17/2022
# plot PCA plot for 18-line ATACseq and RNAseq data, putting cell types together
#in the same plot

# init ####
library(Seurat)
library(Signac)
library(edgeR)

library(sva)
library(FactoMineR)
library(factoextra)

library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(reshape2)


load("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")
RNA_18line <- integrated_labeled
rm(integrated_labeled)
load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/ATAC_new_obj_with_MACS2_new_peaks_and_RNA_labels.RData")
ATAC_new <- subset(ATAC_new, cell.type != "unknown")
ATAC_18line <- ATAC_new
rm(ATAC_new)

# auxiliary functions ####

makeCpmMatCelltypes <- function(obj1, lines_1, types) {
  rownames <- rep_len(NA, (length(lines_1) * length(types)))
  k = 1
  for (i in 1:length(lines_1)) {
    print(lines_1[i])
    rownames[k] <- lines_1[i]
    rownames[k + 1] <- lines_1[i]
    rownames[k + 2] <- lines_1[i]
    obj <- subset(obj1, cell.line.ident == lines_1[i])
    for (l in 1:length(types)) {
      t <- types[l]
      rownames[k] <- paste(rownames[k], t, sep = "_")
      print(rownames[k])
      subtype <- subset(obj, cell.type == t)
      if (i == 1 & l == 1) {
        cpm_mat <- edgeR::cpm(rowSums(subtype@assays$ATAC@counts), normalized.lib.sizes = T, log = F)
        cpm_mat <- cpm_mat[, 1]
      } else {
        to_bind <- edgeR::cpm(rowSums(subtype@assays$ATAC@counts), normalized.lib.sizes = T, log = F)
        to_bind <- to_bind[, 1]
        cpm_mat <- rbind(cpm_mat, to_bind)
      }
      k = k + 1
    }
  }
  rownames(cpm_mat) <- rownames
  return(cpm_mat)
}

makeCpmMatByTime <- function(obj1, lines_1, times) {
  rownames <- rep_len(NA, (length(lines_1) * length(times)))
  k = 1
  for (i in 1:length(lines_1)) {
    print(lines_1[i])
    rownames[k] <- lines_1[i]
    rownames[k + 1] <- lines_1[i]
    rownames[k + 2] <- lines_1[i]
    obj <- subset(obj1, cell.line.ident == lines_1[i])
    for (l in 1:length(times)) {
      t <- times[l]
      rownames[k] <- paste(rownames[k], t, sep = "_")
      print(rownames[k])
      subtime <- subset(obj, time.ident == t)
      if (i == 1 & l == 1) {
        cpm_mat <- edgeR::cpm(rowSums(subtime@assays$ATAC@counts), normalized.lib.sizes = T, log = F)
        cpm_mat <- cpm_mat[, 1]
      } else {
        to_bind <- edgeR::cpm(rowSums(subtime@assays$ATAC@counts), normalized.lib.sizes = T, log = F)
        to_bind <- to_bind[, 1]
        cpm_mat <- rbind(cpm_mat, to_bind)
      }
      k = k + 1
    }
  }
  rownames(cpm_mat) <- rownames
  return(cpm_mat)
}

# RNA ####
lines <- unique(RNA_18line$cell.line.ident)
lines
RNA_18line$cell.type[RNA_18line$cell.type %in% 
                       c("SST_pos_GABA", "GABA", "SEMA3E_pos_GABA")] <-
  "GABA"
RNA_18line <- subset(RNA_18line, cell.type != "unknown")
types <- as.vector(unique(RNA_18line$cell.type))
types
times <- unique(RNA_18line$time.ident)
times


batch <- as.factor(c("12", "12", "12", "8", "8", "8", "8", "9", "9", "9", 
                     "10", "10", "10", "10", "11", "11", "11", "11"))
batch_more <- rep(batch, each = 3)

for (time in times) {
  subtime18 <- subset(RNA_18line, time.ident == time)
  cpm_mat <- makeCpmMatCelltypes(subtime18, lines, types)
  #pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = F)
  com_mat <- ComBat(
    dat = t(cpm_mat),
    batch = batch_more
  )
  pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
  filename <- paste0("ComBat_by_batch_", 
                     time, "_PCA_plot_by_batch.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = T,
                    habillage = batch_more,
                    palette = brewer.pal(name = "Paired", n = 10), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x cell type PCA colored by cocultured batch \n", 
                   time)) 
  print(p)
  dev.off()
  
  filename <- paste0("ComBat_by_batch_", 
                     time, "_PCA_plot_by_exp.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = T,
                    habillage = exp, 
                    palette = brewer.pal(name = "Set1", n = 2), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x cell type PCA colored by 5/18 line \n", time))
  print(p) 
  dev.off()
}

for (type in types) {
  subtype <- subset(RNA_18line, cell.type == type)
  cpm_mat <- makeCpmMatByTime(subtype, lines, times)
  pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = F)
  # com_mat <- ComBat(
  #   dat = t(cpm_mat),
  #   batch = batch_more
  # )
  # pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
  filename <- paste0(#"ComBat_by_batch_", 
                     type, "_PCA_plot_by_batch.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = F,
                    habillage = batch_more,
                    palette = brewer.pal(name = "Paired", n = 10), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x time point PCA colored by cocultured batch \n", 
                   type)) 
  print(p)
  dev.off()
}

# ATAC ####
lines <- sort(as.vector(unique(ATAC_18line$cell.line.ident)))
lines
types <- as.vector(unique(ATAC_18line$cell.type))
types
times <- unique(ATAC_18line$time.ident)
times


batch <- as.factor(c("12", "11", "12", "12", 
                     "8", "8", "9", "9", 
                     "10", "8", "8", "9", 
                     "10", "10", "11", "11", "11", "10"))
batch_more <- rep(batch, each = 3)

for (time in times) {
  subtime18 <- subset(ATAC_18line, time.ident == time)
  cpm_mat <- makeCpmMatCelltypes(subtime18, lines, types)
  #pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = F)
  com_mat <- ComBat(
    dat = t(cpm_mat),
    batch = batch_more
  )
  pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
  filename <- paste0("ComBat_by_batch_", 
    time, "_PCA_plot_by_batch.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = T,
                    habillage = batch_more,
                    palette = brewer.pal(name = "Paired", n = 10), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x cell type PCA colored by cocultured batch \n", 
                   time)) 
  print(p)
  dev.off()
}

for (type in types) {
  subtype <- subset(ATAC_18line, cell.type == type)
  cpm_mat <- makeCpmMatByTime(subtype, lines, times)
  # pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = F)
  com_mat <- ComBat(
    dat = t(cpm_mat),
    batch = batch_more
  )
  pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
  filename <- paste0("ComBat_by_batch_", 
                     type, "_PCA_plot_by_batch_circled.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = T,
                    habillage = batch_more, 
                    palette = brewer.pal(name = "Paired", n = 10), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x time point PCA colored by cocultured batch \n", 
                   type)) 
  print(p)
  dev.off()
}
