# Chuxuan Li 05/17/2023
# plot PCA plot for 030 RNAseq data, putting cell types together in the same plot

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
library(stringr)


load("./integrated_labeled.Rdata")

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
        cpm_mat <- edgeR::cpm(rowSums(subtype@assays$RNA@counts), normalized.lib.sizes = T, log = F)
        cpm_mat <- cpm_mat[, 1]
      } else {
        to_bind <- edgeR::cpm(rowSums(subtype@assays$RNA@counts), normalized.lib.sizes = T, log = F)
        to_bind <- to_bind[, 1]
        cpm_mat <- rbind(cpm_mat, to_bind)
      }
      k = k + 1
      print(k)
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
        cpm_mat <- edgeR::cpm(rowSums(subtime@assays$RNA@counts), normalized.lib.sizes = T, log = F)
        cpm_mat <- cpm_mat[, 1]
      } else {
        to_bind <- edgeR::cpm(rowSums(subtime@assays$RNA@counts), normalized.lib.sizes = T, log = F)
        to_bind <- to_bind[, 1]
        cpm_mat <- rbind(cpm_mat, to_bind)
      }
      k = k + 1
    }
  }
  rownames(cpm_mat) <- rownames
  return(cpm_mat)
}


lines <- sort(unique(integrated_labeled$cell.line.ident))
lines
integrated_labeled <- subset(integrated_labeled, cell.type != "unidentified")
types <- as.vector(unique(integrated_labeled$cell.type))
types
times <- unique(integrated_labeled$time.ident)
times
batch <- c()
for (l in lines) {
  print(l)
  batch <- c(batch,
             unique(str_split(integrated_labeled$orig.ident[integrated_labeled$cell.line.ident == l],
                              pattern = "-", n = 2, simplify = T)[1]))
}
batch <- as.factor(batch)
batch_more <- rep(batch, each = 3)

# RNA - 1. look at how cell types cluster ####
for (time in times) {
  subtime18 <- subset(integrated_labeled, time.ident == time)
  cpm_mat <- makeCpmMatCelltypes(subtime18, lines, types)
  pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = F)
  filename <- paste0("./PCA_plots/raw_data_",
                     time, "_PCA_plot_by_batch.pdf")

  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = F,
                    habillage = batch_more,
                    palette = brewer.pal(name = "Paired", n = 10),
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x cell type PCA colored by cocultured batch \n raw data - ",
                   time))
  print(p)
  dev.off()

  com_mat <- ComBat(
    dat = t(cpm_mat),
    batch = batch_more
  )
  pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
  filename <- paste0("./PCA_plots/ComBat_by_batch_",
                     time, "_PCA_plot_by_batch.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = F,
                    habillage = batch_more,
                    palette = brewer.pal(name = "Paired", n = 10),
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x cell type PCA colored by cocultured batch \n combat-corrected data - ",
                   time))
  print(p)
  dev.off()
}

# RNA - 2. look at how times cluster ####
for (type in types) {
  subtype <- subset(integrated_labeled, cell.type == type)
  cpm_mat <- makeCpmMatByTime(subtype, lines, times)
  pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = F)
  filename <- paste0("./PCA_plots/raw_data_",
                     type, "_PCA_plot_by_batch.pdf")

  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = F,
                    habillage = batch_more,
                    palette = brewer.pal(name = "Paired", n = 10),
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x time PCA colored by cocultured batch \n raw data - ",
                   type))
  print(p)
  dev.off()

  com_mat <- ComBat(
    dat = t(cpm_mat),
    batch = batch_more
  )
  pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
  filename <- paste0("./PCA_plots/ComBat_by_batch_",
                     type, "_PCA_plot_by_batch.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = F,
                    habillage = batch_more,
                    palette = brewer.pal(name = "Paired", n = 10),
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x time PCA colored by cocultured batch \n combat-corrected data - ",
                   type))
  print(p)
  dev.off()
}
