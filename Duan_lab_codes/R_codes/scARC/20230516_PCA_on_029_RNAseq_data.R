# Chuxuan Li 05/16/2023
# plot PCA plot for 029 RNAseq data, putting cell types together
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
library(stringr)


load("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/029_RNA_integrated_labeled.RData")
integrated_labeled_20line <- integrated_labeled
load("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/integrated_labeled_40201.RData")
integrated_labeled_4line <- integrated_labeled
rm(integrated_labeled)

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
      if (t == "GABA" & lines_1[i] == "CW30154") {
        if (i == 1 & l == 1) {
          cpm_mat <- rep_len(0, nrow(obj))
        } else {
          cpm_mat <- cbind(cpm_mat, rep_len(0, nrow(obj)))
        }
      } else {
        subtype <- subset(obj, cell.type == t)
        if (i == 1 & l == 1) {
          cpm_mat <- edgeR::cpm(rowSums(subtype@assays$RNA@counts), normalized.lib.sizes = T, log = F)
          cpm_mat <- cpm_mat[, 1]
        } else {
          to_bind <- edgeR::cpm(rowSums(subtype@assays$RNA@counts), normalized.lib.sizes = T, log = F)
          to_bind <- to_bind[, 1]
          cpm_mat <- cbind(cpm_mat, to_bind)
        }
      }
      k = k + 1
      print(k)
    }
  }
  colnames(cpm_mat) <- rownames
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
        cpm_mat <- cbind(cpm_mat, to_bind)
      }
      k = k + 1
    }
  }
  colnames(cpm_mat) <- rownames
  return(cpm_mat)
}

# RNA ####
lines20 <- sort(unique(integrated_labeled_20line$cell.line.ident))
lines20
lines4 <- sort(unique(integrated_labeled_4line$cell.line.ident))
lines4
integrated_labeled_20line <- subset(integrated_labeled_20line, cell.type != "unidentified")
integrated_labeled_4line <- subset(integrated_labeled_4line, cell.type != "unidentified")
types <- as.vector(unique(integrated_labeled_20line$cell.type))
types
unique(integrated_labeled_4line$cell.type)
times <- unique(integrated_labeled_20line$time.ident)
times
batch20 <- c()
for (l in lines20) {
  print(l)
  batch20 <- c(batch20,
             unique(str_split(integrated_labeled_20line$orig.ident[integrated_labeled_20line$cell.line.ident == l],
                              pattern = "-", n = 2, simplify = T)[1]))
}
batch20 <- as.factor(batch20)
batch4 <- c()
for (l in lines4) {
  print(l)
  batch4 <- c(batch4,
               unique(str_split(integrated_labeled_4line$orig.ident[integrated_labeled_4line$cell.line.ident == l],
                                pattern = "-", n = 2, simplify = T)[1]))
}
batch4 <- as.factor(batch4)
batch <- c(batch20, batch4)
batch_more <- rep(batch, each = 3)

# RNA - 1. look at how cell types cluster ####
for (time in times) {
  subtime20 <- subset(integrated_labeled_20line, time.ident == time)
  subtime4 <- subset(integrated_labeled_4line, time.ident == time)
  cpm_mat <- makeCpmMatCelltypes(subtime20, lines20, types)
  cpm_mat <- cbind(cpm_mat, makeCpmMatCelltypes(subtime4, lines4, types))
  pca_res <- PCA(t(cpm_mat), scale.unit = TRUE, ncp = 5, graph = F)
  filename <- paste0("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/PCA_plots/raw_data_",
                     time, "_PCA_plot_by_batch_added_40201.pdf")

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
    dat = cpm_mat,
    batch = batch_more
  )
  pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
  filename <- paste0("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/PCA_plots/ComBat_by_batch_",
                     time, "_PCA_plot_by_batch_added_40201.pdf")
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
types <- types[types != "GABA"]
for (type in types) {
  subtime20 <- subset(integrated_labeled_20line, cell.type == type)
  subtime4 <- subset(integrated_labeled_4line, cell.type == type)
  cpm_mat <- makeCpmMatByTime(subtime20, lines20, times)
  cpm_mat <- cbind(cpm_mat, makeCpmMatByTime(subtime4, lines4, times))
  
  pca_res <- PCA(t(cpm_mat), scale.unit = TRUE, ncp = 5, graph = F)
  filename <- paste0("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/PCA_plots/by_type-plot_lines_time_added_40201/raw_data_",
                     type, "_PCA_plot_by_batch_added_40201.pdf")
  
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = F,
                    habillage = batch_more,
                    palette = brewer.pal(name = "Paired", n = 10),
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x time  PCA colored by cocultured batch \n raw data - ",
                   type))
  print(p)
  dev.off()
  
  com_mat <- ComBat(
    dat = cpm_mat,
    batch = batch_more
  )
  pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
  filename <- paste0("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/PCA_plots/by_type-plot_lines_time_added_40201/ComBat_by_batch_",
                     type, "_PCA_plot_by_batch_added_40201.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = F,
                    habillage = batch_more,
                    palette = brewer.pal(name = "Paired", n = 10),
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x time point PCA colored by cocultured batch \n combat-corrected data - ",
                   type))
  print(p)
  dev.off()
}

# type == GABA
subtime4 <- subset(integrated_labeled_4line, cell.type == "GABA")
cpm_mat <- makeCpmMatByTime(subtime20, lines20, times)
cpm_mat <- cbind(cpm_mat, makeCpmMatByTime(subtime4, lines4[2:4], times))

pca_res <- PCA(t(cpm_mat), scale.unit = TRUE, ncp = 5, graph = F)
batch_more <- batch_more[1:69]
filename <- paste0("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/PCA_plots/by_type-plot_lines_time_added_40201/raw_data_",
                   "GABA", "_PCA_plot_by_batch_added_40201.pdf")

  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = F,
                    habillage = batch_more,
                    palette = brewer.pal(name = "Paired", n = 10),
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x time  PCA colored by cocultured batch \n raw data - ",
                   "GABA"))
  print(p)
  dev.off()
  
  com_mat <- ComBat(
    dat = cpm_mat,
    batch = batch_more
  )
  pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
  filename <- paste0("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/PCA_plots/by_type-plot_lines_time_added_40201/ComBat_by_batch_",
                     "GABA", "_PCA_plot_by_batch_added_40201.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = F,
                    habillage = batch_more,
                    palette = brewer.pal(name = "Paired", n = 10),
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x time point PCA colored by cocultured batch \n combat-corrected data - ",
                   "GABA"))
  print(p)
  dev.off()
  