# Chuxuan Li 02/13/2023
# PCA combining all cell line, time points, cell types together for 029

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
library(readr)
library(readxl)

CIRM_table <- read_excel("~/NVME/scARC_Duan_018/CIRM control iPSC lines_40_Duan (003).xlsx")

# clean and make metadata table ####
CIRM_table <- CIRM_table[, c(4, 6, 13, 16, 21)]
CIRM_table <- CIRM_table[!is.na(CIRM_table$`Catalog ID`), ]
colnames(CIRM_table) <- c("group", "cell_line", "age", "sex", "aff")
CIRM_table$sex[CIRM_table$sex == "Female"] <- "F"
CIRM_table$sex[CIRM_table$sex == "Male"] <- "M"
CIRM_table$aff <- "control"
CIRM_table$time <- "0hr"
CIRM_table <- rbind(CIRM_table, CIRM_table, CIRM_table)
CIRM_table$time[53:104] <- "1hr"
CIRM_table$time[105:156] <- "6hr"


# RNA ####
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

load("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/029_RNA_integrated_labeled.RData")

integrated_labeled <- subset(integrated_labeled, cell.type != "unidentified")
types <- unique(integrated_labeled$cell.type)
times <- unique(integrated_labeled$time.ident)
lines <- unique(integrated_labeled$cell.line.ident)
batch <- c()
for (l in lines) {
  print(l)
  batch <- c(batch, unique(CIRM_table$group[CIRM_table$cell_line == l]))
}
batch <- as.factor(batch)
batch_more <- rep(batch, each = 3)

setwd("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/")
for (type in types) {
  subtype <- subset(integrated_labeled, cell.type == type)
  cpm_mat <- makeCpmMatByTime(subtype, lines, times)
  pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = F)
  filename <- paste0("./PCA_plots/by_type-plot_lines_time/raw_data_", 
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
  filename <- paste0("./PCA_plots/by_type-plot_lines_time/ComBat_by_batch_",
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

# ATAC ####
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
        cpm_mat <- edgeR::cpm(rowSums(subtime@assays$peaks@counts), normalized.lib.sizes = T, log = F)
        cpm_mat <- cpm_mat[, 1]
      } else {
        to_bind <- edgeR::cpm(rowSums(subtime@assays$peaks@counts), normalized.lib.sizes = T, log = F)
        to_bind <- to_bind[, 1]
        cpm_mat <- rbind(cpm_mat, to_bind)
      }
      k = k + 1
    }
  }
  rownames(cpm_mat) <- rownames
  return(cpm_mat)
}

load("/nvmefs/scARC_Duan_018/Duan_project_029_ATAC/multiomic_obj_with_new_peaks_labeled.RData")
lines <- sort(as.vector(unique(multiomic_obj_new$cell.line.ident)))
lines
types <- as.vector(unique(multiomic_obj_new$RNA.cell.type))
types <- types[types != "unidentified"]
times <- unique(multiomic_obj_new$time.ident)
times
batch <- c()
for (l in lines) {
  print(l)
  batch <- c(batch, unique(CIRM_table$group[CIRM_table$cell_line == l]))
}
batch <- as.factor(batch)
batch_more <- rep(batch, each = 3)

setwd("/nvmefs/scARC_Duan_018/Duan_project_029_ATAC/")
for (type in types) {
  subtype <- subset(multiomic_obj_new, RNA.cell.type == type)
  cpm_mat <- makeCpmMatByTime(subtype, lines, times)
  pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = F)
  
  filename <- paste0("./PCA_plots/by_type-plot_lines_time/raw_", type, "_PCA_plot_by_type.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = T,
                    habillage = batch_more, 
                    palette = brewer.pal(name = "Paired", n = 10), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("029 ATAC-seq data\nCell lines x time point PCA colored by cocultured batch - ", 
                   type)) 
  print(p)
  dev.off()
  
  com_mat <- ComBat(
    dat = t(cpm_mat),
    batch = batch_more
  )
  pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
  filename <- paste0("./PCA_plots/by_type-plot_lines_time/ComBat_by_batch_", type, "_PCA_plot_by_type.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = T,
                    habillage = batch_more, 
                    palette = brewer.pal(name = "Paired", n = 10), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("029 ATAC-seq data\nCell lines x time point PCA colored by cocultured batch \n combat-corrected data ", 
                   type)) 
  print(p)
  dev.off()
}
