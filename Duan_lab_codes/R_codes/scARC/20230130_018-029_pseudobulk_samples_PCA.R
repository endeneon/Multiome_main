# Chuxuan Li 01/30/2023
# look at correlation between all batches of data using PCA

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

load("/nvmefs/scARC_Duan_018/case_control_DE_analysis/018-029_covar_table_final.RData")
load("~/NVME/scARC_Duan_018/case_control_DE_analysis/018-029_by_type_and_time_mat_lst_bytime.RData")

# separated by cell type and time point ####
celllines <- str_remove(colnames(mat_lst_bytime[[1]]), "_[0|1|6]hr$")
for (i in 1:length(celllines)) {
  print(celllines[i])
  if (i == 1) {
    batch = unique(covar_table_final$batch[covar_table_final$cell_line == celllines[i]])
  } else {
    batch = c(batch, unique(covar_table_final$batch[covar_table_final$cell_line == celllines[i]]))
  }
}
batch <- as.factor(batch)

for (i in 1:length(celllines)) {
  print(celllines[i])
  if (i == 1) {
    seqbatch = unique(covar_table_final$seq_batch[covar_table_final$cell_line == celllines[i]])
  } else {
    seqbatch = c(seqbatch, unique(covar_table_final$seq_batch[covar_table_final$cell_line == celllines[i]]))
  }
}
seqbatch <- as.factor(seqbatch)

for (i in 1:length(mat_lst_bytime)) {
  type <- str_split(names(mat_lst_bytime)[i], "_", 2, T)[,1]
  time <- str_split(names(mat_lst_bytime)[i], "_", 2, T)[,2]
  cpm_mat <- cpm(mat_lst_bytime[[i]])
  #pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = F)
  com_mat <- ComBat(
    dat = cpm_mat,
    batch = batch
  )
  pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
  filename <- paste("./PCA_plots/ComBat_by_batch", type, time, "PCA_plot_by_batch.pdf",
                    sep = "_")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = T,
                    habillage = batch,
                    palette = rainbow(20), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines PCA colored by co-cultured batch \n", 
                   type, " - ", time)) 
  print(p)
  dev.off()
  
  filename <- paste("./PCA_plots/ComBat_by_batch", type, time, "PCA_plot_by_sequencing_batch.pdf",
                    sep = "_")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = T,
                    habillage = seqbatch, 
                    palette = brewer.pal(name = "Set1", n = 5), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines PCA colored by sequencing batch \n", 
                   type, " - ", time))
  print(p) 
  dev.off()
}

# all cell type and time point samples put together ####
for (i in 1:length(mat_lst_bytime)) {
  if (i == 1) {
    mat <- mat_lst_bytime[[i]]
    colnames(mat) <- paste0(colnames(mat), "_", str_split(names(mat_lst_bytime)[i], "_", 2, T)[,1])
  } else {
    mat_to_append <- mat_lst_bytime[[i]]
    colnames(mat_to_append) <- paste0(colnames(mat_to_append), "_", str_split(names(mat_lst_bytime)[i], "_", 2, T)[,1])
    mat <- cbind(mat, mat_to_append)
  }
}

celllines <- str_remove(colnames(mat), "_[0|1|6]hr_[a-zA-Z]+")
for (i in 1:length(celllines)) {
  print(celllines[i])
  if (i == 1) {
    batch = unique(covar_table_final$batch[covar_table_final$cell_line == celllines[i]])
  } else {
    batch = c(batch, unique(covar_table_final$batch[covar_table_final$cell_line == celllines[i]]))
  }
}
batch <- as.factor(batch)

for (i in 1:length(celllines)) {
  print(celllines[i])
  if (i == 1) {
    seqbatch = unique(covar_table_final$seq_batch[covar_table_final$cell_line == celllines[i]])
  } else {
    seqbatch = c(seqbatch, unique(covar_table_final$seq_batch[covar_table_final$cell_line == celllines[i]]))
  }
}
seqbatch <- as.factor(seqbatch)

cpm_mat <- cpm(mat)

# 1. combat by batch
# com_mat <- ComBat(
#   dat = cpm_mat,
#   batch = batch
# )
# pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
# pdf("./PCA_plots/ComBat_by_batch_all_types_and_time_PCA_plot_by_batch.pdf", width = 8, height = 8)

# 2. combat by sequencing batch
com_mat <- ComBat(
  dat = cpm_mat,
  batch = seqbatch
)
pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
pdf("./PCA_plots/ComBat_by_seqbatch_all_types_and_time_PCA_plot_by_batch.pdf", width = 8, height = 8)

# 3. no combat, raw
# pca_res <- PCA(t(cpm_mat), scale.unit = TRUE, ncp = 5, graph = F)
# pdf("./PCA_plots/raw_all_types_and_time_PCA_plot_by_batch.pdf", width = 8, height = 8)

p <- fviz_pca_ind(pca_res,
                  repel = T,
                  habillage = batch,
                  palette = rainbow(20), 
                  show.legend = F,
                  invisible = "quali", font.x = 10) +
  theme_light() +
  ggtitle(paste0("All cell lines x time point x cell type\nPCA colored by co-cultured batch")) 
print(p)
dev.off()

# pdf("./PCA_plots/ComBat_by_batch_all_types_and_time_PCA_plot_by_seq_batch.pdf", width = 8, height = 8)
pdf("./PCA_plots/ComBat_by_seqbatch_all_types_and_time_PCA_plot_by_seq_batch.pdf", width = 8, height = 8)
# pdf("./PCA_plots/raw_all_types_and_time_PCA_plot_by_seq_batch.pdf", width = 8, height = 8)
p <- fviz_pca_ind(pca_res,
                  repel = T,
                  habillage = seqbatch, 
                  palette = brewer.pal(name = "Set1", n = 5), 
                  show.legend = F,
                  invisible = "quali", font.x = 10) +
  theme_light() +
  ggtitle(paste0("All cell lines x time point x cell type\nPCA colored by sequencing batch"))
print(p) 
dev.off()

types <- str_extract(colnames(mat), "[A-Za-z]+$")
types <- factor(types)
# pdf("./PCA_plots/ComBat_by_batch_all_types_and_time_PCA_plot_by_type.pdf", width = 8, height = 8)
pdf("./PCA_plots/ComBat_by_seqbatch_all_types_and_time_PCA_plot_by_type.pdf", width = 8, height = 8)
# pdf("./PCA_plots/raw_all_types_and_time_PCA_plot_by_type.pdf", width = 8, height = 8)
p <- fviz_pca_ind(pca_res,
                  repel = T,
                  habillage = types, 
                  palette = brewer.pal(name = "Set2", n = 3), 
                  show.legend = F,
                  invisible = "quali", font.x = 10) +
  theme_light() +
  ggtitle(paste0("All cell lines x time point x cell type\nPCA colored by cell type"))
print(p) 
dev.off()

times <- str_extract(colnames(mat), "[0|1|6]hr")
times <- factor(times)
# pdf("./PCA_plots/ComBat_by_batch_all_types_and_time_PCA_plot_by_time.pdf", width = 8, height = 8)
pdf("./PCA_plots/ComBat_by_seqbatch_all_types_and_time_PCA_plot_by_time.pdf", width = 8, height = 8)
# pdf("./PCA_plots/raw_all_types_and_time_PCA_plot_by_time.pdf", width = 8, height = 8)
p <- fviz_pca_ind(pca_res,
                  repel = T,
                  habillage = times, 
                  palette = rev(brewer.pal(name = "Set1", n = 3)), 
                  show.legend = F,
                  invisible = "quali", font.x = 10) +
  theme_light() +
  ggtitle(paste0("All cell lines x time point x cell type\nPCA colored by time point"))
print(p) 
dev.off()
