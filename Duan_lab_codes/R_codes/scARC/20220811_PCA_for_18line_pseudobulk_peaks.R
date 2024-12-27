
# Chuxuan Li 08/11/2022
# Examine PCA plots for ATACseq pseudobulk peaks, 

# init ####
library(Signac)
library(ggplot2)
library(sva)
library(FactoMineR)
library(factoextra)
library(limma)

library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(reshape2)

load("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/da_peaks_limma/post-combat_count_matrices.RData")

# auxillary functions ####
createDGE <- function(count_matrix){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = celllines)
  A <- rowSums(y$counts)
  #isexpr <- A > 10
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}
filterByCpm <- function(y, cutoff, nsample) {
  cpm <- cpm(y)
  colInd0 <- seq(1, ncol(y), by = 3)
  colInd1 <- seq(2, ncol(y), by = 3)
  colInd6 <- seq(3, ncol(y), by = 3)
  
  cpm_0hr <- cpm[, colInd0]
  cpm_1hr <- cpm[, colInd1]
  cpm_6hr <- cpm[, colInd6]
  #passfilter <- rowSums(zerohr_cpm >= cutoff) >= nsample
  passfilter <- (rowSums(cpm_0hr >= cutoff) >= nsample |
                   rowSums(cpm_1hr >= cutoff) >= nsample |
                   rowSums(cpm_6hr >= cutoff) >= nsample)
  return(passfilter)
}

setwd("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/da_peaks_limma")

# 1. by time- use unfiltered matrix ####
options(ggrepel.max.overlaps = 25)
y_GABA <- cpm(createDGE(GABA_mat_adj))
y_nmglut <- cpm(createDGE(nmglut_mat_adj))
y_npglut <- cpm(createDGE(npglut_mat_adj))
batch = as.factor(rep(c("12", "11", "12", "12", "8", "8", "9", "9", "10", "8", "8", "9", 
              "10", "10", "11", "11", "11", "10"), each = 3))
batch
types <- c("GABA", "nmglut", "npglut")
mat_list <- list(y_GABA, y_nmglut, y_npglut)
for (i in 1:length(types)) {
  print(types[i])
  mat <- mat_list[[i]]
  pca_res <- PCA(t(mat), scale.unit = TRUE, ncp = 5, graph = TRUE) # combat 
  #required transposed matrix, here need to reverse that
  
  filename <- paste0(types[i], "_PCA_plot_by_line_and_time.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = T,
                    habillage = batch, 
                    palette = brewer.pal(name = "Paired", n = 10), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x time point PCA colored by cocultured batch \n", 
                   types[i])) 
  print(p)
  dev.off()
}

# 2. by time - use cpm filtered matrix ####
for (i in 1:length(types)) {
  print(types[i])
  mat <- mat_list[[i]]
  ind.keep <- filterByCpm(mat, 1, 9)
  mat <- mat[ind.keep, ]
  pca_res <- PCA(t(mat), scale.unit = TRUE, ncp = 5, graph = TRUE) # combat 
  #required transposed matrix, here need to reverse that
  
  filename <- paste0(types[i], "_filtered_PCA_plot_by_line_and_time.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = T,
                    habillage = batch, 
                    palette = brewer.pal(name = "Paired", n = 10), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x time point PCA colored by cocultured batch \n", 
                   "filtered by cpm - ", types[i])) 
  print(p)
  dev.off()
}

# 3. by type - use unfiltered matrix ####
load("~/Data/FASTQ/Duan_Project_024/hybrid_output/R_DE_peaks_pseudobulk/frag_count_by_cell_type_specific_peaks_22Jul2022.RData")

GABA <- Seurat_object_list[[2]]
nmglut <- Seurat_object_list[[3]]
npglut <- Seurat_object_list[[4]]
rm(Seurat_object_list)
cellines <- sort(unique(nmglut$cell.line.ident))
types <- sort(unique(nmglut$time.ident))
nmglut$typexline.ident <- NA
for (l in cellines) {
  for (t in times) {
    id <- paste(l, t, sep = "-")
    print(id)
    nmglut$timexline.ident[nmglut$cell.line.ident == l &
                             nmglut$time.ident == t] <- id
  }
}
unique(nmglut$timexline.ident)
lts <- sort(unique(nmglut$timexline.ident))

makeMat4CombatseqType <- function(typeobj, linetimes) {
  rownames <- rep_len(NA, (length(linetimes)))
  for (i in 1:length(linetimes)) {
    print(linetimes[i])
    rownames[i] <- linetimes[i]
    obj <- subset(typeobj, timexline.ident == linetimes[i])
    count_matrix <- obj 
    if (i == 1) {
      mat <- as.array(rowSums(obj@assays$ATAC@counts))
    } else {
      to_bind <- as.array(rowSums(obj@assays$ATAC@counts))
      mat <- rbind(mat, to_bind)
    }
  }
  rownames(mat) <- rownames
  return(mat)
}

GABA_mat <- t(makeMat4CombatseqLine(GABA, lts))
nmglut_mat <- t(makeMat4CombatseqLine(nmglut, lts))
npglut_mat <- t(makeMat4CombatseqLine(npglut, lts))

# combat-seq
batch = rep(c("12", "11", "12", "12", "8", "8", "9", "9", "10", "8", "8", "9", 
              "10", "10", "11", "11", "11", "10"), each = 3)
batch
length(batch) #54
ncol(GABA_mat) #54

GABA_mat_adj <- ComBat_seq(GABA_mat, batch = batch, 
                           group = NULL)
nmglut_mat_adj <- ComBat_seq(nmglut_mat, batch = batch,
                             group = NULL)
npglut_mat_adj <- ComBat_seq(npglut_mat, batch = batch,
                             group = NULL)