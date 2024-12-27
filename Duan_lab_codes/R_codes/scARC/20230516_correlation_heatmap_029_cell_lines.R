# Chuxuan Li 05/16/2023
# plot correlation heatmap of the cpm matrices of 029 RNA data.

# init ####
library(Seurat)
library(Signac)
library(edgeR)

library(FactoMineR)
library(factoextra)
library(sva)

library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(stringr)
library(future)
plan("multisession", workers = 2)
options(future.globals.maxSize= 107374182400)

load("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/029_RNA_integrated_labeled.RData")
integrated_labeled_20line <- integrated_labeled
load("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/integrated_labeled_40201.RData")
integrated_labeled_4line <- integrated_labeled
rm(integrated_labeled)

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# separate data by cell line, time point, celltype, then compute cpm ####
unique(integrated_labeled_20line$cell.type)
integrated_labeled_20line <- subset(integrated_labeled_20line, cell.type != "unidentified")
integrated_labeled_4line <- subset(integrated_labeled_4line, cell.type != "unidentified")
lines20 <- sort(unique(integrated_labeled_20line$cell.line.ident))
lines4 <- sort(unique(integrated_labeled_4line$cell.line.ident))
times <- c("0hr", "1hr", "6hr")
types <- as.vector(sort(unique(integrated_labeled_20line$cell.type)))
samples <- rep_len("", (length(lines20) + length(lines4)) * length(types))
types_rep <- as.factor(rep_len(types, (length(lines20) + length(lines4)) * length(types)))
batch <- c()
for (l in lines20) {
  print(l)
  batch <- c(batch, 
             unique(str_split(integrated_labeled_20line$orig.ident[integrated_labeled_20line$cell.line.ident == l],
                              pattern = "-", n = 2, simplify = T)[1]))
}
names(batch) <- lines20
for (l in lines4) {
  print(l)
  batch <- c(batch, 
             unique(str_split(integrated_labeled_4line$orig.ident[integrated_labeled_4line$cell.line.ident == l],
                              pattern = "-", n = 2, simplify = T)[1]))
}
names(batch)[21:length(batch)] <- lines4
batch_factor <- as.factor(batch)
types <- types[types != "GABA"]
for (k in 1:3){
  t = times[k]
  print(t)
  timed20 <- subset(integrated_labeled_20line, time.ident == t)
  timed4 <- subset(integrated_labeled_4line, time.ident == t)
  for (j in 1:length(types)) {
    print(types[j])
    o <- subset(timed20, cell.type == types[j])
    for (i in 1:length(lines20)) {
      print(lines20[i])
      obj <- subset(o, cell.line.ident == lines20[i])
      if (i == 1) {
        cpm_mat <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
        cpm_mat <- cpm_mat[, 1]
      } else {
        to_bind <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
        to_bind <- to_bind[, 1]
        cpm_mat <- cbind(cpm_mat, to_bind)
      }
    }
    colnames(cpm_mat) <- lines20
    
    o <- subset(timed4, cell.type == types[j])
    for (i in 1:length(lines4)) {
      print(lines4[i])
      obj <- subset(o, cell.line.ident == lines4[i])
      to_bind <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
      to_bind <- to_bind[, 1]
      cpm_mat <- cbind(cpm_mat, to_bind)
    }
    colnames(cpm_mat)[21:24] <- lines4
    
    cor.mat <- cor(cpm_mat)
    # order by group, not hclust
    order <- names(sort((batch)))
    cor.mat <- cor.mat[order, order]
    melted_cor_mat <- melt(cor.mat)
    pdf(file = paste0("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/correlation_heatmaps_with_40201/correlation_heatmap_by_line_and_type_", types[j], "_", t, ".pdf"))
    p <- ggplot(data = melted_cor_mat, 
                aes(x = Var1, y = Var2, fill = value)) + 
      geom_tile() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 0.55, vjust = 0.5, size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      labs(fill = "correlation") +
      scale_fill_gradientn(colours = brewer.pal(5, "Blues")) +
      ggtitle(paste0("Correlation heatmap - ", types[j], " ", t))
    print(p)
    dev.off()
  }
}

# combat-seq batch correction ####
for (k in 1:3){
  t = times[k]
  print(t)
  timed20 <- subset(integrated_labeled_20line, time.ident == t)
  timed4 <- subset(integrated_labeled_4line, time.ident == t)
  for (j in 1:length(types)) {
    print(types[j])
    o <- subset(timed20, cell.type == types[j])
    for (i in 1:length(lines20)) {
      print(lines20[i])
      obj <- subset(o, cell.line.ident == lines20[i])
      if (i == 1) {
        cpm_mat <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
        cpm_mat <- cpm_mat[, 1]
      } else {
        to_bind <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
        to_bind <- to_bind[, 1]
        cpm_mat <- cbind(cpm_mat, to_bind)
      }
    }
    colnames(cpm_mat) <- lines20
    
    o <- subset(timed4, cell.type == types[j])
    for (i in 1:length(lines4)) {
      print(lines4[i])
      obj <- subset(o, cell.line.ident == lines4[i])
      to_bind <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
      to_bind <- to_bind[, 1]
      cpm_mat <- cbind(cpm_mat, to_bind)
    }
    colnames(cpm_mat)[21:24] <- lines4
    com_mat <- ComBat(
      dat = cpm_mat,
      batch = batch_factor
    )
    cor.mat <- cor(com_mat)
    # order by group, not hclust
    order <- names(sort((batch)))
    cor.mat <- cor.mat[order, order]
    melted_cor_mat <- melt(cor.mat)
    pdf(file = paste0("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/correlation_heatmaps_with_40201/after_combat/correlation_heatmap_combat__corrected_", types[j], "_", t, ".pdf"))
    p <- ggplot(data = melted_cor_mat, 
                aes(x = Var1, y = Var2, fill = value)) + 
      geom_tile() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 0.55, vjust = 0.5, size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
      labs(fill = "correlation") +
      scale_fill_gradientn(colours = brewer.pal(5, "Blues")) +
      ggtitle(paste0("Correlation heatmap (after batch correction) - ", types[j], " ", t))
    print(p)
    dev.off()
  }
}

# do GABA separately ####
lines4 <- lines4[2:4]
batch <- batch[c(1:20, 22:24)]
batch_factor <- as.factor(batch)
for (k in 1:3){
  t = times[k]
  print(t)
  timed20 <- subset(integrated_labeled_20line, time.ident == t)
  timed4 <- subset(integrated_labeled_4line, time.ident == t)
  o <- subset(timed20, cell.type == "GABA")
  for (i in 1:length(lines20)) {
    print(lines20[i])
    obj <- subset(o, cell.line.ident == lines20[i])
    if (i == 1) {
      cpm_mat <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
      cpm_mat <- cpm_mat[, 1]
    } else {
      to_bind <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
      to_bind <- to_bind[, 1]
      cpm_mat <- cbind(cpm_mat, to_bind)
    }
  }
  colnames(cpm_mat) <- lines20
    
  o <- subset(timed4, cell.type == "GABA")
  for (i in 1:length(lines4)) {
    print(lines4[i])
    obj <- subset(o, cell.line.ident == lines4[i])
    to_bind <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
    to_bind <- to_bind[, 1]
    cpm_mat <- cbind(cpm_mat, to_bind)
  }
  colnames(cpm_mat)[21:23] <- lines4
  cor.mat <- cor(cpm_mat)
  # order by group, not hclust
  order <- names(sort((batch)))
  cor.mat <- cor.mat[order, order]
  melted_cor_mat <- melt(cor.mat)
  pdf(file = paste0("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/correlation_heatmaps_with_40201/correlation_heatmap_by_line_and_type_GABA_", t, ".pdf"))
  p <- ggplot(data = melted_cor_mat, 
              aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0.55, vjust = 0.5, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(fill = "correlation") +
    scale_fill_gradientn(colours = brewer.pal(5, "Blues")) +
    ggtitle(paste0("Correlation heatmap - ", "GABA", " ", t))
  print(p)
  dev.off()
}

for (k in 1:3){
  t = times[k]
  print(t)
  timed20 <- subset(integrated_labeled_20line, time.ident == t)
  timed4 <- subset(integrated_labeled_4line, time.ident == t)
  o <- subset(timed20, cell.type == "GABA")
  for (i in 1:length(lines20)) {
    print(lines20[i])
    obj <- subset(o, cell.line.ident == lines20[i])
    if (i == 1) {
      cpm_mat <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
      cpm_mat <- cpm_mat[, 1]
    } else {
      to_bind <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
      to_bind <- to_bind[, 1]
      cpm_mat <- cbind(cpm_mat, to_bind)
    }
  }
  colnames(cpm_mat) <- lines20
  
  o <- subset(timed4, cell.type == "GABA")
  for (i in 1:length(lines4)) {
    print(lines4[i])
    obj <- subset(o, cell.line.ident == lines4[i])
    to_bind <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
    to_bind <- to_bind[, 1]
    cpm_mat <- cbind(cpm_mat, to_bind)
  }
  colnames(cpm_mat)[21:23] <- lines4
  com_mat <- ComBat(
    dat = cpm_mat,
    batch = batch_factor
  )
  cor.mat <- cor(com_mat)
  # order by group, not hclust
  order <- names(sort((batch)))
  cor.mat <- cor.mat[order, order]
  melted_cor_mat <- melt(cor.mat)
  pdf(file = paste0("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/correlation_heatmaps_with_40201/after_combat/correlation_heatmap_combat__corrected_", "GABA", "_", t, ".pdf"))
  p <- ggplot(data = melted_cor_mat, 
              aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0.55, vjust = 0.5, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(fill = "correlation") +
    scale_fill_gradientn(colours = brewer.pal(5, "Blues")) +
    ggtitle(paste0("Correlation heatmap (after batch correction) - GABA ", t))
  print(p)
  dev.off()
}

