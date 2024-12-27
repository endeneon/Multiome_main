# Chuxuan Li 04/26/2022
# Check the similarity for cell lines in 20 line RNASeq data, with correlation
#heatmap and PCA plot

# init ####
library(Seurat)
library(Signac)
library(edgeR)

library(FactoMineR)
library(factoextra)

library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(reshape2)


load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/20line_codes_and_rdata/demux_20line_integrated_labeled_obj.RData")
RNA_20line <- integrated_labeled
rm(integrated_labeled)

# separate data by cell line and compute cpm ####
lines_20 <- unique(RNA_20line$cell.line.ident)
for (i in 1:length(lines_20)) {
  print(lines_20[i])
  obj <- subset(RNA_20line, cell.line.ident == lines_20[i])
  if (i == 1) {
    cpm_mat <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
    cpm_mat <- cpm_mat[, 1]
  } else {
    to_bind <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
    to_bind <- to_bind[, 1]
    cpm_mat <- rbind(cpm_mat, to_bind)
  }
}
rownames(cpm_mat) <- c(lines_20)


# calculate and plot coordination ####
cor.mat <- cor(t(cpm_mat))
# corrplot(cor.mat, type = "upper", order = "hclust", 
#          tl.col = "black", tl.srt = 45)
melted_cor_mat <- melt(cor.mat)
ggplot(data = melted_cor_mat, 
       aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.55, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(fill = "correlation") +
  scale_fill_gradientn(colours = brewer.pal(5, "Blues"))

# PCA and plotting ####
pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = TRUE)

batch <- as.factor(c("12", "12", "12",
                     "8", "8", "8", "8", "9", "9", "9", "10", "10", "10", "10",
                     "11", "11", "11", "11"))

fviz_pca_ind(pca_res,
             repel = T,
             habillage = batch,
             palette = brewer.pal(name = "Set1", n = 5), 
             show.legend = F,
             invisible = "quali") +
  theme_light() +
  ggtitle("Cell lines PCA colored by cocultured batch")

# remove batch effect and try again ####
library(sva)
com_mat <- ComBat(
  dat = t(cpm_mat),
  batch = batch
)
pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
fviz_pca_ind(pca_res,
             repel = T,
             habillage = batch,
             palette = brewer.pal(name = "Set1", n = 5), 
             show.legend = F,
             invisible = "quali") +
  theme_light() +
  ggtitle("Cell lines PCA colored by batch after batch correction")


# separate data by cell line and time point ####
lines_20 <- unique(RNA_20line$cell.line.ident)
times <- c("0hr", "1hr", "6hr")
for (k in 1:3){
  t = times[k]
  print(t)
  o <- subset(RNA_20line, time.ident == t)
  for (i in 1:length(lines_20)) {
    print(lines_20[i])
    obj <- subset(o, cell.line.ident == lines_20[i])
    if (i == 1) {
      cpm_mat <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
      cpm_mat <- cpm_mat[, 1]
    } else {
      to_bind <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
      to_bind <- to_bind[, 1]
      cpm_mat <- rbind(cpm_mat, to_bind)
    }
  }
  rownames(cpm_mat) <- c(lines_20)
  cor.mat <- cor(t(cpm_mat))
  melted_cor_mat <- melt(cor.mat)
  jpeg(filename = paste0("correlation_heatmap_", t, ".jpeg"), width = 300, height = 300)
  p <- ggplot(data = melted_cor_mat, 
         aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0.55, vjust = 0.5, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(fill = "correlation") +
    scale_fill_gradientn(colours = brewer.pal(5, "Blues")) +
    ggtitle(t)
  print(p)
  dev.off()
  
  
  pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = TRUE)
  
  batch <- as.factor(c("12", "12", "12",
                       "8", "8", "8", "8", "9", "9", "9", "10", "10", "10", "10",
                       "11", "11", "11", "11"))
  jpeg(filename = paste0("PCA_", t, ".jpeg"), width = 400, height = 350)
  p <- fviz_pca_ind(pca_res,
               repel = T,
               habillage = batch,
               palette = brewer.pal(name = "Set1", n = 5), 
               show.legend = F,
               invisible = "quali") +
    theme_light() +
    ggtitle(paste0("Cell lines PCA colored by cocultured batch - ", t))
  print(p)
  dev.off()
}



# separate data by cell line, time point, cell type ####
lines_20 <- unique(RNA_20line$cell.line.ident)
times <- c("0hr", "1hr", "6hr")
for (k in 1:3){
  t = times[k]
  print(t)
  o <- subset(RNA_20line, time.ident == t)
  for (i in 1:length(lines_20)) {
    print(lines_20[i])
    obj <- subset(o, cell.line.ident == lines_20[i])
    if (i == 1) {
      cpm_mat <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
      cpm_mat <- cpm_mat[, 1]
    } else {
      to_bind <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
      to_bind <- to_bind[, 1]
      cpm_mat <- rbind(cpm_mat, to_bind)
    }
  }
  rownames(cpm_mat) <- c(lines_20)
  cor.mat <- cor(t(cpm_mat))
  melted_cor_mat <- melt(cor.mat)
  jpeg(filename = paste0("correlation_heatmap_", t, ".jpeg"), width = 300, height = 300)
  p <- ggplot(data = melted_cor_mat, 
              aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0.55, vjust = 0.5, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(fill = "correlation") +
    scale_fill_gradientn(colours = brewer.pal(5, "Blues")) +
    ggtitle(t)
  print(p)
  dev.off()
  
  
  pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = TRUE)
  
  batch <- as.factor(c("12", "12", "12",
                       "8", "8", "8", "8", "9", "9", "9", "10", "10", "10", "10",
                       "11", "11", "11", "11"))
  jpeg(filename = paste0("PCA_", t, ".jpeg"), width = 400, height = 350)
  p <- fviz_pca_ind(pca_res,
                    repel = T,
                    habillage = batch,
                    palette = brewer.pal(name = "Set1", n = 5), 
                    show.legend = F,
                    invisible = "quali") +
    theme_light() +
    ggtitle(paste0("Cell lines PCA colored by cocultured batch - ", t))
  print(p)
  dev.off()
}
