# Chuxuan Li 04/19/2022
# QC on 18-line ATACseq data after calling new peaks. 
# For each cell line and each gene, compute cpm/tpm, then plot them in the PCA 
#space to check the distance among cell lines
# also plot heatmap of the cpm distribution.

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

library(future)
plan("multisession", workers = 2)
options(future.globals.maxSize= 107374182400)

load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/ATAC_new_obj_with_MACS2_new_peaks_and_RNA_labels.RData")

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# separate data by cell line and compute cpm ####
lines <- unique(ATAC_new$cell.line.ident)
for (i in 1:length(lines)) {
  print(lines[i])
  obj <- subset(ATAC_new, cell.line.ident == lines[i])
  if (i == 1) {
    cpm_mat <- edgeR::cpm(rowSums(obj@assays$peaks@counts), normalized.lib.sizes = T, log = F)
    cpm_mat <- cpm_mat[, 1]
  } else {
    to_bind <- edgeR::cpm(rowSums(obj@assays$peaks@counts), normalized.lib.sizes = T, log = F)
    to_bind <- to_bind[, 1]
    cpm_mat <- rbind(cpm_mat, to_bind)
  }
}
rownames(cpm_mat) <- lines


# calculate and plot coordination ####
cor.mat <- cor(t(cpm_mat))
# corrplot(cor.mat, type = "upper", order = "hclust", 
#          tl.col = "black", tl.srt = 45)
melted_cor_mat <- melt(cor.mat)
ggplot(data = melted_cor_mat, 
       aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.55, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(fill = "correlation") +
  scale_fill_gradientn(colours = brewer.pal(5, "Blues"))

# PCA and plotting ####
pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = TRUE)
batch <- as.factor(c("10", "9", "12", "9", "8", "9", "12", "11", "10", "8", 
                   "10", "11", "11", "10", "8", "11", "12", "8"))
fviz_pca_ind(pca_res,
             repel = T,
             habillage = batch,
             palette = brewer.pal(name = "Set1", n = 5), 
             show.legend = F,
             invisible = "quali") +
  theme_bw()


# separate data by cell line, time point, celltype, then compute cpm ####
ATAC_new <- subset(ATAC_new, cell.type != "unknown")
lines <- unique(ATAC_new$cell.line.ident)
times <- c("0hr", "1hr", "6hr")
types <- as.vector(unique(ATAC_new$cell.type))
samples <- rep_len("", length(lines) * length(types))
types_rep <- as.factor(rep_len(types, length(lines) * length(types)))

for (k in 1:3){
  t = times[k]
  print(t)
  timed <- subset(ATAC_new, time.ident == t)
  for (j in 1:length(types)) {
    print(types[j])
    o <- subset(timed, cell.type == types[j])
    for (i in 1:length(lines)) {
      print(lines[i])
      obj <- subset(o, cell.line.ident == lines[i])
      if (i == 1) {
        cpm_mat <- edgeR::cpm(rowSums(obj@assays$peaks@counts), normalized.lib.sizes = T, log = F)
        cpm_mat <- cpm_mat[, 1]
      } else {
        to_bind <- edgeR::cpm(rowSums(obj@assays$peaks@counts), normalized.lib.sizes = T, log = F)
        to_bind <- to_bind[, 1]
        cpm_mat <- rbind(cpm_mat, to_bind)
      }
    }
    rownames(cpm_mat) <- lines
    cor.mat <- cor(t(cpm_mat))
    # order by group, not hclust
    batch <- c("10", "9", "12", "9", "8", "9", "12", "11", "10", "8", 
               "10", "11", "11", "10", "8", "11", "12", "8")
    batch <- as.numeric(batch)
    names(batch) <- lines
    sort(batch)
    order <- names(sort((batch)))
    cor.mat <- cor.mat[order, order]
    melted_cor_mat <- melt(cor.mat)
    pdf(file = paste0("correlation_heatmap_by_line_and_type_", types[j], "_", t, ".pdf"))
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
    
    pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = TRUE)
    batch <- as.factor(c("10", "9", "12", "9", "8", "9", "12", "11", "10", "8", 
                         "10", "11", "11", "10", "8", "11", "12", "8"))
    pdf(file = paste0("PCA_by_line_and_type", types[j], "_", t, ".pdf"))
    p <- fviz_pca_ind(pca_res,
                      repel = T,
                      habillage = batch,
                      palette = brewer.pal(name = "Set1", n = 5), 
                      show.legend = F,
                      invisible = "quali") +
      theme_light() +
      ggtitle(paste0("Cell lines PCA colored by cocultured batch - ", types[j], " ", t))
    print(p)
    dev.off()
  }
  
}

# compute top diff accessible peaks for each time point ####
DefaultAssay(ATAC_new) <- "peaks"
Idents(ATAC_new) <- "time.ident"
by_time_markers <- FindAllMarkers(ATAC_new, assay = "peaks", 
                                  test.use = "MAST", random.seed = 10, 
                                  logfc.threshold = 0.05)

# use the diff accessible peaks for correlation matrix ####
lines <- unique(ATAC_new$cell.line.ident)
counts <- GetAssayData(ATAC_new)
counts <- counts[by_time_markers$gene, ]
subATAC <- subset(ATAC_new, features = rownames(counts))
for (i in 1:length(lines)) {
  print(lines[i])
  obj <- subset(subATAC, cell.line.ident == lines[i])
  if (i == 1) {
    cpm_mat <- edgeR::cpm(rowSums(obj@assays$peaks@counts), normalized.lib.sizes = T, log = F)
    cpm_mat <- cpm_mat[, 1]
  } else {
    to_bind <- edgeR::cpm(rowSums(obj@assays$peaks@counts), normalized.lib.sizes = T, log = F)
    to_bind <- to_bind[, 1]
    cpm_mat <- rbind(cpm_mat, to_bind)
  }
}
rownames(cpm_mat) <- lines


# calculate and plot coordination ####
cor.mat <- cor(t(cpm_mat))
# corrplot(cor.mat, type = "upper", order = "hclust", 
#          tl.col = "black", tl.srt = 45)

cor.mat <- reorder_cormat(cor.mat)
melted_cor_mat <- melt(cor.mat)

ggplot(data = melted_cor_mat, 
       aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.2, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(fill = "correlation") +
  scale_fill_gradientn(colours = brewer.pal(6, "Reds"))

# order by group, not hclust
batch <- c("10", "9", "12", "9", "8", "9", "12", "11", "10", "8", 
                     "10", "11", "11", "10", "8", "11", "12", "8")
batch <- as.numeric(batch)
names(batch) <- lines
sort(batch)
order <- names(sort((batch)))
order
cor.mat <- cor.mat[order, order]
melted_cor_mat <- melt(cor.mat)
ggplot(data = melted_cor_mat, 
       aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.58, vjust = 0.6, size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(fill = "correlation") +
  scale_fill_gradientn(colours = brewer.pal(6, "Reds"))


# plot peak tracks for BDNF for each line to check intensity ####
lines <- unique(multiomic_obj$cell.line.ident)

for (i in lines) {
  subobj <- subset(multiomic_obj, cell.line.ident == i)
  Annotation(subobj) <- annotation_ranges
  filename <- paste0("./peak_tract_plots_by_cell_line/", "BDNF_", i, 
                     ".pdf")
  print(filename)
  
  pdf(file = filename, height = 4, width = 8)
  p <- CoveragePlot(
    object = subobj, 
    assay = "peaks", expression.assay = "SCT", 
    group.by = "time.ident", sep = "-", annotation = T, peaks = T,
    region = "BDNF", extend.upstream = 100000, extend.downstream = 100000, 
    region.highlight = GRanges(seqnames = "chr11", strand = "+",
                               ranges = IRanges(start = 27763500, end = 27777500))
  )
  print(p)
  dev.off()
}

