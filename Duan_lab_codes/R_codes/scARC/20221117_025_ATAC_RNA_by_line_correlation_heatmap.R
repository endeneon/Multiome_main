# Chuxuan Li 11/17/2022
# plot correlation heatmap of the cpm matrices of 025 RNA and ATAC data.

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

load("~/NVME/scARC_Duan_018/Duan_project_025_ATAC/ATAC_clean_obj_added_MACS2_peaks.RData")

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# separate data by cell line, time point, celltype, then compute cpm ####
unique(ATAC_clean$RNA.cell.type)
ATAC_clean <- subset(ATAC_clean, RNA.cell.type != "unidentified")
lines <- sort(unique(ATAC_clean$cell.line.ident))
times <- c("0hr", "1hr", "6hr")
types <- as.vector(sort(unique(ATAC_clean$RNA.cell.type)))
samples <- rep_len("", length(lines) * length(types))
types_rep <- as.factor(rep_len(types, length(lines) * length(types)))
batch <- c()
for (l in lines) {
  print(l)
  batch <- c(batch, 
             unique(str_split(ATAC_clean$group.ident[ATAC_clean$cell.line.ident == l],
                              pattern = "-", n = 2, simplify = T)[1]))
}
names(batch) <- lines
batch_factor <- as.factor(batch)

for (k in 1:3){
  t = times[k]
  print(t)
  timed <- subset(ATAC_clean, time.ident == t)
  for (j in 1:length(types)) {
    print(types[j])
    o <- subset(timed, RNA.cell.type == types[j])
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
    
    pdf(file = paste0("PCA_by_line_and_type", types[j], "_", t, ".pdf"))
    p <- fviz_pca_ind(pca_res,
                      repel = T,
                      habillage = batch_factor,
                      palette = brewer.pal(name = "Set1", n = 5), 
                      show.legend = F,
                      invisible = "quali") +
      theme_light() +
      ggtitle(paste0("Cell lines PCA colored by cocultured batch - ", types[j], " ", t))
    print(p)
    dev.off()
  }
  
}
