# Chuxuan Li 05/09/2022
# compare the ATACseq data of 5-line and 18-line experiments by PCA plot
#for each of the time points and cell types, and correlation heatmap for each
#cell type at 0hr

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

load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/ATAC_new_obj_with_MACS2_new_peaks_and_RNA_labels.RData")

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
