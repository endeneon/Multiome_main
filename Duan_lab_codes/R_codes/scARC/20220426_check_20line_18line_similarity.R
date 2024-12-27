# Chuxuan Li 04/26/2022
# Check the similarity between 20-line and 18-line RNASeq data, with correlation
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
library(ggrepel)
library(reshape2)


load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/20line_codes_and_rdata/demux_20line_integrated_labeled_obj.RData")
RNA_20line <- integrated_labeled
load("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")
RNA_18line <- integrated_labeled
rm(integrated_labeled)

# prepare for analysis ####
genes.use <- rownames(RNA_18line@assays$RNA)
counts <- GetAssayData(RNA_20line, assay = "RNA")
counts <- counts[which(rownames(counts) %in% genes.use), ]
RNA_20line <- subset(RNA_20line, features = rownames(counts))

# auxiliary functions
# separate data by cell line and compute cpm
makeCPMMat <- function(obj1, obj2, lines_1, lines_2){
  j = 1
  for (i in 1:(length(lines_1) + length(lines_2))) {
    if (i <= length(lines_1)) {
      print(lines_1[i])
      obj <- subset(obj1, cell.line.ident == lines_1[i])
    } else {
      print(lines_2[j])
      obj <- subset(obj2, cell.line.ident == lines_2[j])
      j = j + 1
    }
    if (i == 1) {
      cpm_mat <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
      cpm_mat <- cpm_mat[, 1]
    } else {
      to_bind <- edgeR::cpm(rowSums(obj@assays$RNA@counts), normalized.lib.sizes = T, log = F)
      to_bind <- to_bind[, 1]
      cpm_mat <- rbind(cpm_mat, to_bind)
    }
  }
  rownames(cpm_mat) <- c(lines_1, lines_2)
  return(cpm_mat)
}

# calculate and plot correlation
plotCorHeatmap <- function(cpm_mat, title){
  cor.mat <- cor(t(cpm_mat))
  melted_cor_mat <- melt(cor.mat)
  p <- ggplot(data = melted_cor_mat, 
         aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0.55, vjust = 0.5, size = 6),
          axis.text.y = element_text(angle = 45, hjust = 0.55, vjust = 0.5, size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(fill = "correlation") +
    scale_fill_gradientn(colours = brewer.pal(5, "Blues")) +
    ggtitle(title)
  return(p)
}


lines_20 <- unique(RNA_20line$cell.line.ident)
lines_18 <- unique(RNA_18line$cell.line.ident)
batch <- as.factor(c("3", "3", "3", "3", "4", "4", "4", "4", "5", "5", "5", "5",
                     "6", "6", "6", "6", "7", "7", "7", "7", "12", "12", "12",
                     "8", "8", "8", "8", "9", "9", "9", "10", "10", "10", "10",
                     "11", "11", "11", "11"))
exp <- as.factor(c(rep_len("20 line", 20),
                   rep_len("18 line", 18)))
RNA_18line$cell.type[RNA_18line$cell.type %in% 
                       c("SST_pos_GABA", "GABA", "SEMA3E_pos_GABA")] <-
  "GABA"
RNA_18line <- subset(RNA_18line, cell.type != "unknown")
types <- unique(RNA_18line$cell.type)
types
unique(RNA_20line$cell.type)
RNA_20line$cell.type[RNA_20line$cell.type %in% 
                       c("GABA", "SEMA3E_pos_GABA")] <- "GABA"
RNA_20line <- subset(RNA_20line, cell.type != "unknown")
times <- unique(RNA_18line$time.ident)
times


# PCA for cell type, time point separated ####
for (i in types) {
  print(i)
  subtype_18 <- subset(RNA_18line, cell.type == i)
  subtype_20 <- subset(RNA_20line, cell.type == i)
  for (j in times) {
    subtime18 <- subset(subtype_18, time.ident == j)
    subtime20 <- subset(subtype_20, time.ident == j)
    cpm_mat <- makeCPMMat(subtime20, subtime18, lines_20, lines_18)
    pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = F)
    filename <- paste0(i, "_", j, "_PCA_plot_by_batch.pdf")
    pdf(filename)
    p <- fviz_pca_ind(pca_res,
                 repel = T,
                 habillage = batch,
                 palette = brewer.pal(name = "Paired", n = 10), 
                 show.legend = F,
                 invisible = "quali") +
      theme_light() +
      ggtitle(paste0("Cell lines PCA colored by cocultured batch \n", 
                     i, " - ", j))
    print(p)
    dev.off()
   
    filename <- paste0(i, "_", j, "_PCA_plot_by_exp.pdf")
    pdf(filename)
    p <- fviz_pca_ind(pca_res,
                 repel = T,
                 habillage = exp,
                 palette = brewer.pal(name = "Set1", n = 2), 
                 show.legend = F,
                 invisible = "quali") +
      theme_light() +
      ggtitle(paste0("Cell lines PCA colored by 20/18 line \n", 
                     i, " - ", j))
    print(p) 
    dev.off()
  }
}

# correlation heatmap for cell type, time point separated ####
for (i in types) {
  print(i)
  subtype_18 <- subset(RNA_18line, cell.type == i)
  subtype_20 <- subset(RNA_20line, cell.type == i)
  for (j in times) {
    subtime18 <- subset(subtype_18, time.ident == j)
    subtime20 <- subset(subtype_20, time.ident == j)
    cpm_mat <- makeCPMMat(subtime18, subtime20)
    title <- paste0("correlation heatmap \n", 
                    i, " - ", j)
    cor_heatmap <- plotCorHeatmap(cpm_mat, title = title)
    filename <- paste0(i, "_", j, "_corr_heatmap.pdf")
    pdf(file = filename)
    print(cor_heatmap)
    dev.off()
  }
}

# remove batch effect and try again ####
library(sva)

for (i in types) {
  print(i)
  subtype_18 <- subset(RNA_18line, cell.type == i)
  subtype_20 <- subset(RNA_20line, cell.type == i)
  for (j in times) {
    subtime18 <- subset(subtype_18, time.ident == j)
    subtime20 <- subset(subtype_20, time.ident == j)
    cpm_mat <- makeCPMMat(subtime20, subtime18, lines_20, lines_18)
    com_mat <- ComBat(
      dat = t(cpm_mat),
      batch = batch
      #batch = exp
    )
    pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
    filename <- paste0("ComBat_by_batch_", i, "_", j, "_PCA_plot_by_batch.pdf")
    pdf(filename)
    p <- fviz_pca_ind(pca_res,
                      repel = T,
                      habillage = batch,
                      palette = brewer.pal(name = "Paired", n = 10), 
                      show.legend = F,
                      invisible = "quali") +
      theme_light() +
      ggtitle(paste0("Cell lines PCA colored by cocultured batch \n", 
                     i, " - ", j))
    print(p)
    dev.off()
    
    filename <- paste0("ComBat_by_batch_", i, "_", j, "_PCA_plot_by_exp.pdf")
    pdf(filename)
    p <- fviz_pca_ind(pca_res,
                      repel = T,
                      habillage = exp,
                      palette = brewer.pal(name = "Set1", n = 2), 
                      show.legend = F,
                      invisible = "quali") +
      theme_light() +
      ggtitle(paste0("Cell lines PCA colored by 20/18 line \n", 
                     i, " - ", j))
    print(p) 
    dev.off()
  }
}


# put cell types into the same PCA plot ####
makeCpmMatCelltypes <- function(obj1, obj2, lines_1, lines_2) {
  rownames <- rep_len(NA, (length(lines_1) + length(lines_2)) * length(types))
  j = 1
  k = 1
  for (i in 1:(length(lines_1) + length(lines_2))) {
    if (i <= length(lines_1)) {
      print(lines_1[i])
      rownames[k] <- lines_1[i]
      rownames[k + 1] <- lines_1[i]
      rownames[k + 2] <- lines_1[i]
      obj <- subset(obj1, cell.line.ident == lines_1[i])
    } else {
      print(lines_2[j])
      rownames[k] <- lines_2[j]
      rownames[k + 1] <- lines_2[j]
      rownames[k + 2] <- lines_2[j]
      obj <- subset(obj2, cell.line.ident == lines_2[j])
      j = j + 1
    }
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
    }
  }
  rownames(cpm_mat) <- rownames
  return(cpm_mat)
}
batch <- as.factor(c("3", "3", "3", "3", "4", "4", "4", "4", "5", "5", "5", "5",
                     "6", "6", "6", "6", "7", "7", "7", "7", "12", "12", "12",
                     "8", "8", "8", "8", "9", "9", "9", "10", "10", "10", "10",
                     "11", "11", "11", "11"))
batch_more <- rep(batch, each = 3)
exp <- as.factor(c(rep_len("20 line", 60),
                   rep_len("18 line", 54)))
for (time in times) {
  subtime18 <- subset(RNA_18line, time.ident == time)
  subtime20 <- subset(RNA_20line, time.ident == time)
  cpm_mat <- makeCpmMatCelltypes(subtime20, subtime18, lines_20, lines_18)
  #pca_res <- PCA(cpm_mat, scale.unit = TRUE, ncp = 5, graph = F)
  com_mat <- ComBat(
    dat = t(cpm_mat),
    #batch = batch_more
    batch = exp
  )
  pca_res <- PCA(t(com_mat), scale.unit = TRUE, ncp = 5, graph = TRUE)
  filename <- paste0("ComBat_by_exp_", 
                     #"ComBat_by_batch_", 
    time, "_PCA_plot_by_batch.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = F,
                    habillage = batch_more,
                    palette = brewer.pal(name = "Paired", n = 10), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x cell type PCA colored by cocultured batch \n", 
                   time)) 
  print(p)
  dev.off()
  
  filename <- paste0("ComBat_by_exp_", 
    #"ComBat_by_batch_", 
    time, "_PCA_plot_by_exp.pdf")
  pdf(filename, width = 8, height = 8)
  p <- fviz_pca_ind(pca_res,
                    repel = F,
                    habillage = exp, 
                    palette = brewer.pal(name = "Set1", n = 2), 
                    show.legend = F,
                    invisible = "quali", font.x = 10) +
    theme_light() +
    ggtitle(paste0("Cell lines x cell type PCA colored by 20/18 line \n", time))
  print(p) 
  dev.off()
}
