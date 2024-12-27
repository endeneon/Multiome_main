# Chuxuan Li 06/20/2022
# use 5-line RNA and ATAC data, and pseudobulk-generated DEG list, plot heatmaps
#of gene activity score and gene expression for the top 1000 DEGs using log
#normalized counts exponentiated and then averaged (Seurat method)

# init ####
library(Seurat)
library(Signac)
library(GenomicRanges)
library(gplots)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(future)
library(colorRamps)

plan("multisession", workers = 1)

# load data ####
load("./5line_multiomic_obj_complete_added_geneactivities.RData")
DEG_paths <- list.files("../Analysis_RNAseq_pseudobulk/5line_DE_results/filtered_by_basemean/significant/gene_name_only/", 
           full.names = T, recursive = F)
DEG_paths <- DEG_paths[str_detect(DEG_paths, "NPC", T)]
DEG_paths_nmglut_up <- DEG_paths[str_detect(DEG_paths, "nmglut", F)]
DEG_paths_nmglut_up <- DEG_paths_nmglut_up[str_detect(DEG_paths_nmglut_up, "up", F)]
DEG_paths_npglut_up <- DEG_paths[str_detect(DEG_paths, "npglut", F)]
DEG_paths_npglut_up <- DEG_paths_npglut_up[str_detect(DEG_paths_npglut_up, "up", F)]
DEG_paths_GABA_up <- DEG_paths[str_detect(DEG_paths, "GABA", F)]
DEG_paths_GABA_up <- DEG_paths_GABA_up[str_detect(DEG_paths_GABA_up, "up", F)]
#DEG_lists <- vector("list", length(DEG_paths_nmglut_up))
for (i in 1:length(DEG_lists)) {
  DEG_paths_nmglut_up[[i]] <- read_delim(DEG_paths_nmglut_up[i], delim = "\t", col_names = F)
  DEG_paths_npglut_up[[i]] <- read_delim(DEG_paths_npglut_up[i], delim = "\t", col_names = F)
  DEG_paths_GABA_up[[i]] <- read_delim(DEG_paths_GABA_up[i], delim = "\t", col_names = F)
}
genelist <- intersect(DEG_paths_nmglut_up[[1]]$X1, DEG_paths_nmglut_up[[2]]$X1)
genelist <- intersect(genelist, DEG_paths_npglut_up[[1]]$X1)
genelist <- intersect(genelist, DEG_paths_npglut_up[[2]]$X1)
genelist <- intersect(genelist, DEG_paths_GABA_up[[1]]$X1)
genelist <- intersect(genelist, DEG_paths_GABA_up[[2]]$X1)
#genelist <- genelist[1:1000]

# compute average gex and gact matrices ####
unique(obj_complete$fine.cell.type)
nmglut <- subset(obj_complete, fine.cell.type == "NEFM-/CUX2+ glut")
npglut <- subset(obj_complete, fine.cell.type == "NEFM+/CUX2- glut")
GABA <- subset(obj_complete, fine.cell.type == "GABA")

times <- sort(unique(nmglut$time.ident))
gex_df <- data.frame(GABA_zero = genelist, GABA_one = genelist, GABA_six = genelist,
                     nmglut_zero = genelist, nmglut_one = genelist, nmglut_six = genelist,
                     npglut_zero = genelist, npglut_one = genelist, npglut_six = genelist)
gact_df <- data.frame(GABA_zero = genelist, GABA_one = genelist, GABA_six = genelist,
                      nmglut_zero = genelist, nmglut_one = genelist, nmglut_six = genelist,
                      npglut_zero = genelist, npglut_one = genelist, npglut_six = genelist)
subobjlst <- list(GABA, nmglut, npglut)
for (k in 1:length(subobjlst)) {
  ctobj <- subobjlst[[k]]
  for (j in 1:length(times)) {
    obj <- subset(ctobj, time.ident == times[j])
    for (i in 1:length(genelist)) {
      g <- genelist[i]
      if (i == 1) {
        DefaultAssay(obj) <- "SCT"
        gex_mat <- sum(exp(obj@assays$SCT@data[rownames(obj) == g, ])) /
          ncol(obj)
        DefaultAssay(obj) <- "GACT"
        gact_mat <- sum(exp(obj@assays$GACT@data[rownames(obj) == g, ])) /
          ncol(obj)
      } else {
        DefaultAssay(obj) <- "SCT"
        gex_to_append <- sum(exp(obj@assays$SCT@data[rownames(obj) == g, ])) /
          ncol(obj)
        gex_mat <- rbind(gex_mat, gex_to_append)
        DefaultAssay(obj) <- "GACT"
        gact_to_append <- sum(exp(obj@assays$GACT@data[rownames(obj) == g, ])) /
          ncol(obj)
        gact_mat <- rbind(gact_mat, gact_to_append)
      }
    }
    rownames(gex_mat) <- genelist
    rownames(gact_mat) <- genelist
    gex_df[, 3 * k - 3 + j] <- as.vector(gex_mat)
    gact_df[, 3 * k - 3 + j] <- as.vector(gact_mat)
  }
}


gex_clust <- hclust(dist(gex_df, method = "maximum"))
heatmap.2(x = as.matrix(gex_df),
          scale = "row", 
          trace = "none", 
          Rowv = gex_clust$order, 
          Colv = T, 
          dendrogram = "none",
          main = "gene expression")
heatmap.2(x = as.matrix(gact_df),
          scale = "row", 
          trace = "none", 
          Rowv = gex_clust$order, 
          Colv = T,
          #col = "bluered",
          dendrogram = "none",
          main = "gene activity")

gact_df_nonzero <- gact_df[rowSums(gact_df) != 0, ]
gex_df_nonzero <- gex_df[rowSums(gact_df) != 0, ]
gex_clust <- hclust(dist(gex_df_nonzero, method = "minkowski"))
heatmap.2(x = as.matrix(gex_df_nonzero),
          scale = "row", 
          trace = "none", 
          Rowv = gex_clust$order, 
          Colv = F, 
          dendrogram = "none",
          main = "gene expression")
heatmap.2(x = as.matrix(gact_df_nonzero),
          scale = "row", 
          trace = "none", 
          Rowv = gex_clust$order, 
          Colv = F, 
          dendrogram = "none",
          main = "gene activity")


# Check by cluster columns ####
obj_complete$clusterxtime.ident <- "others"
cls <- sort(unique(obj_complete$seurat_clusters))
for (k in 1:length(cls)) {
  cl <- cls[k]
  for (l in 1:length(times)) {
    t <- times[l]
    obj_complete$clusterxtime.ident[obj_complete$seurat_clusters == cl & 
                                      obj_complete$time.ident == t] <- paste(cl, t, sep = "_")
  }
}

cts <- sort(unique(obj_complete$clusterxtime.ident))

for (j in 1:length(cts)) {
  obj <- subset(obj_complete, clusterxtime.ident == cts[j])
  for (i in 1:length(genelist)) {
    g <- genelist[i]
    if (i == 1) {
      DefaultAssay(obj) <- "SCT"
      gex_col <- sum(exp(obj@assays$SCT@data[rownames(obj) == g, ])) /
        ncol(obj)
      DefaultAssay(obj) <- "GACT"
      gact_col <- sum(exp(obj@assays$GACT@data[rownames(obj) == g, ])) /
        ncol(obj)
    } else {
      DefaultAssay(obj) <- "SCT"
      gex_col_app <- sum(exp(obj@assays$SCT@data[rownames(obj) == g, ])) /
        ncol(obj)
      gex_col <- rbind(gex_col, gex_col_app)
      DefaultAssay(obj) <- "GACT"
      gact_col_app <- sum(exp(obj@assays$GACT@data[rownames(obj) == g, ])) /
        ncol(obj)
      gact_col <- rbind(gact_col, gact_col_app)
    }
  }
  rownames(gex_col) <- genelist
  rownames(gact_col) <- genelist
  if (j == 1) {
    print("entered if")
    gex_mat <- gex_col
    gact_mat <- gact_col
  } else {
    gex_mat <- cbind(gex_mat, gex_col)
    gact_mat <- cbind(gact_mat, gact_col)
  }
  print(j)
}
gact_df_nonzero <- gact_mat[rowSums(gact_mat) != 0, ]
gex_df_nonzero <- gex_mat[rowSums(gex_mat) != 0, ]
gex_clust <- hclust(dist(gex_df_nonzero, method = "euclidean"))
heatmap.2(x = as.matrix(gex_df_nonzero),
          scale = "row", 
          trace = "none", 
          #Rowv = T,
          Rowv = gex_clust$order, 
          Colv = T, 
          dendrogram = "column",
          col = "bluered",
          main = "gene expression")
heatmap.2(x = as.matrix(gact_df_nonzero),
          scale = "row", 
          trace = "none", 
          Rowv = gex_clust$order, 
          Colv = T, 
          dendrogram = "none",
          col = "bluered",
          main = "gene activity")

# change gene list
genelist <- union(DEG_paths_nmglut_up[[1]]$X1, DEG_paths_nmglut_up[[2]]$X1)
genelist <- union(genelist, DEG_paths_npglut_up[[1]]$X1)
genelist <- union(genelist,DEG_paths_npglut_up[[2]]$X1)
genelist <- union(genelist, DEG_paths_GABA_up[[1]]$X1)
genelist <- union(genelist, DEG_paths_GABA_up[[2]]$X1)
genelist <- genelist[sample.int(length(genelist), size = 2000, replace = F)]

# test within one time point
zero <- subset(obj_complete, time.ident == "0hr")
cs <- unique(zero$seurat_clusters)
for (j in 1:length(cs)) {
  obj <- subset(zero, seurat_clusters == cs[j])
  for (i in 1:length(genelist)) {
    g <- genelist[i]
    if (i == 1) {
      DefaultAssay(obj) <- "SCT"
      gex_col <- sum(exp(obj@assays$SCT@data[rownames(obj) == g, ])) /
        ncol(obj)
      DefaultAssay(obj) <- "GACT"
      gact_col <- sum(exp(obj@assays$GACT@data[rownames(obj) == g, ])) /
        ncol(obj)
    } else {
      DefaultAssay(obj) <- "SCT"
      gex_col_app <- sum(exp(obj@assays$SCT@data[rownames(obj) == g, ])) /
        ncol(obj)
      gex_col <- rbind(gex_col, gex_col_app)
      DefaultAssay(obj) <- "GACT"
      gact_col_app <- sum(exp(obj@assays$GACT@data[rownames(obj) == g, ])) /
        ncol(obj)
      gact_col <- rbind(gact_col, gact_col_app)
    }
  }
  rownames(gex_col) <- genelist
  rownames(gact_col) <- genelist
  if (j == 1) {
    print("entered if")
    gex_mat <- gex_col
    gact_mat <- gact_col
  } else {
    gex_mat <- cbind(gex_mat, gex_col)
    gact_mat <- cbind(gact_mat, gact_col)
  }
  print(j)
}
gact_df_nonzero <- gact_mat[rowSums(gact_mat) != 0, ]
gex_df_nonzero <- gex_mat[rownames(gex_mat) %in% rownames(gact_df_nonzero), ]
colnames(gex_df_nonzero) <- cs
colnames(gact_df_nonzero) <- cs
gex_clust_row <- hclust(dist(gex_df_nonzero, method = "euclidean"))
gex_clust_col <- hclust(dist(t(gex_df_nonzero), method = "euclidean"))
p<- heatmap.2(x = as.matrix(gex_df_nonzero),
          scale = "row", 
          trace = "none", 
          #Rowv = T,
          Rowv = as.dendrogram(gex_clust_row), 
          Colv = as.dendrogram(gex_clust_col), 
          dendrogram = "none",
          #col = "blue2yellow",
          col = c("blue", "black", "yellow"),
          breaks = c(-5, -0.5, 0.5, 5),
          main = "gene expression")
heatmap.2(x = as.matrix(gact_df_nonzero),
          scale = "row", 
          trace = "none", 
          Rowv = as.dendrogram(gex_clust_row), 
          Colv = as.dendrogram(gex_clust_col), 
          dendrogram = "none",
          col = c("blue", "black", "yellow"),
          breaks = c(-5, -0.5, 0.5, 5),
          main = "gene activity")

# use up and down regulated genes ####
DEG_paths <- list.files("../Analysis_RNAseq_pseudobulk/5line_DE_results/filtered_by_basemean/significant/", 
                        full.names = T, pattern = ".csv")
DEG_paths <- DEG_paths[str_detect(DEG_paths, "NPC", T)]
DEG_paths_up <- DEG_paths[str_detect(DEG_paths, "upregulated")]
DEG_paths_do <- DEG_paths[str_detect(DEG_paths, "downregulated")]
DEG_lists_up <- vector("list", length(DEG_paths_up))
DEG_lists_do <- vector("list", length(DEG_paths_do))

for (i in 1:length(DEG_lists_up)) {
  DEG_lists_up[[i]] <- read.csv(DEG_paths_up[i], row.names = 1)
  DEG_lists_do[[i]] <- read.csv(DEG_paths_do[i], row.names = 1)
  DEG_lists_up[[i]] <- arrange(DEG_lists_up[[i]], q_value)
  DEG_lists_do[[i]] <- arrange(DEG_lists_do[[i]], q_value)
}
for (i in 1:length(DEG_lists_up)) {
  if (i == 1) {
    genelist_up <- row.names(DEG_lists_up[[i]])
    genelist_do <- row.names(DEG_lists_do[[i]])
  } else {
    genelist_up <- c(genelist_up, row.names(DEG_lists_up[[i]]))
    genelist_do <- c(genelist_do, row.names(DEG_lists_do[[i]]))
  }
}
genelist_up <- unique(genelist_up)
genelist_do <- unique(genelist_do)
genelist_up <- genelist_up[1:1000]
genelist_do <- genelist_do[1:1000]
genelist <- c(genelist_up, genelist_do)

# Seurat heatmap ####
DoHeatmap(obj_complete, genelist, group.by = "clusterxtime.ident", slot = "data", assay = "SCT")
DoHeatmap(obj_complete, genelist, group.by = "clusterxtime.ident", slot = "data", assay = "GACT")
