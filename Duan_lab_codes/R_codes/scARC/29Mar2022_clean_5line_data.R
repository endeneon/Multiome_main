# Chuxuan Li 03/29/2022
# clean 5-line raw data to prepare for combination with 18-line data

# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)

library(stringr)
library(readr)
library(dplyr)
library(future)

library(ggplot2)
library(RColorBrewer)
library(viridis)
library(graphics)
library(patchwork)
library(cowplot)
library(gplots)

# set threads and parallelization
plan("multisession", workers = 1)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 107374182400)

# load 5 line data ####
# note this h5 file contains both atac-seq and gex information
mat_lst <- vector(mode = "list", length = 6L)
pathlist <- list.files(path = "/home/cli/NVME/scARC_Duan_018/hg38_Rnor6_mixed/",
                       pattern = "filtered_feature_bc_matrix.h5", 
                       full.names = T, recursive = T)
transformed_lst <- vector(mode = "list", length = 6L)

for (i in 1:length(pathlist)){
  h5file <- Read10X_h5(pathlist[i])
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(pathlist[i], "[2|8]_[0|1|6]"))
  print(paste0(i, " number of genes: ", nrow(obj), 
               ", number of cells: ", ncol(obj)))
  obj <- 
    PercentageFeatureSet(obj,
                         pattern = c("^MT-"),
                         col.name = "percent.mt")
  transformed_lst[[i]] <- SCTransform(obj, 
                                      vars.to.regress = "percent.mt",
                                      method = "glmGamPoi",
                                      variable.features.n = 8000,
                                      seed.use = 115,
                                      verbose = T)
}

# dimensional reduction and clustering ####
dimreduclust <- function(obj){
  obj <- ScaleData(obj)
  obj <- RunPCA(obj,
                verbose = T,
                seed.use = 11)
  obj <- RunUMAP(obj,
                 reduction = "pca",
                 dims = 1:30,
                 seed.use = 11)
  obj <- FindNeighbors(obj,
                       reduction = "pca",
                       dims = 1:30)
  obj <- FindClusters(obj,
                      resolution = 0.5,
                      random.seed = 11)
  return(obj)
  
}
for (i in 1:length(transformed_lst)){
  transformed_lst[[i]] <- dimreduclust(transformed_lst[[i]])
}

# check and remove rat cells ####
cleanobj_lst <- vector(mode = "list", length = length(transformed_lst))

setwd("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v2_combine_5line_18line/5line_QC_plots")
for(i in 1:length(pathlist)){
  jpeg(paste0(str_extract(pathlist[i], "[2|8]_[0|1|6]"), "_dimplot_by_cluster.jpeg"))
  p <- DimPlot(transformed_lst[[i]],
               label = T,
               group.by = "seurat_clusters") +
    NoLegend() +
    ggtitle("2_0") +
    theme(text = element_text(size = 12))
  print(p)
  dev.off()
  
  jpeg(paste0(str_extract(pathlist[i], "[2|8]_[0|1|6]"), "_Gfap.jpeg"))
  p <- FeaturePlot(transformed_lst[[i]], features = "Gfap") +
    ggtitle(str_extract(pathlist[i], "[2|8]_[0|1|6]"))
  print(p)
  dev.off()
  
  jpeg(paste0(str_extract(pathlist[i], "[2|8]_[0|1|6]"), "_GAD1.jpeg"))
  p <- FeaturePlot(transformed_lst[[i]], features = "GAD1") +
    ggtitle(str_extract(pathlist[i], "[2|8]_[0|1|6]"))
  print(p)
  dev.off()
  
  jpeg(paste0(str_extract(pathlist[i], "[2|8]_[0|1|6]"), "_SLC17A6.jpeg"))
  p <- FeaturePlot(transformed_lst[[i]], features = "SLC17A6") +
    ggtitle(str_extract(pathlist[i], "[2|8]_[0|1|6]"))
  print(p)
  dev.off()
  
  jpeg(paste0(str_extract(pathlist[i], "[2|8]_[0|1|6]"), "_Gfap_vln.jpeg"))
  p <- VlnPlot(transformed_lst[[i]], "Gfap") +
    ggtitle(paste0(str_extract(pathlist[i], "[2|8]_[0|1|6]"), " Gfap expression by cluster"))
  print(p)
  dev.off()
}

# remove clusters that are rat cells, count
transformed_lst[[6]]$rat.ident <- "human"
transformed_lst[[1]]$rat.ident[transformed_lst[[1]]$seurat_clusters %in% c(3, 5, 9, 17)] <- "rat"
transformed_lst[[2]]$rat.ident[transformed_lst[[2]]$seurat_clusters %in% c(5, 6, 12)] <- "rat"
transformed_lst[[3]]$rat.ident[transformed_lst[[3]]$seurat_clusters %in% c(4, 5, 10)] <- "rat"
transformed_lst[[4]]$rat.ident[transformed_lst[[4]]$seurat_clusters %in% c(3)] <- "rat"
transformed_lst[[5]]$rat.ident[transformed_lst[[5]]$seurat_clusters %in% c(5, 13, 15, 16)] <- "rat"
transformed_lst[[6]]$rat.ident[transformed_lst[[6]]$seurat_clusters %in% c(3)] <- "rat"
sum(transformed_lst[[6]]$rat.ident ==  "rat")
cleanobj_lst[[6]] <- subset(transformed_lst[[6]], subset = rat.ident != "rat")
length(unique(cleanobj_lst[[6]]$seurat_clusters))

setwd("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v2_combine_5line_18line")
save(cleanobj_lst, file = "5line_removed_rat_list.RData")

# remove rat genes ####
# load 18-line objects to match the gene list
load("../Analysis_v1_normalize_by_lib_then_integrate/demuxed_obj_removed_rat_gene_after_integration_nfeat_3000_ftoi_features_updated.RData")
use.genes <- integrated@assays$RNA@counts@Dimnames[[1]]
finalobj_lst <- vector(mode = "list", length = length(cleanobj_lst))
# import the 36601 genes from previous analysis
for(i in 1:length(finalobj_lst)){
  prev.genes <- finalobj_lst[[i]]@assays$RNA@counts@Dimnames[[1]]
  use.genes <- intersect(human.genes, prev.genes)
  counts <- GetAssayData(finalobj_lst[[i]], assay = "RNA")
  counts <- counts[(which(rownames(counts) %in% use.genes)), ]
  finalobj_lst[[i]] <- subset(finalobj_lst[[i]], features = rownames(counts))
}

save(finalobj_lst, file = "5line_removed_rat_genes_list.RData")

# map to demuxed barcodes
# assign cell line identities
times <- c("0", "1", "6")
lines_g2 <- c("CD_27", "CD_54")
lines_g8 <- c("CD_08", "CD_25", "CD_26")
g_2_0_bc_lst <- vector(mode = "list", length = 2)
g_2_1_bc_lst <- vector(mode = "list", length = 2)
g_2_6_bc_lst <- vector(mode = "list", length = 2)
g_8_0_bc_lst <- vector(mode = "list", length = 3)
g_8_1_bc_lst <- vector(mode = "list", length = 3)
g_8_6_bc_lst <- vector(mode = "list", length = 3)

# import cell ident barcodes
for (t in times){
  for (i in 1:length(lines_g2)){
    l <- lines_g2[i]
    file_name <- paste0("/nvmefs/scARC_Duan_018/GRCh38_mapped_only/common_barcodes/g_2_",
                        t,
                        "_common_CD27_CD54.best.",
                        l)
    print(file_name)
    bc <- read_csv(file_name, col_names = FALSE)
    bc <- unlist(bc)
    bc <- str_sub(bc, end = -3L)
    if (t %in% "0"){
      g_2_0_bc_lst[[i]] <- bc
    } else if (t %in% "1"){
      g_2_1_bc_lst[[i]] <- bc
    } else {
      g_2_6_bc_lst[[i]] <- bc
    }
  }
}

for (t in times){
  for (i in 1:length(lines_g8)){
    l <- lines_g8[i]
    file_name <- paste0("/nvmefs/scARC_Duan_018/GRCh38_mapped_only/common_barcodes/g_8_",
                        t,
                        "_common_CD08_CD25_CD26.best.",
                        l)
    print(file_name)
    bc <- read_csv(file_name, col_names = FALSE)
    bc <- unlist(bc)
    bc <- str_sub(bc, end = -3L)
    if (t %in% "0"){
      g_8_0_bc_lst[[i]] <- bc
    } else if (t %in% "1"){
      g_8_1_bc_lst[[i]] <- bc
    } else {
      g_8_6_bc_lst[[i]] <- bc
    }
  }
}

# assign cell line identities to the objects
for (i in 1:length(finalobj_lst)){
  finalobj_lst[[i]]$cell.line.ident <- "unmatched"
  finalobj_lst[[i]]$cell.line.ident[str_sub(finalobj_lst[[i]]@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                                              c(g_2_0_bc_lst[[1]], g_2_1_bc_lst[[1]], g_2_6_bc_lst[[1]])] <- "CD_27"
  finalobj_lst[[i]]$cell.line.ident[str_sub(finalobj_lst[[i]]@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                                              c(g_2_0_bc_lst[[2]], g_2_1_bc_lst[[2]], g_2_6_bc_lst[[2]])] <- "CD_54"
  
  finalobj_lst[[i]]$cell.line.ident[str_sub(finalobj_lst[[i]]@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                                              c(g_8_0_bc_lst[[1]], g_8_1_bc_lst[[1]], g_8_6_bc_lst[[1]])] <- "CD_08"
  finalobj_lst[[i]]$cell.line.ident[str_sub(finalobj_lst[[i]]@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                                              c(g_8_0_bc_lst[[2]], g_8_1_bc_lst[[2]], g_8_6_bc_lst[[2]])] <- "CD_25"
  finalobj_lst[[i]]$cell.line.ident[str_sub(finalobj_lst[[i]]@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                                              c(g_8_0_bc_lst[[3]], g_8_1_bc_lst[[3]], g_8_6_bc_lst[[3]])] <- "CD_26"
  print(sum(finalobj_lst[[i]]$cell.line.ident == "unmatched"))
  finalobj_lst[[i]] <- subset(finalobj_lst[[i]], cell.line.ident != "unmatched")
}

save(finalobj_lst, file = "5line_mapped_to_demuxed_barcodes_human_mapped_only_list.RData")
