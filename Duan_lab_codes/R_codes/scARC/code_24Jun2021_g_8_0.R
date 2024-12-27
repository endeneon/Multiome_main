# Siwei 24 Jun 2021
# import GEX data of g_8_0 (3 samples) for t-sne plotting

# init
library(Seurat)
library(sctransform)

library(data.table)
library(ggplot2)

# library(Matrix)
library(Rfast)
library(plyr)
library(dplyr)
library(stringr)
library(future)

library(readr)

# set threads and parallelization
# this one will cause Cstack error if used with SCTransform() or NormalizeData()
plan("multisession", workers = 16)
plan()
options(expressions = 20000)

# load data
# load data use read10x_h5
# note this h5 file contains both atac-seq and gex information
g_8_0_raw <- 
  Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_8_0/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_8_0_raw_gex <- 
  CreateSeuratObject(counts = g_8_0_raw$`Gene Expression`,
                     project = "g_8_0_raw_gex")

nrow(g_8_0_raw_gex) # 91461 genes (hg38 + Rnor6)
ncol(g_8_0_raw_gex) # 15026 cells

# make a seurat object for seubsequent operation
# Note: all genes with total counts = 0 have been pre-removed
g_8_0_gex <- g_8_0_raw_gex

# store mitochondrial percentage in object meta data
# human starts with "MT-", rat starts with "Mt-"
g_8_0_gex <- 
  PercentageFeatureSet(g_8_0_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")

# pre-treat data 
# use SCTransform()
# Note: SCTransform has to be downloaded and installed directly from CRAN
# since it has to be compiled locally, do not use the anaconda distribution
# including NormalizeData(), ScaleData(), and FindVariableFeatures
g_8_0_gex <- 
  SCTransform(g_8_0_gex, 
              vars.to.regress = "percent.mt",
              method = "glmGamPoi",
              verbose = T)

g_8_0_gex <- RunPCA(g_8_0_gex, 
                    verbose = T)
g_8_0_gex <- RunUMAP(g_8_0_gex,
                     dims = 1:30,
                     verbose = T)
g_8_0_gex <- FindNeighbors(g_8_0_gex,
                           dims = 1:30,
                           verbose = T)
g_8_0_gex <- FindClusters(g_8_0_gex,
                          verbose = T)

# make dimension plot
DimPlot(g_8_0_gex, 
        label = T) +
  NoLegend()

# check specific genes
FeaturePlot(g_8_0_gex, 
            features = c("GAD1", "SLC17A6", 
                         "Gfap", "SOX2"), 
            pt.size = 0.2,
            ncol = 2)

# assign cell line identities
## add one column of cell line identities to meta.data
g_8_0_gex@meta.data$cell.line.ident <- "unknown"

## import cell ident barcodes
CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/backup/g_8_0_CD08_barcodes.txt", 
           col_names = FALSE)
CD_08_barcodes <- unlist(CD_08_barcodes)

CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/backup/g_8_0_CD25_barcodes.txt", 
           col_names = FALSE)
CD_25_barcodes <- unlist(CD_25_barcodes)

CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/backup/g_8_0_CD26_barcodes.txt", 
           col_names = FALSE)
CD_26_barcodes <- unlist(CD_26_barcodes)

## assign cell line identity
g_8_0_gex@meta.data$cell.line.ident[g_8_0_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      CD_08_barcodes] <- "CD_08"
g_8_0_gex@meta.data$cell.line.ident[g_8_0_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      CD_25_barcodes] <- "CD_25"
g_8_0_gex@meta.data$cell.line.ident[g_8_0_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      CD_26_barcodes] <- "CD_26"
g_8_0_gex@meta.data$cell.line.ident <-
  factor(g_8_0_gex@meta.data$cell.line.ident)

# check cell distribution
# check specific genes
DimPlot(g_8_0_gex, 
        label = T) +
  NoLegend()

FeaturePlot(g_8_0_gex, 
            features = c("GAD1", "SLC17A6", 
                         "Gfap", "SOX2"), 
            pt.size = 0.2,
            ncol = 2)

FeaturePlot(g_8_0_gex, 
            features = c("SLC17A6", "SLC17A7"), 
            # split.by = "cell.line.ident",
            pt.size = 0.2,
            ncol = 2)


DimPlot(g_8_0_gex,
        cols = c("darkred", "green", "blue", "black"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_8_0_gex@meta.data$cell.line.ident,
        group.by = "cell.line.ident",
        label = F)# +
  NoLegend()

  
sum(g_8_0_gex@meta.data$seurat_clusters %in% 
      c(21, 9, 0, 10, 13, 14, 8, 4, 7))
sum(g_8_0_gex@meta.data$seurat_clusters %in% 
      c(2, 6, 5, 20))

sum((g_8_0_gex@meta.data$seurat_clusters %in% 
       c(21, 9, 0, 10, 13, 14,
         18, 4, 7, 
         6, 5, 20, 2)) & 
    (g_8_0_gex@meta.data$cell.line.ident %in% c("unknown")))

sum((g_8_0_gex@meta.data$seurat_clusters %in% 
       c(21, 9, 0, 10, 13, 14,
         18, 4, 7, 
         6, 5, 20, 2)))
    
sum((g_8_0_gex@meta.data$seurat_clusters %in% 
       c(1)))
