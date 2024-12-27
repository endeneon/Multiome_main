# Chuxuan 9/14/2021
# Using GEX data from group 2 and group 8 to generate barcodes by line, then extract the barcodes of bad cells

rm(list=ls())
# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
### There might be conflicts between packages,
### if running SCTransform, only load the three libraries above
##############################################################



library(readr)

# library(data.table)
library(ggplot2)

# library(Matrix)
# library(Rfast) # This one is incompatible with Seurat SCTransform(), load separately!
# library(plyr)
# library(dplyr)
# library(stringr)


library(future)
# set threads and parallelization

plan("multisession", workers = 8)
# plan("sequential")
plan()
options(expressions = 20000)

# load data
# load data use read10x_h5
# note this h5 file contains both atac-seq and gex information
g_2_0_raw <- 
  Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_2_0/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_2_0_raw_gex <- 
  CreateSeuratObject(counts = g_2_0_raw$`Gene Expression`,
                     project = "g_2_0_raw_gex")

g_2_1_raw <- 
  Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_2_1/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_2_1_raw_gex <- 
  CreateSeuratObject(counts = g_2_1_raw$`Gene Expression`,
                     project = "g_2_1_raw_gex")

g_2_6_raw <- 
  Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_2_6/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_2_6_raw_gex <- 
  CreateSeuratObject(counts = g_2_6_raw$`Gene Expression`,
                     project = "g_2_6_raw_gex")


g_8_0_raw <- 
  Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_8_0/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_8_0_raw_gex <- 
  CreateSeuratObject(counts = g_8_0_raw$`Gene Expression`,
                     project = "g_8_0_raw_gex")

g_8_1_raw <- 
  Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_8_1/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_8_1_raw_gex <- 
  CreateSeuratObject(counts = g_8_1_raw$`Gene Expression`,
                     project = "g_8_1_raw_gex")

g_8_6_raw <- 
  Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_8_6/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_8_6_raw_gex <- 
  CreateSeuratObject(counts = g_8_6_raw$`Gene Expression`,
                     project = "g_8_6_raw_gex")


# make a seurat object for subsequent operation
# Note: all genes with total counts = 0 have been pre-removed
g_2_0_gex <- g_2_0_raw_gex
g_2_1_gex <- g_2_1_raw_gex
g_2_6_gex <- g_2_6_raw_gex
g_8_0_gex <- g_8_0_raw_gex
g_8_1_gex <- g_8_1_raw_gex
g_8_6_gex <- g_8_6_raw_gex


# store mitochondrial percentage in object meta data
# human starts with "MT-", rat starts with "Mt-"
g_2_0_gex <- 
  PercentageFeatureSet(g_2_0_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")
g_2_1_gex <- 
  PercentageFeatureSet(g_2_1_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")
g_2_6_gex <- 
  PercentageFeatureSet(g_2_6_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")

g_8_0_gex <- 
  PercentageFeatureSet(g_8_0_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")
g_8_1_gex <- 
  PercentageFeatureSet(g_8_1_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")
g_8_6_gex <- 
  PercentageFeatureSet(g_8_6_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")

VlnPlot(g_2_0_gex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(g_2_1_gex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(g_2_6_gex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(g_8_0_gex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(g_8_1_gex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(g_8_6_gex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

sum(g_2_0_gex@meta.data$percent.mt > 10) # 476 cells
sum(g_2_1_gex@meta.data$percent.mt > 10) # 347 cells
sum(g_2_6_gex@meta.data$percent.mt > 10) # 353 cells
sum(g_8_0_gex@meta.data$percent.mt > 10) # 576 cells
sum(g_8_1_gex@meta.data$percent.mt > 10) # 548 cells
sum(g_8_6_gex@meta.data$percent.mt > 10) # 290 cells

# exclude cells with percent.mt greater than 10
g_2_0_gex <- subset(g_2_0_gex, subset = percent.mt < 10)
g_2_1_gex <- subset(g_2_1_gex, subset = percent.mt < 10)
g_2_6_gex <- subset(g_2_6_gex, subset = percent.mt < 10)
g_8_0_gex <- subset(g_8_0_gex, subset = percent.mt < 10)
g_8_1_gex <- subset(g_8_1_gex, subset = percent.mt < 10)
g_8_6_gex <- subset(g_8_6_gex, subset = percent.mt < 10)



# pre-treat data 
# use SCTransform()
# Note: SCTransform has to be downloaded and installed directly from CRAN
# since it has to be compiled locally, do not use the anaconda distribution
# including NormalizeData(), ScaleData(), and FindVariableFeatures
g_2_0_gex <- 
  SCTransform(g_2_0_gex, 
              vars.to.regress = "percent.mt", 
              variable.features.n = 8000,
              method = "glmGamPoi",
              verbose = T)

g_2_0_gex <- RunPCA(g_2_0_gex, 
                    verbose = T)
g_2_0_gex <- RunUMAP(g_2_0_gex,
                     dims = 1:30,
                     verbose = T)
g_2_0_gex <- FindNeighbors(g_2_0_gex,
                           dims = 1:30,
                           verbose = T)
g_2_0_gex <- FindClusters(g_2_0_gex,
                          verbose = T)


g_2_1_gex <- 
  SCTransform(g_2_1_gex, 
              vars.to.regress = "percent.mt", 
              variable.features.n = 8000,
              method = "glmGamPoi",
              verbose = T)

g_2_1_gex <- RunPCA(g_2_1_gex, 
                    verbose = T)
g_2_1_gex <- RunUMAP(g_2_1_gex,
                     dims = 1:30,
                     verbose = T)
g_2_1_gex <- FindNeighbors(g_2_1_gex,
                           dims = 1:30,
                           verbose = T)
g_2_1_gex <- FindClusters(g_2_1_gex,
                          verbose = T)


g_2_6_gex <- 
  SCTransform(g_2_6_gex, 
              vars.to.regress = "percent.mt", 
              variable.features.n = 8000,
              method = "glmGamPoi",
              verbose = T)

g_2_6_gex <- RunPCA(g_2_6_gex, 
                    verbose = T)
g_2_6_gex <- RunUMAP(g_2_6_gex,
                     dims = 1:30,
                     verbose = T)
g_2_6_gex <- FindNeighbors(g_2_6_gex,
                           dims = 1:30,
                           verbose = T)
g_2_6_gex <- FindClusters(g_2_6_gex,
                          verbose = T)


g_8_0_gex <- 
  SCTransform(g_8_0_gex, 
              vars.to.regress = "percent.mt", 
              variable.features.n = 8000,
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


g_8_1_gex <- 
  SCTransform(g_8_1_gex, 
              vars.to.regress = "percent.mt", 
              variable.features.n = 8000,
              method = "glmGamPoi",
              verbose = T)

g_8_1_gex <- RunPCA(g_8_1_gex, 
                    verbose = T)
g_8_1_gex <- RunUMAP(g_8_1_gex,
                     dims = 1:30,
                     verbose = T)
g_8_1_gex <- FindNeighbors(g_8_1_gex,
                           dims = 1:30,
                           verbose = T)
g_8_1_gex <- FindClusters(g_8_1_gex,
                          verbose = T)


g_8_6_gex <- 
  SCTransform(g_8_6_gex, 
              vars.to.regress = "percent.mt", 
              variable.features.n = 8000,
              method = "glmGamPoi",
              verbose = T)

g_8_6_gex <- RunPCA(g_8_6_gex, 
                    verbose = T)
g_8_6_gex <- RunUMAP(g_8_6_gex,
                     dims = 1:30,
                     verbose = T)
g_8_6_gex <- FindNeighbors(g_8_6_gex,
                           dims = 1:30,
                           verbose = T)
g_8_6_gex <- FindClusters(g_8_6_gex,
                          verbose = T)



# assign cell line identities
## add one column of cell line identities to meta.data
g_2_0_gex@meta.data$cell.line.ident <- "unknown"

# assign cell line identities
## add one column of cell line identities to meta.data
g_2_1_gex@meta.data$cell.line.ident <- "unknown"

# assign cell line identities
## add one column of cell line identities to meta.data
g_2_6_gex@meta.data$cell.line.ident <- "unknown"

# assign cell line identities
## add one column of cell line identities to meta.data
g_8_0_gex@meta.data$cell.line.ident <- "unknown"

# assign cell line identities
## add one column of cell line identities to meta.data
g_8_1_gex@meta.data$cell.line.ident <- "unknown"

# assign cell line identities
## add one column of cell line identities to meta.data
g_8_6_gex@meta.data$cell.line.ident <- "unknown"


## import cell ident barcodes
g_2_0_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_0_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
g_2_0_CD_27_barcodes <- unlist(g_2_0_CD_27_barcodes)
g_2_0_CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_0_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
g_2_0_CD_54_barcodes <- unlist(g_2_0_CD_54_barcodes)


g_2_1_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_1_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
g_2_1_CD_27_barcodes <- unlist(g_2_1_CD_27_barcodes)
g_2_1_CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_1_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
g_2_1_CD_54_barcodes <- unlist(g_2_1_CD_54_barcodes)

g_2_6_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_6_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
g_2_6_CD_27_barcodes <- unlist(g_2_6_CD_27_barcodes)
g_2_6_CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_6_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
g_2_6_CD_54_barcodes <- unlist(g_2_6_CD_54_barcodes)


# import group 8 cell line barcodes
g_8_0_CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_0_common_CD08_CD25_CD26.best.CD_08", 
           col_names = FALSE)
g_8_0_CD_08_barcodes <- unlist(g_8_0_CD_08_barcodes)

g_8_0_CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_0_common_CD08_CD25_CD26.best.CD_25", 
           col_names = FALSE)
g_8_0_CD_25_barcodes <- unlist(g_8_0_CD_25_barcodes)

g_8_0_CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_0_common_CD08_CD25_CD26.best.CD_26", 
           col_names = FALSE)
g_8_0_CD_26_barcodes <- unlist(g_8_0_CD_26_barcodes)


g_8_1_CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_1_common_CD08_CD25_CD26.best.CD_08", 
           col_names = FALSE)
g_8_1_CD_08_barcodes <- unlist(g_8_1_CD_08_barcodes)

g_8_1_CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_1_common_CD08_CD25_CD26.best.CD_25", 
           col_names = FALSE)
g_8_1_CD_25_barcodes <- unlist(g_8_1_CD_25_barcodes)

g_8_1_CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_1_common_CD08_CD25_CD26.best.CD_26", 
           col_names = FALSE)
g_8_1_CD_26_barcodes <- unlist(g_8_1_CD_26_barcodes)


g_8_6_CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_6_common_CD08_CD25_CD26.best.CD_08", 
           col_names = FALSE)
g_8_6_CD_08_barcodes <- unlist(g_8_6_CD_08_barcodes)

g_8_6_CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_6_common_CD08_CD25_CD26.best.CD_25", 
           col_names = FALSE)
g_8_6_CD_25_barcodes <- unlist(g_8_6_CD_25_barcodes)

g_8_6_CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_6_common_CD08_CD25_CD26.best.CD_26", 
           col_names = FALSE)
g_8_6_CD_26_barcodes <- unlist(g_8_6_CD_26_barcodes)


# 
# sum(g_2_6_gex@assays$SCT@counts@Dimnames[[2]] %in%
#       CD_27_barcodes)
#sum(g_2_6_gex@assays$SCT@counts@Dimnames[[2]] %in%
#      CD_54_barcodes)

## assign cell line identity
g_2_0_gex@meta.data$cell.line.ident[g_2_0_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_2_0_CD_27_barcodes] <- "CD_27"
g_2_0_gex@meta.data$cell.line.ident[g_2_0_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_2_0_CD_54_barcodes] <- "CD_54"

g_2_0_gex@meta.data$cell.line.ident <-
  factor(g_2_0_gex@meta.data$cell.line.ident)

summary(g_2_0_gex@meta.data$cell.line.ident)


g_2_1_gex@meta.data$cell.line.ident[g_2_1_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_2_1_CD_27_barcodes] <- "CD_27"
g_2_1_gex@meta.data$cell.line.ident[g_2_1_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_2_1_CD_54_barcodes] <- "CD_54"

g_2_1_gex@meta.data$cell.line.ident <-
  factor(g_2_1_gex@meta.data$cell.line.ident)

summary(g_2_1_gex@meta.data$cell.line.ident)


g_2_6_gex@meta.data$cell.line.ident[g_2_6_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_2_6_CD_27_barcodes] <- "CD_27"
g_2_6_gex@meta.data$cell.line.ident[g_2_6_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_2_6_CD_54_barcodes] <- "CD_54"

g_2_6_gex@meta.data$cell.line.ident <-
  factor(g_2_6_gex@meta.data$cell.line.ident)

summary(g_2_6_gex@meta.data$cell.line.ident)



g_8_0_gex@meta.data$cell.line.ident[g_8_0_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_8_0_CD_08_barcodes] <- "CD_08"
g_8_0_gex@meta.data$cell.line.ident[g_8_0_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_8_0_CD_25_barcodes] <- "CD_25"
g_8_0_gex@meta.data$cell.line.ident[g_8_0_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_8_0_CD_26_barcodes] <- "CD_26"
g_8_0_gex@meta.data$cell.line.ident <-
  factor(g_8_0_gex@meta.data$cell.line.ident)

summary(g_8_0_gex@meta.data$cell.line.ident)

g_8_1_gex@meta.data$cell.line.ident[g_8_1_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_8_1_CD_08_barcodes] <- "CD_08"
g_8_1_gex@meta.data$cell.line.ident[g_8_1_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_8_1_CD_25_barcodes] <- "CD_25"
g_8_1_gex@meta.data$cell.line.ident[g_8_1_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_8_1_CD_26_barcodes] <- "CD_26"
g_8_1_gex@meta.data$cell.line.ident <-
  factor(g_8_1_gex@meta.data$cell.line.ident)

summary(g_8_1_gex@meta.data$cell.line.ident)

g_8_6_gex@meta.data$cell.line.ident[g_8_6_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_8_6_CD_08_barcodes] <- "CD_08"
g_8_6_gex@meta.data$cell.line.ident[g_8_6_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_8_6_CD_25_barcodes] <- "CD_25"
g_8_6_gex@meta.data$cell.line.ident[g_8_6_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      g_8_6_CD_26_barcodes] <- "CD_26"
g_8_6_gex@meta.data$cell.line.ident <-
  factor(g_8_6_gex@meta.data$cell.line.ident)

summary(g_8_6_gex@meta.data$cell.line.ident)


# check cell distribution
# check specific genes
DimPlot(g_2_0_gex, 
        label = T) +
  NoLegend()+
  ggtitle("group_2_0hr")

FeaturePlot(g_2_0_gex,
            features = c("GAD1", "SLC17A6",
                         "Gfap", "SOX2"),
            pt.size = 0.2,
            ncol = 2)

DimPlot(g_2_0_gex,
        cols = c("darkred", "green", "black"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "cell.line.ident",
        label = F) +
  ggtitle("group_2_0hr")

# confirm some clusters have low counts and low features 
g_2_0_gex@meta.data$log_nCount_RNA <- log(g_2_0_gex@meta.data$nCount_RNA)
g_2_0_gex@meta.data$log_nFeature_RNA <- log(g_2_0_gex@meta.data$nFeature_RNA)

FeaturePlot(g_2_0_gex, 
            features = "log_nCount_RNA", 
            pt.size = 0.2)

FeaturePlot(g_2_0_gex, 
            features = "log_nFeature_RNA", 
            pt.size = 0.2)

FeaturePlot(g_2_0_gex, 
            features = "percent.mt", 
            pt.size = 0.2)

# write out barcodes for rat astrocytes and bad quality cells
g_2_0_bad_barcodes <- g_2_0_gex@assays$RNA@counts@Dimnames[[2]][g_2_0_gex@meta.data$seurat_clusters %in% c("10", "20", "16", "18", 
                                                                                                           "21", "5", "13", "3")]
write.table(g_2_0_bad_barcodes,
            file = "QC_barcodes/g_2_0_exclude_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

# repeat for g21, g26, group 8
DimPlot(g_2_1_gex, 
        label = T) +
  NoLegend() +
  ggtitle("group_2_1hr")

FeaturePlot(g_2_1_gex,
            features = c("GAD1", "SLC17A6",
                         "Gfap", "SOX2"),
            pt.size = 0.2,
            ncol = 2)

DimPlot(g_2_1_gex,
        cols = c("darkred", "green", "black"),
        group.by = "cell.line.ident",
        label = F) +
  ggtitle("group_2_1hr")

# confirm some clusters have low counts and low features 
g_2_1_gex@meta.data$log_nCount_RNA <- log(g_2_1_gex@meta.data$nCount_RNA)
g_2_1_gex@meta.data$log_nFeature_RNA <- log(g_2_1_gex@meta.data$nFeature_RNA)

FeaturePlot(g_2_1_gex, 
            features = "log_nCount_RNA", 
            pt.size = 0.2)

FeaturePlot(g_2_1_gex, 
            features = "log_nFeature_RNA", 
            pt.size = 0.2)

FeaturePlot(g_2_1_gex, 
            features = "percent.mt", 
            pt.size = 0.2)

# write out barcodes for rat astrocytes and bad quality cells
g_2_1_bad_barcodes <- g_2_1_gex@assays$SCT@counts@Dimnames[[2]][g_2_1_gex@meta.data$seurat_clusters %in% c("19", "6", "13", "3", 
                                                                                                           "14", "12")]
write.table(g_2_1_bad_barcodes,
            file = "QC_barcodes/g_2_1_exclude_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)


DimPlot(g_2_6_gex, 
        label = T) +
  NoLegend() +
  ggtitle("group_2_6hr")

FeaturePlot(g_2_6_gex,
            features = c("GAD1", "SLC17A6",
                         "Gfap", "SOX2"),
            pt.size = 0.2,
            ncol = 2)

DimPlot(g_2_6_gex,
        cols = c("darkred", "green", "black"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "cell.line.ident",
        label = F) +
  ggtitle("group_2_6hr")

# confirm some clusters have low counts and low features 
g_2_6_gex@meta.data$log_nCount_RNA <- log(g_2_6_gex@meta.data$nCount_RNA)
g_2_6_gex@meta.data$log_nFeature_RNA <- log(g_2_6_gex@meta.data$nFeature_RNA)

FeaturePlot(g_2_6_gex, 
            features = "log_nCount_RNA", 
            pt.size = 0.2)

FeaturePlot(g_2_6_gex, 
            features = "log_nFeature_RNA", 
            pt.size = 0.2)

FeaturePlot(g_2_6_gex, 
            features = "percent.mt", 
            pt.size = 0.2)

# write out barcodes for rat astrocytes and bad quality cells
g_2_6_bad_barcodes <- g_2_6_gex@assays$SCT@counts@Dimnames[[2]][g_2_6_gex@meta.data$seurat_clusters %in% c("20", "19", "16", "21",
                                                                                                           "2", "11", "5", "13")]
write.table(g_2_6_bad_barcodes,
            file = "QC_barcodes/g_2_6_exclude_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)



# check cell distribution
# check specific genes
DimPlot(g_8_0_gex, 
        label = T) +
  NoLegend()+
  ggtitle("group_8_0hr")

FeaturePlot(g_8_0_gex,
            features = c("GAD1", "SLC17A6",
                         "Gfap", "SOX2"),
            pt.size = 0.2,
            ncol = 2)

DimPlot(g_8_0_gex,
        cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "cell.line.ident",
        label = F) +
  ggtitle("group_8_0hr")

# confirm some clusters have low counts and low features 
g_8_0_gex@meta.data$log_nCount_RNA <- log(g_8_0_gex@meta.data$nCount_RNA)
g_8_0_gex@meta.data$log_nFeature_RNA <- log(g_8_0_gex@meta.data$nFeature_RNA)

FeaturePlot(g_8_0_gex, 
            features = "log_nCount_RNA", 
            pt.size = 0.2)

FeaturePlot(g_8_0_gex, 
            features = "log_nFeature_RNA", 
            pt.size = 0.2)

FeaturePlot(g_8_0_gex, 
            features = "percent.mt", 
            pt.size = 0.2)

# write out barcodes for rat astrocytes and bad quality cells
g_8_0_bad_barcodes <- g_8_0_gex@assays$SCT@counts@Dimnames[[2]][g_8_0_gex@meta.data$seurat_clusters %in% c("0", "6", "11", "18", "19")]
write.table(g_8_0_bad_barcodes,
            file = "QC_barcodes/g_8_0_exclude_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

# repeat for g21, g26, group 8
DimPlot(g_8_1_gex, 
        label = T) +
  NoLegend()+
  ggtitle("group_8_1hr")

FeaturePlot(g_8_1_gex,
            features = c("GAD1", "SLC17A6",
                         "Gfap", "SOX2"),
            pt.size = 0.2,
            ncol = 2)

# FeaturePlot(g_8_1_gex,
#             features = c("SLC17A7"),
#             pt.size = 0.2)

DimPlot(g_8_1_gex,
        cols = c("cyan", "magenta", "yellow", "black"),
        group.by = "cell.line.ident",
        label = F) +
  ggtitle("group_8_1hr")

# confirm some clusters have low counts and low features 
g_8_1_gex@meta.data$log_nCount_RNA <- log(g_8_1_gex@meta.data$nCount_RNA)
g_8_1_gex@meta.data$log_nFeature_RNA <- log(g_8_1_gex@meta.data$nFeature_RNA)

FeaturePlot(g_8_1_gex, 
            features = "log_nCount_RNA", 
            pt.size = 0.2)

FeaturePlot(g_8_1_gex, 
            features = "log_nFeature_RNA", 
            pt.size = 0.2)

FeaturePlot(g_8_1_gex, 
            features = "percent.mt", 
            pt.size = 0.2)

# write out barcodes for rat astrocytes and bad quality cells
g_8_1_bad_barcodes <- g_8_1_gex@assays$SCT@counts@Dimnames[[2]][g_8_1_gex@meta.data$seurat_clusters %in% c("5", "17", "3", "19",
                                                                                                           "10", "20")]
write.table(g_8_1_bad_barcodes,
            file = "QC_barcodes/g_8_1_exclude_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)


DimPlot(g_8_6_gex, 
        label = T) +
  NoLegend()+
  ggtitle("group_8_6hr")

FeaturePlot(g_8_6_gex,
            features = c("GAD1", "SLC17A6",
                         "Gfap", "SOX2"),
            pt.size = 0.2,
            ncol = 2)

DimPlot(g_8_6_gex,
        cols = c("cyan", "magenta", "yellow", "black"),
        group.by = "cell.line.ident",
        label = F) +
  ggtitle("group_8_6hr")

# confirm some clusters have low counts and low features 
g_8_6_gex@meta.data$log_nCount_RNA <- log(g_8_6_gex@meta.data$nCount_RNA)
g_8_6_gex@meta.data$log_nFeature_RNA <- log(g_8_6_gex@meta.data$nFeature_RNA)

FeaturePlot(g_8_6_gex, 
            features = "log_nCount_RNA", 
            pt.size = 0.2)

FeaturePlot(g_8_6_gex, 
            features = "log_nFeature_RNA", 
            pt.size = 0.2)

FeaturePlot(g_8_6_gex, 
            features = "percent.mt", 
            pt.size = 0.2)

# write out barcodes for rat astrocytes and bad quality cells
g_8_6_bad_barcodes <- g_8_6_gex@assays$SCT@counts@Dimnames[[2]][g_8_6_gex@meta.data$seurat_clusters %in% c("0", "8", "17", "20", "19", 
                                                                                                           "7", "16", "18", "10")]
write.table(g_8_6_bad_barcodes,
            file = "QC_barcodes/g_8_6_exclude_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)


StackedVlnPlot(obj = g_8_6_gex, 
               features = c("GAD1", "GAD2", "PNOC", 
                            "SLC17A6", "DLG4", "GLS", 
                            "VIM", "NES", "SOX2", 
                            "S100b", "Slc1a3", "Gfap")) +
  coord_flip()




