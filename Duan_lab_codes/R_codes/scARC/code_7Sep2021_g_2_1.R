# Chuxuan 9/7/2021
# Using GEX data of g_2_1 to generate barcodes by line and cell type

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
g_2_1_raw <- 
  Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_2_1/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_2_1_raw_gex <- 
  CreateSeuratObject(counts = g_2_1_raw$`Gene Expression`,
                     project = "g_2_1_raw_gex")

nrow(g_2_1_raw_gex) # 91461 genes (hg38 + Rnor6)
ncol(g_2_1_raw_gex) # 16068 cells

# make a seurat object for subsequent operation
# Note: all genes with total counts = 0 have been pre-removed
g_2_1_gex <- g_2_1_raw_gex

# store mitochondrial percentage in object meta data
# human starts with "MT-", rat starts with "Mt-"
# create separate Seurat objects to test the difference between including
# and not including rat astrocytes, check which one is better

g_2_1_gex_wo_rat <- 
  PercentageFeatureSet(g_2_1_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")

# g_2_1_gex_w_rat <- 
#   PercentageFeatureSet(g_2_1_gex,
#                        pattern = c("^M[Tt]-"),
#                        col.name = "percent.mt")

VlnPlot(g_2_1_gex_wo_rat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(g_2_1_gex_w_rat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# pre-treat data 
# use SCTransform()
# Note: SCTransform has to be downloaded and installed directly from CRAN
# since it has to be compiled locally, do not use the anaconda distribution
# including NormalizeData(), ScaleData(), and FindVariableFeatures
g_2_1_gex_wo_rat <- 
  SCTransform(g_2_1_gex_wo_rat, 
              vars.to.regress = "percent.mt",
              method = "glmGamPoi",
              verbose = T)

g_2_1_gex_wo_rat <- RunPCA(g_2_1_gex_wo_rat, 
                    verbose = T)
g_2_1_gex_wo_rat <- RunUMAP(g_2_1_gex_wo_rat,
                     dims = 1:30,
                     verbose = T)
g_2_1_gex_wo_rat <- FindNeighbors(g_2_1_gex_wo_rat,
                           dims = 1:30,
                           verbose = T)
g_2_1_gex_wo_rat <- FindClusters(g_2_1_gex_wo_rat,
                          verbose = T)


# g_2_1_gex_w_rat <- 
#   SCTransform(g_2_1_gex_w_rat, 
#               vars.to.regress = "percent.mt",
#               method = "glmGamPoi",
#               verbose = T)
# 
# g_2_1_gex_w_rat <- RunPCA(g_2_1_gex_w_rat, 
#                            verbose = T)
# g_2_1_gex_w_rat <- RunUMAP(g_2_1_gex_w_rat,
#                             dims = 1:30,
#                             verbose = T)
# g_2_1_gex_w_rat <- FindNeighbors(g_2_1_gex_w_rat,
#                                   dims = 1:30,
#                                   verbose = T)
# g_2_1_gex_w_rat <- FindClusters(g_2_1_gex_w_rat,
#                                  verbose = T)

# make dimension plot
DimPlot(g_2_1_gex_wo_rat, 
        label = T) +
  NoLegend()

# DimPlot(g_2_1_gex_w_rat, 
#         label = T) +
#   NoLegend()
# check specific genes
FeaturePlot(g_2_1_gex_wo_rat, 
            features = c("GAD1", "SLC17A6", 
                         "Gfap", "SOX2"), 
            pt.size = 0.2,
            ncol = 2)

g_2_1_gex_wo_rat@meta.data$log_nCount_RNA <- log(g_2_1_gex_wo_rat@meta.data$nCount_RNA)
g_2_1_gex_wo_rat@meta.data$log_nFeature_RNA <- log(g_2_1_gex_wo_rat@meta.data$nFeature_RNA)


FeaturePlot(g_2_1_gex_wo_rat, 
            features = "log_nCount_RNA", 
            pt.size = 0.2)

FeaturePlot(g_2_1_gex_wo_rat, 
            features = "log_nFeature_RNA", 
            pt.size = 0.2)

# FeaturePlot(g_2_1_gex_w_rat, 
#             features = c("GAD1", "SLC17A6", 
#                          "Gfap", "SOX2"), 
#             pt.size = 0.2,
#             ncol = 2)
# for consistency, use transformed data excluding rat from this point on

# assign cell line identities
## add one column of cell line identities to meta.data
g_2_1_gex_wo_rat@meta.data$cell.line.ident <- "unknown"

## import cell ident barcodes
CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_1_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
CD_27_barcodes <- unlist(CD_27_barcodes)

CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_1_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
CD_54_barcodes <- unlist(CD_54_barcodes)


## assign cell line identity
g_2_1_gex_wo_rat@meta.data$cell.line.ident[g_2_1_gex_wo_rat@assays$SCT@counts@Dimnames[[2]] %in%
                                      CD_27_barcodes] <- "CD_27"
g_2_1_gex_wo_rat@meta.data$cell.line.ident[g_2_1_gex_wo_rat@assays$SCT@counts@Dimnames[[2]] %in%
                                      CD_54_barcodes] <- "CD_54"

g_2_1_gex_wo_rat@meta.data$cell.line.ident <-
  factor(g_2_1_gex_wo_rat@meta.data$cell.line.ident)

summary(g_2_1_gex_wo_rat@meta.data$cell.line.ident)

# summary(g_2_1_gex_wo_rat@meta.data$cell.line.ident <-
#           factor(g_2_1_gex_wo_rat@meta.data$cell.line.ident))


StackedVlnPlot(obj = g_2_1_gex_wo_rat, 
               features = c("GAD1", "GAD2", "PNOC", 
                            "SLC17A6", "SLC17A7", "DLG4", "GLS", 
                            "VIM", "NES", "SOX2", 
                            "S100b", "Slc1a3", "Gfap")) +
  coord_flip()


## separate cell types
Glut_cell_barcodes <- # these clusters are glutamatergic cells
  g_2_1_gex_wo_rat@assays$SCT@counts@Dimnames[[2]][g_2_1_gex_wo_rat@meta.data$seurat_clusters %in% 
                                              c("0", "2", "6", "8", "11", "13", "18", "21")]
GABA_cell_barcodes <- # these clusters are GABAergic cells
  g_2_1_gex_wo_rat@assays$SCT@counts@Dimnames[[2]][g_2_1_gex_wo_rat@meta.data$seurat_clusters %in% 
                                              c("1", "3", "15", "16", "22")]

## write out line-type barcodes
CD_27_Glut_barcodes <- CD_27_barcodes[CD_27_barcodes %in% Glut_cell_barcodes]
write.table(CD_27_Glut_barcodes,
            file = "line_type_barcodes/2_1_CD27_Glut_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

CD_27_GABA_barcodes <- CD_27_barcodes[CD_27_barcodes %in% GABA_cell_barcodes]
write.table(CD_27_GABA_barcodes,
            file = "line_type_barcodes/2_1_CD27_GABA_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

CD_54_Glut_barcodes <- CD_54_barcodes[CD_54_barcodes %in% Glut_cell_barcodes]
write.table(CD_54_Glut_barcodes,
            file = "line_type_barcodes/2_1_CD54_Glut_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

CD_54_GABA_barcodes <- CD_54_barcodes[CD_54_barcodes %in% GABA_cell_barcodes]
write.table(CD_54_GABA_barcodes,
            file = "line_type_barcodes/2_1_CD54_GABA_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

# check cell distribution
# check specific genes
DimPlot(g_2_1_gex_wo_rat, 
        label = T) +
  NoLegend()

FeaturePlot(g_2_1_gex_wo_rat, 
            features = c("GAD1", "SLC17A6", 
                         "Gfap", "SOX2"), 
            pt.size = 0.2,
            ncol = 2) 

FeaturePlot(g_2_1_gex_wo_rat, 
            features = c("GLS"), 
            pt.size = 0.2,
            ncol = 1)

FeaturePlot(g_2_1_gex_wo_rat, 
            features = c("GAD1", "GAD2"), # check markers for GABA
            # split.by = "cell.line.ident",
            pt.size = 0.2,
            blend = T,
            ncol = 1)

FeaturePlot(g_2_1_gex_wo_rat, 
            features = c("SLC17A6", "SLC17A7"), # check markers for Glut
            # split.by = "cell.line.ident",
            blend = T,
            pt.size = 0.2,
            ncol = 1)



DimPlot(g_2_1_gex_wo_rat,
        cols = c("darkred", "green", "black"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_1_gex_wo_rat@meta.data$cell.line.ident,
        group.by = "cell.line.ident",
        label = F) +
  ggtitle("group_2_1hr")

head(FindMarkers(g_2_1_gex_wo_rat,
                 ident.1 = c("1", "3", "15", "16", "22"),
                 ident.2 = c("0", "2", "6", "8", "11", "13", "18", "21")))


# check the number of cells in each type of cells 
sum(g_2_1_gex_wo_rat@meta.data$seurat_clusters %in% 
      c(1, 3, 15, 16, 22))  # GABAergic cells

sum(g_2_1_gex_wo_rat@meta.data$seurat_clusters %in% 
      c(0, 2, 6, 8, 11, 13, 18, 21)) # glutamatergic cells

sum(g_2_1_gex_wo_rat@meta.data$seurat_clusters %in% 
      c(4, 10, 5))  # rat astrocytes

sum(g_2_1_gex_wo_rat@meta.data$seurat_clusters %in% 
      c(10, 14)) # proliferating progenitor cells

sum((g_2_1_gex_wo_rat@meta.data$seurat_clusters %in% 
       c(1, 3, 15, 16, 22, 0, 2, 6, 8, 11, 13, 18, 21)) & 
      is.na(g_2_1_gex@meta.data$cell.line.ident)) # all human cells that are not
                                                  # assigned any cell line

sum(g_2_1_gex_wo_rat@meta.data$seurat_clusters %in% 
      c(1, 3, 15, 16, 22, 0, 2, 6, 8, 11, 13, 18, 21)) # all human cells 
