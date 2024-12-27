# Chuxuan 9/7/2021
# Using GEX data of g_2_6 to generate barcodes by line and cell type

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
options(future.globals.maxSize = 4294967295) # set future.globals.maxSize = 4GB

# load data
# load data use read10x_h5
# note this h5 file contains both atac-seq and gex information
g_2_6_raw <- 
  Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_2_6/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_2_6_raw_gex <- 
  CreateSeuratObject(counts = g_2_6_raw$`Gene Expression`,
                     project = "g_2_6_raw_gex")

# check the number of genes and cells
nrow(g_2_6_raw_gex) # 91461 genes (hg38 + Rnor6)
ncol(g_2_6_raw_gex) # 16370 cells

# make a seurat object for subsequent operation
# Note: all genes with total counts = 0 have been pre-removed
g_2_6_gex <- g_2_6_raw_gex

# store mitochondrial percentage in object meta data
# human starts with "MT-", rat starts with "Mt-"

g_2_6_gex <- 
  PercentageFeatureSet(g_2_6_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")

VlnPlot(g_2_6_gex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# pre-treat data 
# use SCTransform()
# Note: SCTransform has to be downloaded and installed directly from CRAN
# since it has to be compiled locally, do not use the anaconda distribution
# including NormalizeData(), ScaleData(), and FindVariableFeatures
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

# make dimension plot
DimPlot(g_2_6_gex, 
        label = T) +
  NoLegend()

# check specific genes
FeaturePlot(g_2_6_gex, 
            features = c("GAD1", "SLC17A6", 
                         "Gfap", "SOX2"), 
            pt.size = 0.2,
            ncol = 2)

# assign cell line identities
## add one column of cell line identities to meta.data
g_2_6_gex@meta.data$cell.line.ident <- "unknown"

## import cell ident barcodes
CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_6_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
CD_27_barcodes <- unlist(CD_27_barcodes)

# read CD 27 barcodes from an ATAC-Seq-only list made from the following bash command:
# cat g_2_6_atac_CD27_CD54.best | cut -f 1-6 | grep "SNG" | grep "CD_27" | cut -f 1 > g_2_6_atac_seq_cd_27.barcodes
## !! code added !!
# CD_27_barcodes <- 
#   read_csv("/home/cli/NVME/scARC_Duan_018/hg38_Rnor6_mixed/group_2/atac_output/g_2_6_atac_seq_cd_27.barcodes", 
#            col_names = FALSE)
# CD_27_barcodes <- unlist(CD_27_barcodes)
# 
# sum(g_2_6_gex@assays$SCT@counts@Dimnames[[2]] %in%
#       CD_27_barcodes)
## !! code added ends here !!

CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_6_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
CD_54_barcodes <- unlist(CD_54_barcodes)

sum(g_2_6_gex@assays$SCT@counts@Dimnames[[2]] %in%
      CD_54_barcodes)

## assign cell line identity
g_2_6_gex@meta.data$cell.line.ident[g_2_6_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      CD_27_barcodes] <- "CD_27"
g_2_6_gex@meta.data$cell.line.ident[g_2_6_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      CD_54_barcodes] <- "CD_54"

g_2_6_gex@meta.data$cell.line.ident <-
  factor(g_2_6_gex@meta.data$cell.line.ident)

summary(g_2_6_gex@meta.data$cell.line.ident)

# summary(g_2_6_gex@meta.data$cell.line.ident <-
#           factor(g_2_6_gex@meta.data$cell.line.ident))


StackedVlnPlot(obj = g_2_6_gex, 
               features = c("GAD1", "GAD2", "PNOC", 
                            "SLC17A6", "DLG4", "GLS", 
                            "VIM", "NES", "SOX2", 
                            "S100b", "Slc1a3", "Gfap")) +
  coord_flip()

## separate cell types
Glut_cell_barcodes <- 
  g_2_6_gex@assays$SCT@counts@Dimnames[[2]][g_2_6_gex@meta.data$seurat_clusters %in% 
                                              c("0", "9", "10", "12", "14", "22", "4", "5", "6", "18", "23")]
GABA_cell_barcodes <- 
  g_2_6_gex@assays$SCT@counts@Dimnames[[2]][g_2_6_gex@meta.data$seurat_clusters %in% 
                                              c("2", "7", "13", "15", "17")]

## write out line-type barcodes
CD_27_Glut_barcodes <- CD_27_barcodes[CD_27_barcodes %in% Glut_cell_barcodes]
write.table(CD_27_Glut_barcodes,
            file = "line_type_barcodes/2_6_CD27_Glut_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

CD_27_GABA_barcodes <- CD_27_barcodes[CD_27_barcodes %in% GABA_cell_barcodes]
write.table(CD_27_GABA_barcodes,
            file = "line_type_barcodes/2_6_CD27_GABA_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

CD_54_Glut_barcodes <- CD_54_barcodes[CD_54_barcodes %in% Glut_cell_barcodes]
write.table(CD_54_Glut_barcodes,
            file = "line_type_barcodes/2_6_CD54_Glut_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

CD_54_GABA_barcodes <- CD_54_barcodes[CD_54_barcodes %in% GABA_cell_barcodes]
write.table(CD_54_GABA_barcodes,
            file = "line_type_barcodes/2_6_CD54_GABA_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

# check cell distribution
# check specific genes
DimPlot(g_2_6_gex, 
        label = T) +
  NoLegend()


FeaturePlot(g_2_6_gex,
            features = c("GAD1", "SLC17A7",
                         "Gfap", "SOX2"),
            pt.size = 0.2,
            ncol = 2)

FeaturePlot(g_2_6_gex, 
            features = c("GLS"), 
            pt.size = 0.2,
            ncol = 1)

FeaturePlot(g_2_6_gex, 
            features = c("GAD1", "GAD2"), 
            # split.by = "cell.line.ident",
            pt.size = 0.2,
            blend = T,
            ncol = 1)

FeaturePlot(g_2_6_gex, 
            features = "SLC17A6", 
            # split.by = "cell.line.ident",
            pt.size = 0.2,
            ncol = 1)
FeaturePlot(g_2_6_gex, 
            features = c("SLC17A6", "SLC17A7"), 
            # split.by = "cell.line.ident",
            blend = T,
            pt.size = 0.2,
            ncol = 1)



DimPlot(g_2_6_gex,
        cols = c("darkred", "green", "black"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "cell.line.ident",
        label = F) +
  ggtitle("group_2_6hr")# +
#NoLegend()

# head(FindMarkers(g_2_6_gex,
#                  ident.1 = c("0", "9", "10", "12", "14", "22", "4", "5", "6", "18", "23"),
#                  ident.2 = c("2", "7", "13", "15", "17")))


# check the number of cells in each type of cells 
sum(g_2_6_gex@meta.data$seurat_clusters %in% 
      c(0, 9, 10, 12, 14, 22, 4, 5, 6, 18, 23))  # glutamatergic cells (SLC17A6)

sum(g_2_6_gex@meta.data$seurat_clusters %in% 
      c(2, 7, 13, 15, 17)) # GABAergic cells (GAD1)

sum(g_2_6_gex@meta.data$seurat_clusters %in% 
      c(1, 3, 8, 11))  # rat astrocytes (Gfap)

sum(g_2_6_gex@meta.data$seurat_clusters %in% 
      c(16, 19, 20)) # proliferating progenitor cells (SOX2)

sum((g_2_6_gex@meta.data$seurat_clusters %in% 
       c(0, 9, 10, 12, 14, 22, 4, 5, 6, 18, 23, 2, 7, 13, 15, 17)) & 
      is.na(g_2_6_gex@meta.data$cell.line.ident)) # all human cells that are not
                                                  # assigned any cell line

sum(g_2_6_gex@meta.data$seurat_clusters %in% 
      c(0, 9, 10, 12, 14, 22, 4, 5, 6, 18, 23, 2, 7, 13, 15, 17)) # all human cells 

