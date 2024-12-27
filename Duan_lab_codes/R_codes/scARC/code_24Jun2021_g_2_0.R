# Siwei 24 Jun 2021
# import GEX data of g_2_0 (2 samples) for t-sne plotting

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

nrow(g_2_0_raw_gex) # 91461 genes (hg38 + Rnor6)
ncol(g_2_0_raw_gex) # 16068 cells

# make a seurat object for subsequent operation
# Note: all genes with total counts = 0 have been pre-removed
g_2_0_gex <- g_2_0_raw_gex

# store mitochondrial percentage in object meta data
# human starts with "MT-", rat starts with "Mt-"
# g_2_0_gex <- 
#   PercentageFeatureSet(g_2_0_gex,
#                        pattern = c("^M[Tt]-"),
#                        col.name = "percent.mt")

g_2_0_gex <- 
  PercentageFeatureSet(g_2_0_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")

# pre-treat data 
# use SCTransform()
# Note: SCTransform has to be downloaded and installed directly from CRAN
# since it has to be compiled locally, do not use the anaconda distribution
# including NormalizeData(), ScaleData(), and FindVariableFeatures
# just compatibility issue with RFast
g_2_0_gex <- 
  SCTransform(g_2_0_gex, 
              vars.to.regress = "percent.mt",
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

# make dimension plot
DimPlot(g_2_0_gex, 
        label = T) +
  NoLegend()

# check specific genes
FeaturePlot(g_2_0_gex, 
            features = c("GAD1", "SLC17A6", 
                         "Gfap", "SOX2"), 
            pt.size = 0.2,
            ncol = 2)

# assign cell line identities
## add one column of cell line identities to meta.data
g_2_0_gex@meta.data$cell.line.ident <- "unknown"

## import cell ident barcodes
CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/backup/g_2_0_CD27_barcodes.txt", 
           col_names = FALSE)
CD_27_barcodes <- unlist(CD_27_barcodes)

CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/backup/g_2_0_CD54_barcodes.txt", 
           col_names = FALSE)
CD_54_barcodes <- unlist(CD_54_barcodes)


## assign cell line identity
g_2_0_gex@meta.data$cell.line.ident[g_2_0_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      CD_27_barcodes] <- "CD_27"
g_2_0_gex@meta.data$cell.line.ident[g_2_0_gex@assays$SCT@counts@Dimnames[[2]] %in%
                                      CD_54_barcodes] <- "CD_54"

g_2_0_gex@meta.data$cell.line.ident <-
  factor(g_2_0_gex@meta.data$cell.line.ident)

summary(g_2_0_gex@meta.data$cell.line.ident)

# summary(g_2_0_gex@meta.data$cell.line.ident <-
#           factor(g_2_0_gex@meta.data$cell.line.ident))

# create stacked violin plots to check clusters' cell type identity
StackedVlnPlot(obj = g_2_0_gex, 
               features = c("GAD1", "GAD2", "PNOC", 
                            "SLC17A6", "DLG4", "GLS", 
                            "VIM", "NES", "SOX2", 
                            "S100b", "Slc1a3", "Gfap")) +
  coord_flip()

## separate cell types
Glut_cell_barcodes <- 
  g_2_0_gex@assays$SCT@counts@Dimnames[[2]][g_2_0_gex@meta.data$seurat_clusters %in% 
                                              c("0", "1", "6", "7", "12", "19")]
GABA_cell_barcodes <- 
  g_2_0_gex@assays$SCT@counts@Dimnames[[2]][g_2_0_gex@meta.data$seurat_clusters %in% 
                                              c("3", "4", "9", "18")]

## write out line-type barcodes
CD_27_Glut_barcodes <- CD_27_barcodes[CD_27_barcodes %in% Glut_cell_barcodes]
write.table(CD_27_Glut_barcodes,
            file = "line_type_barcodes/2_0_CD27_Glut_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

CD_27_GABA_barcodes <- CD_27_barcodes[CD_27_barcodes %in% GABA_cell_barcodes]
write.table(CD_27_GABA_barcodes,
            file = "line_type_barcodes/2_0_CD27_GABA_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

CD_54_Glut_barcodes <- CD_54_barcodes[CD_54_barcodes %in% Glut_cell_barcodes]
write.table(CD_54_Glut_barcodes,
            file = "line_type_barcodes/2_0_CD54_Glut_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

CD_54_GABA_barcodes <- CD_54_barcodes[CD_54_barcodes %in% GABA_cell_barcodes]
write.table(CD_54_GABA_barcodes,
            file = "line_type_barcodes/2_0_CD54_GABA_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

# check cell distribution
# check specific genes
DimPlot(g_2_0_gex, 
        label = T) +
  NoLegend()

FeaturePlot(g_2_0_gex, 
            features = c("GAD1", "SLC17A7", 
                         "Gfap", "SOX2"), 
            pt.size = 0.2,
            ncol = 2) 

FeaturePlot(g_2_0_gex, 
            features = c("GLS"), 
            pt.size = 0.2,
            ncol = 1)

FeaturePlot(g_2_0_gex, 
            features = "GAD1", 
            # split.by = "cell.line.ident",
            pt.size = 0.2,
            ncol = 1)
FeaturePlot(g_2_0_gex, 
            features = "SLC17A6", 
            # split.by = "cell.line.ident",
            pt.size = 0.2,
            ncol = 1)
FeaturePlot(g_2_0_gex, 
            features = c("SLC17A6", "SLC17A7"), 
            # split.by = "cell.line.ident",
            blend = T,
            pt.size = 0.2,
            ncol = 1)



DimPlot(g_2_0_gex,
        cols = c("darkred", "green", "black"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_0_gex@meta.data$cell.line.ident,
        group.by = "cell.line.ident",
        label = F) +
  ggtitle("group_2_0hr")# +
NoLegend()

head(FindMarkers(g_2_0_gex,
                 ident.1 = c("0", "6", "19", "1", "7", "12"),
                 ident.2 = c("10", "14")))


sum(g_2_0_gex@meta.data$seurat_clusters %in% 
      c(21, 9, 0, 10, 13, 14, 8, 4, 7))
sum(g_2_0_gex@meta.data$seurat_clusters %in% 
      c(2, 6, 5, 20))
sum(g_2_0_gex@meta.data$seurat_clusters %in% 
      c(2, 11, 5))
sum(g_2_0_gex@meta.data$seurat_clusters %in% 
      c(10, 14))

sum((g_2_0_gex@meta.data$seurat_clusters %in% 
      c(18, 3, 4, 9, 6, 0, 19, 17, 1, 7, 12, 20)) & 
      is.na(g_2_0_gex@meta.data$cell.line.ident))

sum(g_2_0_gex@meta.data$seurat_clusters %in% 
       c(18, 3, 4, 9, 6, 0, 19, 17, 1, 7, 12, 20))


