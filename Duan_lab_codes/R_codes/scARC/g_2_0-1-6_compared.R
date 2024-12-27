# Lexi Li 9/8/2021
# combining data for GABAergic and Glutamatergic cells

rm(list=ls())
# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)


plan("multisession", workers = 8)
# plan("sequential")
plan()
options(expressions = 20000)
options(future.globals.maxSize = 4294967295) # set future.globals.maxSize = 4GB


# read human-only filtered 10x data
human_2_0_raw <- 
  Read10X_h5(filename = "../GRCh38_mapped_only/libraries_2_0/outs/filtered_feature_bc_matrix.h5")
# select gene expression data only
human_2_0_raw_gex <- 
  CreateSeuratObject(counts = human_2_0_raw$`Gene Expression`,
                     project = "human_2_0_raw_gex")

human_2_1_raw <- 
  Read10X_h5(filename = "../GRCh38_mapped_only/libraries_2_1/outs/filtered_feature_bc_matrix.h5")
# select gene expression data only
human_2_1_raw_gex <- 
  CreateSeuratObject(counts = human_2_1_raw$`Gene Expression`,
                     project = "human_2_1_raw_gex")


human_2_6_raw <- 
  Read10X_h5(filename = "../GRCh38_mapped_only/libraries_2_6/outs/filtered_feature_bc_matrix.h5")
# select gene expression data only
human_2_6_raw_gex <- 
  CreateSeuratObject(counts = human_2_6_raw$`Gene Expression`,
                     project = "human_2_6_raw_gex")


# check the number of genes and cells
nrow(human_2_0_raw_gex) # 36601 genes (hg38 only)
ncol(human_2_0_raw_gex) # 14572 cells
nrow(human_2_1_raw_gex) # 36601 genes (hg38 only)
ncol(human_2_1_raw_gex) # 13560 cells
nrow(human_2_6_raw_gex) # 36601 genes (hg38 only)
ncol(human_2_6_raw_gex) # 14721 cells

# make a separate seurat object for subsequent operation
human_2_0_gex <- human_2_0_raw_gex
human_2_1_gex <- human_2_1_raw_gex
human_2_6_gex <- human_2_6_raw_gex


# store mitochondrial percentage in object meta data
# human starts with "MT-", rat starts with "Mt-"
human_2_0_gex <- 
  PercentageFeatureSet(human_2_0_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")
human_2_1_gex <- 
  PercentageFeatureSet(human_2_1_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")
human_2_6_gex <- 
  PercentageFeatureSet(human_2_6_gex,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")


# read CD27 barcodes
CD27_g20_glut_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/line_type_barcodes/2_0_CD27_Glut_barcodes.txt",
                                     header = F)
CD27_g20_GABA_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/line_type_barcodes/2_0_CD27_GABA_barcodes.txt",
                                     header = F)
CD27_g21_glut_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/line_type_barcodes/2_1_CD27_Glut_barcodes.txt",
                                     header = F)
CD27_g21_GABA_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/line_type_barcodes/2_1_CD27_GABA_barcodes.txt",
                                     header = F)
CD27_g26_glut_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/line_type_barcodes/2_6_CD27_Glut_barcodes.txt",
                                     header = F)
CD27_g26_GABA_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/line_type_barcodes/2_6_CD27_GABA_barcodes.txt",
                                     header = F)
CD27_g20_glut_barcodes <- unlist(CD27_g20_glut_barcodes)
CD27_g20_GABA_barcodes <- unlist(CD27_g20_GABA_barcodes)
CD27_g21_glut_barcodes <- unlist(CD27_g21_glut_barcodes)
CD27_g21_GABA_barcodes <- unlist(CD27_g21_GABA_barcodes)
CD27_g26_glut_barcodes <- unlist(CD27_g26_glut_barcodes)
CD27_g26_GABA_barcodes <- unlist(CD27_g26_GABA_barcodes)

# read CD54 barcodes
CD54_g20_glut_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/line_type_barcodes/2_0_CD54_Glut_barcodes.txt",
                                     header = F)
CD54_g20_GABA_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/line_type_barcodes/2_0_CD54_GABA_barcodes.txt",
                                     header = F)
CD54_g21_glut_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/line_type_barcodes/2_1_CD54_Glut_barcodes.txt",
                                     header = F)
CD54_g21_GABA_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/line_type_barcodes/2_1_CD54_GABA_barcodes.txt",
                                     header = F)
CD54_g26_glut_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/line_type_barcodes/2_6_CD54_Glut_barcodes.txt",
                                     header = F)
CD54_g26_GABA_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/line_type_barcodes/2_6_CD54_GABA_barcodes.txt",
                                     header = F)
CD54_g20_glut_barcodes <- unlist(CD54_g20_glut_barcodes)
CD54_g20_GABA_barcodes <- unlist(CD54_g20_GABA_barcodes)
CD54_g21_glut_barcodes <- unlist(CD54_g21_glut_barcodes)
CD54_g21_GABA_barcodes <- unlist(CD54_g21_GABA_barcodes)
CD54_g26_glut_barcodes <- unlist(CD54_g26_glut_barcodes)
CD54_g26_GABA_barcodes <- unlist(CD54_g26_GABA_barcodes)


# assign cell line and cell type information for data from every time point
# g_2_0
human_2_0_gex@meta.data$triple.ident[human_2_0_gex@assays$RNA@counts@Dimnames[[2]] %in%
                                          CD27_g20_glut_barcodes] <- "g_2_0_CD_27_glut"

human_2_0_gex@meta.data$triple.ident[human_2_0_gex@assays$RNA@counts@Dimnames[[2]] %in%
                                          CD27_g20_GABA_barcodes] <- "g_2_0_CD_27_GABA"

human_2_0_gex@meta.data$triple.ident[human_2_0_gex@assays$RNA@counts@Dimnames[[2]] %in%
                                          CD54_g20_glut_barcodes] <- "g_2_0_CD_54_glut"

human_2_0_gex@meta.data$triple.ident[human_2_0_gex@assays$RNA@counts@Dimnames[[2]] %in%
                                          CD54_g20_GABA_barcodes] <- "g_2_0_CD_54_GABA"

# g_2_1
human_2_1_gex@meta.data$triple.ident[human_2_1_gex@assays$RNA@counts@Dimnames[[2]] %in%
                                          CD27_g21_glut_barcodes] <- "g_2_1_CD_27_glut"

human_2_1_gex@meta.data$triple.ident[human_2_1_gex@assays$RNA@counts@Dimnames[[2]] %in%
                                          CD27_g21_GABA_barcodes] <- "g_2_1_CD_27_GABA"

human_2_1_gex@meta.data$triple.ident[human_2_1_gex@assays$RNA@counts@Dimnames[[2]] %in%
                                          CD54_g21_glut_barcodes] <- "g_2_1_CD_54_glut"

human_2_1_gex@meta.data$triple.ident[human_2_1_gex@assays$RNA@counts@Dimnames[[2]] %in%
                                          CD54_g21_GABA_barcodes] <- "g_2_1_CD_54_GABA"

# g_2_6
human_2_6_gex@meta.data$triple.ident[human_2_6_gex@assays$RNA@counts@Dimnames[[2]] %in%
                                          CD27_g26_glut_barcodes] <- "g_2_6_CD_27_glut"

human_2_6_gex@meta.data$triple.ident[human_2_6_gex@assays$RNA@counts@Dimnames[[2]] %in%
                                          CD27_g26_GABA_barcodes] <- "g_2_6_CD_27_GABA"

human_2_6_gex@meta.data$triple.ident[human_2_6_gex@assays$RNA@counts@Dimnames[[2]] %in%
                                          CD54_g26_glut_barcodes] <- "g_2_6_CD_54_glut"

human_2_6_gex@meta.data$triple.ident[human_2_6_gex@assays$RNA@counts@Dimnames[[2]] %in%
                                          CD54_g26_GABA_barcodes] <- "g_2_6_CD_54_GABA"


# merge three time points for CD_27 cell line
# glutamatergic cells
g20_CD_27_glut <- subset(x = human_2_0_gex, subset = triple.ident == "g_2_0_CD_27_glut")
g21_CD_27_glut <- subset(x = human_2_1_gex, subset = triple.ident == "g_2_1_CD_27_glut")
g26_CD_27_glut <- subset(x = human_2_6_gex, subset = triple.ident == "g_2_6_CD_27_glut")

group_2_CD_27_glut = merge(g20_CD_27_glut, c(g21_CD_27_glut, g26_CD_27_glut))

# GABAergic cells
g20_CD_27_GABA <- subset(x = human_2_0_gex, subset = triple.ident == "g_2_0_CD_27_GABA")
g21_CD_27_GABA <- subset(x = human_2_1_gex, subset = triple.ident == "g_2_1_CD_27_GABA")
g26_CD_27_GABA <- subset(x = human_2_6_gex, subset = triple.ident == "g_2_6_CD_27_GABA")

group_2_CD_27_GABA = merge(g20_CD_27_GABA, c(g21_CD_27_GABA, g26_CD_27_GABA))

# do the same for CD_54 cell line
# glutamatergic cells
g20_CD_54_glut <- subset(x = human_2_0_gex, subset = triple.ident == "g_2_0_CD_54_glut")
g21_CD_54_glut <- subset(x  = human_2_1_gex, subset = triple.ident == "g_2_1_CD_54_glut")
g26_CD_54_glut <- subset(x = human_2_6_gex, subset = triple.ident == "g_2_6_CD_54_glut")

group_2_CD_54_glut = merge(g20_CD_54_glut, c(g21_CD_54_glut, g26_CD_54_glut))

# GABAergic cells
g20_CD_54_GABA <- subset(x = human_2_0_gex, subset = triple.ident == "g_2_0_CD_54_GABA")
g21_CD_54_GABA <- subset(x = human_2_1_gex, subset = triple.ident == "g_2_1_CD_54_GABA")
g26_CD_54_GABA <- subset(x = human_2_6_gex, subset = triple.ident == "g_2_6_CD_54_GABA")

group_2_CD_54_GABA = merge(g20_CD_54_GABA, c(g21_CD_54_GABA, g26_CD_54_GABA))


# check common features in CD_27 and CD_54 for combining them
# glut
# common.features <- intersect(rownames(group_2_CD_27_glut), rownames(group_2_CD_54_glut))
# length(x = common.features)

# combine CD_27 and CD_54
group_2_combined_glut <- merge(group_2_CD_27_glut, y = group_2_CD_54_glut, project = "g_2_glut_combined")
group_2_combined_GABA <- merge(group_2_CD_27_GABA, y = group_2_CD_54_GABA, project = "g_2_GABA_combined")


### comparison between 0&1, 0&6

## CD_27 only

# glutamatergic

# perform SCTransformation
group_2_CD_27_glut <- 
  SCTransform(group_2_CD_27_glut, 
              vars.to.regress = "percent.mt", 
              variable.features.n = 30000,
              method = "glmGamPoi",
              verbose = T)

# test what happens if PCA and subsequent analyses are run
g2_CD27_glut <- RunPCA(group_2_CD_27_glut, verbose = FALSE)
g2_CD27_glut <- RunUMAP(g2_CD27_glut, dims = 1:30, verbose = FALSE)
g2_CD27_glut <- FindNeighbors(g2_CD27_glut, dims = 1:30, verbose = FALSE)
g2_CD27_glut <- FindClusters(g2_CD27_glut, verbose = FALSE)

# check if the three time points are well-separated
DimPlot(g2_CD27_glut,
        cols = c("darkred", "green", "black"),
        group.by = "triple.ident",
        label = F) +
  ggtitle("CD 27 glutamatergic cells at 0, 1, 6hrs")

# rename the identities to make plots with the correct and succinct labels
g2_CD27_glut@meta.data$triple.ident[g2_CD27_glut@meta.data$triple.ident == "g_2_0_CD_27_glut"] <- "0hr"
g2_CD27_glut@meta.data$triple.ident[g2_CD27_glut@meta.data$triple.ident == "g_2_1_CD_27_glut"] <- "1hr"
g2_CD27_glut@meta.data$triple.ident[g2_CD27_glut@meta.data$triple.ident == "g_2_6_CD_27_glut"] <- "6hr"
# re-subset cd27 glut data to separate into 0vs1, 0vs6 for plotting
g2_CD27_glut_01 <- subset(x = g2_CD27_glut, subset = triple.ident %in% c("0hr", "1hr"))
g2_CD27_glut_06 <- subset(x = g2_CD27_glut, subset = triple.ident %in% c("0hr", "6hr"))


# VlnPlot(g2_CD27_glut_06, 
#         features = c(),
#         group.by = "triple.ident",
#         slot = "scale.data",
#         pt.size = 0.2) + 
#   geom_boxplot(width=0.1, fill="white")

# check distribution of the data
hist(g2_CD27_glut_01@assays$SCT@scale.data[rownames(g2_CD27_glut_01@assays$SCT@scale.data) %in% "FOS", ],
     breaks = 100)  
DefaultAssay(g2_CD27_glut_01)


# write a function to plot repeatedly
fgplot_draw <- function(gene_lst, seuratobj_lst){
  count = 0
  for (sobj in seuratobj_lst){
    for(gene in gene_lst){
      count = count + 1
      p <- VlnPlot(sobj, 
                   features = gene,
                   group.by = "triple.ident",
                   slot = "scale.data",
                   pt.size = 0) + 
        geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1)
      print(p)
    }
  }
}

four_genes <- c("FOS", "NPAS4", "BDNF", "PNOC")
fgplot_draw(four_genes, c(g2_CD27_glut_01, g2_CD27_glut_06))


# GABAergic
group_2_CD_27_GABA <- 
  SCTransform(group_2_CD_27_GABA, 
              vars.to.regress = "percent.mt", 
              variable.features.n = 30000,
              method = "glmGamPoi",
              verbose = T)

g2_CD27_GABA <- RunPCA(group_2_CD_27_GABA, verbose = FALSE)
g2_CD27_GABA <- RunUMAP(g2_CD27_GABA, dims = 1:30, verbose = FALSE)
g2_CD27_GABA <- FindNeighbors(g2_CD27_GABA, dims = 1:30, verbose = FALSE)
g2_CD27_GABA <- FindClusters(g2_CD27_GABA, verbose = FALSE)

# check if the three time points are well-separated
DimPlot(g2_CD27_GABA,
        cols = c("darkred", "green", "black"),
        group.by = "triple.ident",
        label = F) +
  ggtitle("CD 27 GABAergic cells at 0, 1, 6hrs")

# rename the identities to make plots with the correct and succinct labels
g2_CD27_GABA@meta.data$triple.ident[g2_CD27_GABA@meta.data$triple.ident == "g_2_0_CD_27_GABA"] <- "0hr"
g2_CD27_GABA@meta.data$triple.ident[g2_CD27_GABA@meta.data$triple.ident == "g_2_1_CD_27_GABA"] <- "1hr"
g2_CD27_GABA@meta.data$triple.ident[g2_CD27_GABA@meta.data$triple.ident == "g_2_6_CD_27_GABA"] <- "6hr"
# re-subset cd27 GABA data to separate into 0vs1, 0vs6 for plotting
g2_CD27_GABA_01 <- subset(x = g2_CD27_GABA, subset = triple.ident %in% c("0hr", "1hr"))
g2_CD27_GABA_06 <- subset(x = g2_CD27_GABA, subset = triple.ident %in% c("0hr", "6hr"))

fgplot_draw(four_genes, c(g2_CD27_GABA_01, g2_CD27_GABA_06))


## CD_54 only

# glutamatergic
group_2_CD_54_glut <- 
  SCTransform(group_2_CD_54_glut, 
              vars.to.regress = "percent.mt", 
              variable.features.n = 30000,
              method = "glmGamPoi",
              verbose = T)

g2_CD54_glut <- RunPCA(group_2_CD_54_glut, verbose = FALSE)
g2_CD54_glut <- RunUMAP(g2_CD54_glut, dims = 1:30, verbose = FALSE)
g2_CD54_glut <- FindNeighbors(g2_CD54_glut, dims = 1:30, verbose = FALSE)
g2_CD54_glut <- FindClusters(g2_CD54_glut, verbose = FALSE)

# check if the three time points are well-separated
DimPlot(g2_CD54_glut,
        cols = c("darkred", "green", "black"),
        group.by = "triple.ident",
        label = F) +
  ggtitle("CD 54 glutamatergic cells at 0, 1, 6hrs")

# rename the identities to make plots with the correct and succinct labels
g2_CD54_glut@meta.data$triple.ident[g2_CD54_glut@meta.data$triple.ident == "g_2_0_CD_54_glut"] <- "0hr"
g2_CD54_glut@meta.data$triple.ident[g2_CD54_glut@meta.data$triple.ident == "g_2_1_CD_54_glut"] <- "1hr"
g2_CD54_glut@meta.data$triple.ident[g2_CD54_glut@meta.data$triple.ident == "g_2_6_CD_54_glut"] <- "6hr"
# re-subset cd54 glut data to separate into 0vs1, 0vs6 for plotting
g2_CD54_glut_01 <- subset(x = g2_CD54_glut, subset = triple.ident %in% c("0hr", "1hr"))
g2_CD54_glut_06 <- subset(x = g2_CD54_glut, subset = triple.ident %in% c("0hr", "6hr"))

fgplot_draw(four_genes, c(g2_CD54_glut_01, g2_CD54_glut_06))

# GABAergic
group_2_CD_54_GABA <- 
  SCTransform(group_2_CD_54_GABA, 
              vars.to.regress = "percent.mt", 
              variable.features.n = 30000,
              method = "glmGamPoi",
              verbose = T)

g2_CD54_GABA <- RunPCA(group_2_CD_54_GABA, verbose = FALSE)
g2_CD54_GABA <- RunUMAP(g2_CD54_GABA, dims = 1:30, verbose = FALSE)
g2_CD54_GABA <- FindNeighbors(g2_CD54_GABA, dims = 1:30, verbose = FALSE)
g2_CD54_GABA <- FindClusters(g2_CD54_GABA, verbose = FALSE)

# check if the three time points are well-separated
DimPlot(g2_CD54_GABA,
        cols = c("darkred", "green", "black"),
        group.by = "triple.ident",
        label = F) +
  ggtitle("CD 27 GABAergic cells at 0, 1, 6hrs")

# rename the identities to make plots with the correct and succinct labels
g2_CD54_GABA@meta.data$triple.ident[g2_CD54_GABA@meta.data$triple.ident == "g_2_0_CD_54_GABA"] <- "0hr"
g2_CD54_GABA@meta.data$triple.ident[g2_CD54_GABA@meta.data$triple.ident == "g_2_1_CD_54_GABA"] <- "1hr"
g2_CD54_GABA@meta.data$triple.ident[g2_CD54_GABA@meta.data$triple.ident == "g_2_6_CD_54_GABA"] <- "6hr"
# re-subset CD54 GABA data to separate into 0vs1, 0vs6 for plotting
g2_CD54_GABA_01 <- subset(x = g2_CD54_GABA, subset = triple.ident %in% c("0hr", "1hr"))
g2_CD54_GABA_06 <- subset(x = g2_CD54_GABA, subset = triple.ident %in% c("0hr", "6hr"))

fgplot_draw(four_genes, c(g2_CD54_GABA_01, g2_CD54_GABA_06))

## CD_27,54

# glutamatergic
group_2_combined_glut <- 
  SCTransform(group_2_combined_glut, 
              vars.to.regress = "percent.mt", 
              variable.features.n = 30000,
              method = "glmGamPoi",
              verbose = T)

g2_combined_glut <- RunPCA(group_2_combined_glut, verbose = FALSE)
g2_combined_glut <- RunUMAP(g2_combined_glut, dims = 1:30, verbose = FALSE)
g2_combined_glut <- FindNeighbors(g2_combined_glut, dims = 1:30, verbose = FALSE)
g2_combined_glut <- FindClusters(g2_combined_glut, verbose = FALSE)

# rename the identities to make plots with the correct and succinct labels
g2_combined_glut@meta.data$triple.ident[g2_combined_glut@meta.data$triple.ident == "g_2_0_CD_27_glut"] <- "0hr"
g2_combined_glut@meta.data$triple.ident[g2_combined_glut@meta.data$triple.ident == "g_2_0_CD_54_glut"] <- "0hr"
g2_combined_glut@meta.data$triple.ident[g2_combined_glut@meta.data$triple.ident == "g_2_1_CD_27_glut"] <- "1hr"
g2_combined_glut@meta.data$triple.ident[g2_combined_glut@meta.data$triple.ident == "g_2_1_CD_54_glut"] <- "1hr"
g2_combined_glut@meta.data$triple.ident[g2_combined_glut@meta.data$triple.ident == "g_2_6_CD_27_glut"] <- "6hr"
g2_combined_glut@meta.data$triple.ident[g2_combined_glut@meta.data$triple.ident == "g_2_6_CD_54_glut"] <- "6hr"
# re-subset combined glut data to separate into 0vs1, 0vs6 for plotting
g2_combined_glut_01 <- subset(x = g2_combined_glut, subset = triple.ident %in% c("0hr", "1hr"))
g2_combined_glut_06 <- subset(x = g2_combined_glut, subset = triple.ident %in% c("0hr", "6hr"))

# check if the three time points are well-separated
DimPlot(g2_combined_glut,
        cols = c("darkred", "green", "black"),
        group.by = "triple.ident",
        label = F) +
  ggtitle("CD 27+54 glutamatergic cells at 0, 1, 6hrs")

fgplot_draw(four_genes, c(g2_combined_glut_01, g2_combined_glut_06))


# GABAergic
group_2_combined_GABA <- 
  SCTransform(group_2_combined_GABA, 
              vars.to.regress = "percent.mt", 
              variable.features.n = 30000,
              method = "glmGamPoi",
              verbose = T)

g2_combined_GABA <- RunPCA(group_2_combined_GABA, verbose = FALSE)
g2_combined_GABA <- RunUMAP(g2_combined_GABA, dims = 1:30, verbose = FALSE)
g2_combined_GABA <- FindNeighbors(g2_combined_GABA, dims = 1:30, verbose = FALSE)
g2_combined_GABA <- FindClusters(g2_combined_GABA, verbose = FALSE)


# rename the identities to make plots with the correct and succinct labels
g2_combined_GABA@meta.data$triple.ident[g2_combined_GABA@meta.data$triple.ident %in% c("g_2_0_CD_27_GABA", "g_2_0_CD_54_GABA")] <- "0hr"
g2_combined_GABA@meta.data$triple.ident[g2_combined_GABA@meta.data$triple.ident %in% c("g_2_1_CD_27_GABA", "g_2_1_CD_54_GABA")] <- "1hr"
g2_combined_GABA@meta.data$triple.ident[g2_combined_GABA@meta.data$triple.ident %in% c("g_2_6_CD_27_GABA", "g_2_6_CD_54_GABA")] <- "6hr"
# re-subset combined GABA data to separate into 0vs1, 0vs6 for plotting
g2_combined_GABA_01 <- subset(x = g2_combined_GABA, subset = triple.ident %in% c("0hr", "1hr"))
g2_combined_GABA_06 <- subset(x = g2_combined_GABA, subset = triple.ident %in% c("0hr", "6hr"))

# check if the three time points are well-separated
DimPlot(g2_combined_GABA,
        cols = c("darkred", "green", "black"),
        group.by = "triple.ident",
        label = F) +
  ggtitle("CD 27+54 GABAergic cells at 0, 1, 6hrs")

fgplot_draw(four_genes, c(g2_combined_GABA_01, g2_combined_GABA_06))


# export Seurat objects
saveRDS(g2_CD27_glut_01, file = "g2_CD27_glut_01.rds")
saveRDS(g2_CD27_glut_06, file = "g2_CD27_glut_06.rds")
saveRDS(g2_CD27_GABA_01, file = "g2_CD27_GABA_01.rds")
saveRDS(g2_CD27_GABA_06, file = "g2_CD27_GABA_06.rds")
saveRDS(g2_CD54_glut_01, file = "g2_CD54_glut_01.rds")
saveRDS(g2_CD54_glut_06, file = "g2_CD54_glut_06.rds")
saveRDS(g2_CD54_GABA_01, file = "g2_CD54_GABA_01.rds")
saveRDS(g2_CD54_GABA_06, file = "g2_CD54_GABA_06.rds")

write.table(as.matrix(GetAssayData(object = g2_CD27_glut_01, slot = "scale.data")), 
            'CD27_glut_01_SCT_matrix.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)


# perform t test

