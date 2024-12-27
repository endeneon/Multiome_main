# Chuxuan 9/14/2021
# Siwei 16 Sept 2021
# Run use Rscript standalone
# Already obtained barcodes for cells to exclude, now subset human only data with the barcodes, 
#merge cell lines, and preprocess the data; finally, check marker gene expression

# rm(list=ls())
# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)
library(future)
# set threads and parallelization
plan("multisession", workers = 6)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)


# load human only data, create Seurat objects
# load data use read10x_h5
# # note this h5 file contains both atac-seq and gex information
# g_2_0_read <- 
#   Read10X_h5(filename = "../GRCh38_mapped_only/libraries_2_0/outs/filtered_feature_bc_matrix.h5")
# # load gex information first
# g_2_0_read_gex <- 
#   CreateSeuratObject(counts = g_2_0_read$`Gene Expression`,
#                      project = "g_2_0_read_gex")
# 
# g_2_1_read <- 
#   Read10X_h5(filename = "../GRCh38_mapped_only/libraries_2_1/outs/filtered_feature_bc_matrix.h5")
# # load gex information first
# g_2_1_read_gex <- 
#   CreateSeuratObject(counts = g_2_1_read$`Gene Expression`,
#                      project = "g_2_1_read_gex")
# 
# g_2_6_read <- 
#   Read10X_h5(filename = "../GRCh38_mapped_only/libraries_2_6/outs/filtered_feature_bc_matrix.h5")
# # load gex information first
# g_2_6_read_gex <- 
#   CreateSeuratObject(counts = g_2_6_read$`Gene Expression`,
#                      project = "g_2_6_read_gex")
# 
# 
# g_8_0_read <- 
#   Read10X_h5(filename = "../GRCh38_mapped_only/libraries_8_0/outs/filtered_feature_bc_matrix.h5")
# # load gex information first
# g_8_0_read_gex <- 
#   CreateSeuratObject(counts = g_8_0_read$`Gene Expression`,
#                      project = "g_8_0_read_gex")
# 
# g_8_1_read <- 
#   Read10X_h5(filename = "../GRCh38_mapped_only/libraries_8_1/outs/filtered_feature_bc_matrix.h5")
# # load gex information first
# g_8_1_read_gex <- 
#   CreateSeuratObject(counts = g_8_1_read$`Gene Expression`,
#                      project = "g_8_1_read_gex")
# 
# g_8_6_read <- 
#   Read10X_h5(filename = "../GRCh38_mapped_only/libraries_8_6/outs/filtered_feature_bc_matrix.h5")
# # load gex information first
# g_8_6_read_gex <- 
#   CreateSeuratObject(counts = g_8_6_read$`Gene Expression`,
#                      project = "g_8_6_read_gex")
# 
# # make a separate seurat object for subsequent operation
# # Note: all genes with total counts = 0 have been pre-removed
# g_2_0 <- g_2_0_read_gex
# g_2_1 <- g_2_1_read_gex
# g_2_6 <- g_2_6_read_gex
# g_8_0 <- g_8_0_read_gex
# g_8_1 <- g_8_1_read_gex
# g_8_6 <- g_8_6_read_gex
# 
# # read barcode files
# g_2_0_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_2_0_exclude_barcodes.txt",
#                              header = F)
# g_2_0_barcodes <- unlist(g_2_0_barcodes)
# g_2_1_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_2_1_exclude_barcodes.txt",
#                              header = F)
# g_2_1_barcodes <- unlist(g_2_1_barcodes)
# g_2_6_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_2_6_exclude_barcodes.txt",
#                              header = F)
# g_2_6_barcodes <- unlist(g_2_6_barcodes)
# g_8_0_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_8_0_exclude_barcodes.txt",
#                              header = F)
# g_8_0_barcodes <- unlist(g_8_0_barcodes)
# g_8_1_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_8_1_exclude_barcodes.txt",
#                              header = F)
# g_8_1_barcodes <- unlist(g_8_1_barcodes)
# g_8_6_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_8_6_exclude_barcodes.txt",
#                              header = F)
# g_8_6_barcodes <- unlist(g_8_6_barcodes)
# 
# # select the data with barcodes not in the barcode files
# # the subset() function does not take counts@dimnames, has to give cells an ident first
# g_2_0@meta.data$exclude.status <- "include"
# g_2_0@meta.data$exclude.status[colnames(g_2_0@assays$RNA@counts) %in% g_2_0_barcodes] <- "exclude"
# g_2_0 <- subset(x = g_2_0, subset = exclude.status %in% "include")
# 
# g_2_1@meta.data$exclude.status <- "include"
# g_2_1@meta.data$exclude.status[colnames(g_2_1@assays$RNA@counts) %in% g_2_1_barcodes] <- "exclude"
# g_2_1 <- subset(x = g_2_1, subset = exclude.status %in% "include")
# 
# g_2_6@meta.data$exclude.status <- "include"
# g_2_6@meta.data$exclude.status[colnames(g_2_6@assays$RNA@counts) %in% g_2_0_barcodes] <- "exclude"
# g_2_6 <- subset(x = g_2_6, subset = exclude.status %in% "include")
# 
# g_8_0@meta.data$exclude.status <- "include"
# g_8_0@meta.data$exclude.status[colnames(g_8_0@assays$RNA@counts) %in% g_8_0_barcodes] <- "exclude"
# g_8_0 <- subset(x = g_8_0, subset = exclude.status %in% "include")
# 
# g_8_1@meta.data$exclude.status <- "include"
# g_8_1@meta.data$exclude.status[colnames(g_8_1@assays$RNA@counts) %in% g_8_1_barcodes] <- "exclude"
# g_8_1 <- subset(x = g_8_1, subset = exclude.status %in% "include")
# 
# g_8_6@meta.data$exclude.status <- "include"
# g_8_6@meta.data$exclude.status[colnames(g_8_6@assays$RNA@counts) %in% g_8_6_barcodes] <- "exclude"
# g_8_6 <- subset(x = g_8_6, subset = exclude.status %in% "include")
# 
# 
# # merge the cell lines at 0, 1, 6
# time_0 = merge(g_2_0, g_8_0)
# time_1 = merge(g_2_1, g_8_1)
# time_6 = merge(g_2_6, g_8_6)
# 
# # assign time point ident
# time_0@meta.data$time.ident <- "0hr"
# time_1@meta.data$time.ident <- "1hr"
# time_6@meta.data$time.ident <- "6hr"
# 
# # merge 0, 1, 6 time points into one object
# combined_gex_read <- merge(time_0, c(time_1, time_6))
# 
# # analyze gene expression of marker genes
# # pre-treat data 
# 
# # store mitochondrial percentage in object meta data
# # human starts with "MT-", rat starts with "Mt-"
# combined_gex <- 
#   PercentageFeatureSet(combined_gex_read,
#                        pattern = c("^MT-"),
#                        col.name = "percent.mt")

# save.image(file = "combined.gex.b4.sct.RData")

load(file = "combined.gex.b4.sct.RData")

cat("!!! I changed back the object name to combined_gex !!!\n
      Check how many .RData files in the path and use the largest one!!!\n -- Siwei")
# use SCTransform()
combined_gex_raw <- combined_gex

combined_gex <- 
  SCTransform(combined_gex_raw, 
              vars.to.regress = "percent.mt",
              method = "glmGamPoi", 
              variable.features.n = 5000,
              verbose = T)
# print("!!! PCA, UMAP, etc. set")

combined_gex <- RunPCA(combined_gex,
                       verbose = T)
combined_gex <- RunUMAP(combined_gex,
                        dims = 1:30,
                        verbose = T)
combined_gex <- FindNeighbors(combined_gex,
                              dims = 1:30,
                              verbose = T)
combined_gex <- FindClusters(combined_gex,
                             verbose = T,
                             resolution = 0.7)

save.image(file = "combined.gex.after.sct.5000.RData")

# use SCTransform()
combined_gex <- 
  SCTransform(combined_gex_raw, 
              vars.to.regress = "percent.mt",
              method = "glmGamPoi", 
              variable.features.n = 10000,
              verbose = T)

combined_gex <- RunPCA(combined_gex,
                       verbose = T)
combined_gex <- RunUMAP(combined_gex,
                        dims = 1:30,
                        verbose = T)
combined_gex <- FindNeighbors(combined_gex,
                              dims = 1:30,
                              verbose = T)
combined_gex <- FindClusters(combined_gex,
                             verbose = T,
                             resolution = 0.7)

save.image(file = "combined.gex.after.sct.10000.RData")

# use SCTransform()
combined_gex <- 
  SCTransform(combined_gex_raw, 
              vars.to.regress = "percent.mt",
              method = "glmGamPoi", 
              variable.features.n = 20000,
              verbose = T)

combined_gex <- RunPCA(combined_gex,
                       verbose = T)
combined_gex <- RunUMAP(combined_gex,
                        dims = 1:30,
                        verbose = T)
combined_gex <- FindNeighbors(combined_gex,
                              dims = 1:30,
                              verbose = T)
combined_gex <- FindClusters(combined_gex,
                             verbose = T,
                             resolution = 0.7)

save.image(file = "combined.gex.after.sct.20000.RData")
