# 28 Jul 2020 Siwei
# Run RNASeq use SID0 and observe SP4 expression related to cluster

# init
library(Seurat)
library(future)
library(readr)
library(sctransform)
library(stringr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(hdf5r)

# set parameters
options(future.globals.maxSize = 6442450944)
options(future.fork.enable = T) # force enable forks


# setup multithread
plan("multicore", workers = 6) # force multicore

# load corresponding RNASeq data
SID0.RNA.object <- readRDS(file = "SID0.RNASeq.RData")
SID0.RNA.object$tech <- "rna"
DefaultAssay(SID0.RNA.object)

# recalculate dimensions
SID0.RNA.object <- RunPCA(SID0.RNA.object, verbose = T)
SID0.RNA.object <- RunUMAP(SID0.RNA.object, 
                           dims = 1:50, 
                           verbose = T)
SID0.RNA.object <- FindNeighbors(SID0.RNA.object, 
                                 dims = 1:50,
                                 verbose = T)
SID0.RNA.object <- FindClusters(SID0.RNA.object,
                                resolution = 0.39,
                                verbose = T)

DimPlot(SID0.RNA.object, label = T) +
  NoLegend()
Idents(SID0.RNA.object) <- SID0.RNA.object$seurat_clusters
DimPlot(SID0.RNA.object, 
        cols = SID0.RNA.object$orig.ident,
        label = T) +
  NoLegend()

FeaturePlot(SID0.RNA.object, 
            features = "SP4")
FeaturePlot(SID0.RNA.object, 
            features = "VPS45")
FeaturePlot(SID0.RNA.object, 
            features = "NEUROD1")


FeaturePlot(SID0.RNA.object, 
            features = "SLC17A6")
FeaturePlot(SID0.RNA.object, 
            features = "SLC17A7")

FeaturePlot(SID0.RNA.object, 
            features = "NES")

# assign SP4 (+/-) identity
hist(SID0.RNA.object@assays$SCT@data[rownames(SID0.RNA.object@assays$SCT@data) %in% "SP4", ],
     breaks = 50)
sum(SID0.RNA.object@assays$SCT@data[rownames(SID0.RNA.object@assays$SCT@data) %in% "SP4", ] > 0)
SP4.exp <- SID0.RNA.object@assays$SCT@data[rownames(SID0.RNA.object@assays$SCT@data) %in% "SP4", ]
SP4.exp <- SP4.exp[SP4.exp > 0]
hist(SP4.exp, breaks = 50)
# sum(rownames(SID0.RNA.object@assays$SCT@data) == "SP4")

SID0.RNA.object$SP4_exp <- "negative"
SID0.RNA.object$SP4_exp[SID0.RNA.object@assays$SCT@data[rownames(SID0.RNA.object@assays$SCT@data) %in% "SP4", ] > 0] <- "positive"
Idents(SID0.RNA.object) <- SID0.RNA.object$SP4_exp

# subset cell types
Glut.SID0.RNA.object <- SID0.RNA.object[, SID0.RNA.object$seurat_clusters %in% c("1", "2", "8")]
GABA.SID0.RNA.object <- SID0.RNA.object[, SID0.RNA.object$seurat_clusters %in% c("0", "6")]
NPC.SID0.RNA.object <- SID0.RNA.object[, SID0.RNA.object$seurat_clusters %in% c("3", "4")]

# DE analysis in SP4 (+/-) cells
DefaultAssay(Glut.SID0.RNA.object)
Glut.SP4.markers <- FindMarkers(Glut.SID0.RNA.object,
                                ident.1 = "positive", 
                                ident.2 = "negative",
                                min.pct = 0.1,
                                logfc.threshold = 0.1,
                                assay = "RNA",
                                test.use = "DESeq2",
                                verbose = T)
GABA.SP4.markers <- FindMarkers(GABA.SID0.RNA.object,
                                ident.1 = "positive", 
                                ident.2 = "negative",
                                min.pct = 0.02,
                                logfc.threshold = 0.1,
                                test.use = "bimod",
                                verbose = T)
NPC.SP4.markers <- FindMarkers(NPC.SID0.RNA.object,
                               ident.1 = "positive", 
                               ident.2 = "negative",
                               min.pct = 0.02,
                               logfc.threshold = 0.1,
                               test.use = "bimod",
                               verbose = T)
