# 04 Aug 2020 Siwei
# Run RNASeq use SID1 and observe SP4 expression related to cluster

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
SID1.RNA.object <- readRDS(file = "SID1.RNA.object.RData")
SID1.RNA.object$tech <- "rna"
DefaultAssay(SID1.RNA.object)

# recalculate dimensions
SID1.RNA.object <- RunPCA(SID1.RNA.object, verbose = T)
SID1.RNA.object <- RunUMAP(SID1.RNA.object, 
                           dims = 1:50, 
                           verbose = T)
SID1.RNA.object <- FindNeighbors(SID1.RNA.object, 
                                 dims = 1:50,
                                 verbose = T)
SID1.RNA.object <- FindClusters(SID1.RNA.object,
                                resolution = 0.394,
                                verbose = T)

DimPlot(SID1.RNA.object, label = T) +
  NoLegend()
Idents(SID1.RNA.object) <- SID1.RNA.object$seurat_clusters
DimPlot(SID1.RNA.object, 
        cols = SID1.RNA.object$orig.ident,
        label = T) +
  NoLegend()

FeaturePlot(SID1.RNA.object, 
            features = "SP4")
FeaturePlot(SID1.RNA.object, 
            features = "VPS45")
FeaturePlot(SID1.RNA.object, 
            features = "NEUROD1")


FeaturePlot(SID1.RNA.object, 
            features = "SLC17A6")
FeaturePlot(SID1.RNA.object, 
            features = "SLC17A7")

FeaturePlot(SID1.RNA.object, 
            features = "NES")

# assign SP4 (+/-) identity
hist(SID1.RNA.object@assays$SCT@data[rownames(SID1.RNA.object@assays$SCT@data) %in% "SP4", ],
     breaks = 50)
sum(SID1.RNA.object@assays$SCT@data[rownames(SID1.RNA.object@assays$SCT@data) %in% "SP4", ] > 0)
SP4.exp <- SID1.RNA.object@assays$SCT@data[rownames(SID1.RNA.object@assays$SCT@data) %in% "SP4", ]
SP4.exp <- SP4.exp[SP4.exp > 0]
hist(SP4.exp, breaks = 50)
# sum(rownames(SID1.RNA.object@assays$SCT@data) == "SP4")

SID1.RNA.object$SP4_exp <- "negative"
SID1.RNA.object$SP4_exp[SID1.RNA.object@assays$SCT@data[rownames(SID1.RNA.object@assays$SCT@data) %in% "SP4", ] > 0] <- "positive"
Idents(SID1.RNA.object) <- SID1.RNA.object$SP4_exp

# subset cell types
Glut.SID1.RNA.object <- SID1.RNA.object[, SID1.RNA.object$seurat_clusters %in% c("1", "2", "8")]
GABA.SID1.RNA.object <- SID1.RNA.object[, SID1.RNA.object$seurat_clusters %in% c("0", "7")]
NPC.SID1.RNA.object <- SID1.RNA.object[, SID1.RNA.object$seurat_clusters %in% c("3", "4")]

# DE analysis in SP4 (+/-) cells
DefaultAssay(Glut.SID1.RNA.object)
Glut.SP4.markers <- FindMarkers(Glut.SID1.RNA.object,
                                ident.1 = "positive", 
                                ident.2 = "negative",
                                min.pct = 0.02,
                                logfc.threshold = 0.1,
                                # assay = "RNA",
                                test.use = "bimod",
                                verbose = T)
GABA.SP4.markers <- FindMarkers(GABA.SID1.RNA.object,
                                ident.1 = "positive", 
                                ident.2 = "negative",
                                min.pct = 0.02,
                                logfc.threshold = 0.1,
                                test.use = "bimod",
                                verbose = T)
NPC.SP4.markers <- FindMarkers(NPC.SID1.RNA.object,
                               ident.1 = "positive", 
                               ident.2 = "negative",
                               min.pct = 0.02,
                               logfc.threshold = 0.1,
                               test.use = "bimod",
                               verbose = T)

write.table(Glut.SP4.markers, 
            file = "SID1.Glut.SP4.markers.txt",
            sep = "\t")
write.table(GABA.SP4.markers,
            file = "SID1.GABA.SP4.markers.txt",
            sep = "\t")
write.table(NPC.SP4.markers,
            file = "SID1.NPC.SP4.markers.txt",
            sep = "\t")
