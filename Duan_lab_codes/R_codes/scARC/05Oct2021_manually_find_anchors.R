# Chuxuan Li 10/05/2021
# Find anchor genes from time = 0hr data by looking up the top markers, then use
# these markers to check if the clusters can be linked. 

rm(list=ls())
# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)
library(future)
library(RColorBrewer)
library(dplyr)

# set threads and parallelization
plan("multisession", workers = 3)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

load("raw_group_2_group_8_seurat_objects.RData")

# merge and get the objects for each time point
time_0 <- merge(g_2_0, g_8_0)
time_1 <- merge(g_2_1, g_8_1)
time_6 <- merge(g_2_6, g_8_6)

# use this function to do the steps before clustering
process_data <- function(seurat_obj){
  # store mitochondrial percentage in object meta data
  # human starts with "MT-", rat starts with "Mt-"
  seurat_obj <- PercentageFeatureSet(seurat_obj,
                                     pattern = c("^MT-"),
                                     col.name = "percent.mt",
                                     assay = 'RNA')
  #use SCTransform()
  seurat_obj <-
    SCTransform(seurat_obj,
                vars.to.regress = "percent.mt",
                method = "glmGamPoi",
                variable.features.n = 8000,
                verbose = T,
                seed.use = 42)
  seurat_obj <- RunPCA(seurat_obj,
                       seed.use = 42,
                       verbose = T)
  seurat_obj <- RunUMAP(seurat_obj,
                        dims = 1:30,
                        seed.use = 42,
                        verbose = T)
  seurat_obj <- FindNeighbors(seurat_obj,
                              dims = 1:30,
                              verbose = T)
  seurat_obj <- FindClusters(seurat_obj,
                             verbose = T, 
                             random.seed = 42,
                             resolution = 0.14)
  
  return(seurat_obj)
}

time_0_p <- process_data(time_0)

# adjust number of clusters
time_0_p <- FindClusters(time_0_p,
                         verbose = T, 
                         random.seed = 42,
                         resolution = 0.12)
# take a look at time 0 
DimPlot(time_0_p, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of merged 0hr")

time_0_markers <- FindAllMarkers(time_0_p)
top_markers <- time_0_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.table(top_markers, sep = "," , file = "top_markers_0hr.csv")
cluster_0_top <- top_markers[1:10, ]$gene
cluster_1_top <- top_markers[11:20, ]$gene
cluster_2_top <- top_markers[21:30, ]$gene
cluster_3_top <- top_markers[31:40, ]$gene
cluster_4_top <- top_markers[41:50, ]$gene
cluster_5_top <- top_markers[51:60, ]$gene
cluster_6_top <- top_markers[61:70, ]$gene
cluster_7_top <- top_markers[71:80, ]$gene
cluster_8_top <- top_markers[81:90, ]$gene
cluster_9_top <- top_markers[91:100, ]$gene



# now normalize and PCA for time 1 and time 6, then look at the data
time_1_p <- process_data(time_1)
# # PCA failed, redo PCA and following steps
# time_1_p <- RunPCA(time_1_p,
#                      seed.use = 42,
#                      verbose = T)
# time_1_p <- RunUMAP(time_1_p,
#                       dims = 1:30,
#                       seed.use = 42,
#                       verbose = T)
# time_1_p <- FindNeighbors(time_1_p,
#                             dims = 1:30,
#                             verbose = T)
# time_1_p <- FindClusters(time_1_p,
#                            verbose = T, 
#                            random.seed = 42,
#                            resolution = 0.14)

DimPlot(time_1_p, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of merged 1hr")

time_6_p <- process_data(time_6)


time_6_p <- FindClusters(time_6_p,
                         verbose = T, 
                         random.seed = 42,
                         resolution = 0.16)
DimPlot(time_6_p, 
        label = T,
        repel = T) +
  NoLegend() +
  ggtitle("clustering of merged 6hr")

# check the expression of top markers in 1hr and 6hr
FeaturePlot(time_1_p, 
            features = cluster_0_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_1_p, 
            features = cluster_1_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_1_p, 
            features = cluster_2_top, 
            pt.size = 0.2,
            ncol = 5,  
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_1_p, 
            features = cluster_3_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_1_p, 
            features = cluster_4_top, 
            pt.size = 0.2,
            ncol = 5,  
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_1_p, 
            features = cluster_5_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_1_p, 
            features = cluster_6_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_1_p, 
            features = cluster_7_top, 
            pt.size = 0.2,
            ncol = 5,  
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_1_p, 
            features = cluster_8_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_1_p, 
            features = cluster_9_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

# 6hr feature plots
FeaturePlot(time_6_p, 
            features = cluster_0_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_6_p, 
            features = cluster_1_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_6_p, 
            features = cluster_2_top, 
            pt.size = 0.2,
            ncol = 5,  
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_6_p, 
            features = cluster_3_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_6_p, 
            features = cluster_4_top, 
            pt.size = 0.2,
            ncol = 5,  
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_6_p, 
            features = cluster_5_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_6_p, 
            features = cluster_6_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_6_p, 
            features = cluster_7_top, 
            pt.size = 0.2,
            ncol = 5,  
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_6_p, 
            features = cluster_8_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_6_p, 
            features = cluster_9_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))


# 0hr feature plots
FeaturePlot(time_0_p, 
            features = cluster_0_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_0_p, 
            features = cluster_1_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_0_p, 
            features = cluster_2_top, 
            pt.size = 0.2,
            ncol = 5,  
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_0_p, 
            features = cluster_3_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_0_p, 
            features = cluster_4_top, 
            pt.size = 0.2,
            ncol = 5,  
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_0_p, 
            features = cluster_5_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_0_p, 
            features = cluster_6_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_0_p, 
            features = cluster_7_top, 
            pt.size = 0.2,
            ncol = 5,  
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_0_p, 
            features = cluster_8_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(time_0_p, 
            features = cluster_9_top, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))
