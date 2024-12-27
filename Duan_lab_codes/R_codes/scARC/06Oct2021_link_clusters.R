# Chuxuan Li
# Clustering the datasets from individual time points, use cell type markers
# to link the clusters

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
                             resolution = 0.5)
  
  return(seurat_obj)
}


time_0_p <- process_data(time_0)
time_1_p <- process_data(time_1)
time_6_p <- process_data(time_6)

markers_for_vln <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
                     "EBF1", "SERTAD4", # striatal
                     "FOXG1", "DLX2",  # forebrain 
                     "POU3F2", # cortical, excitatory
                     "NPY", "DLX5", "SST", "CALB1", "CALB2", "PVALB", # inhibitory
                     "SATB2", "CBLN2", "CUX2", "NEFM", # excitatory
                     "TBR1", # pallial glutamatergic
                     "VIM", "SOX2", "NES", #NPC
                     "PAX6", "PPP1R1B")


StackedVlnPlot(obj = time_0_p, features = markers_for_vln) +
  coord_flip()

StackedVlnPlot(obj = time_1_p, features = markers_for_vln) +
  coord_flip()

StackedVlnPlot(obj = time_6_p, features = markers_for_vln) +
  coord_flip()



transparent_0hr_dimplot <- DimPlot(time_0_p, 
                     label = T) +
  ggtitle("clustering at 0hr")

transparent_0hr_dimplot[[1]]$layers[[1]]$aes_params$alpha = .1
transparent_0hr_dimplot

transparent_1hr_dimplot <- DimPlot(time_1_p, 
                                   label = T) +
  ggtitle("clustering at 1hr")

transparent_1hr_dimplot[[1]]$layers[[1]]$aes_params$alpha = .1
transparent_1hr_dimplot

transparent_6hr_dimplot <- DimPlot(time_6_p, 
                                   label = T) +
  ggtitle("clustering at 6hr")

transparent_6hr_dimplot[[1]]$layers[[1]]$aes_params$alpha = .1
transparent_6hr_dimplot

