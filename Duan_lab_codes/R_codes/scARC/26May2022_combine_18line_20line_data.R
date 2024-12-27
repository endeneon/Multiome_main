#Chuxuan li 05/26/2022
#Combine RNAseq data of 18line and 20 line, normalize with Seurat independently
#for each library, then integrate, and use Harmony to minimize batch effect

# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)

library(stringr)
library(readr)
library(dplyr)
library(future)

library(ggplot2)
library(RColorBrewer)
library(viridis)
library(graphics)
library(patchwork)
library(cowplot)
library(gplots)

# set threads and parallelization
plan("multisession", workers = 1)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 107374182400)

# load data ####
rm(list = ls())
load("../Analysis_v1_normalize_by_lib_then_integrate/removed_rat_genes_list.RData")
load("./20line_demuxed_removed_nonhuman_transformed_lst.RData")
list_18 <- finalobj_lst
rm(finalobj_lst)
list_20 <- transformed_lst
rm(transformed_lst)
gc()

DefaultAssay(list_18[[1]])
DefaultAssay(list_20[[1]])
# integration using reciprocal PCA ####
# find anchors
features <- SelectIntegrationFeatures(object.list = c(list_18, list_20),
                                      nfeatures = 5000,
                                      fvf.nfeatures = 5000)
combined_lst <- PrepSCTIntegration(c(list_18, list_20), anchor.features = features)
combined_lst <- lapply(X = combined_lst, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(object.list = combined_lst,
                                  anchor.features = features,
                                  reference = c(1, 2, 3),
                                  reduction = "rpca",
                                  normalization.method = "SCT",
                                  scale = F, 
                                  dims = 1:50)
#Found 4314 anchors
# integrate
integrated <- IntegrateData(anchorset = anchors,
                            verbose = T)

save("integrated", file = "18line_20line_combined_after_integration_nfeat_5000.RData")
