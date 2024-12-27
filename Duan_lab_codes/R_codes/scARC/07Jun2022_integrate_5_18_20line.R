# Chuxuan Li 06/07/2022
# Combine 5-line, 18-line, and 20-line RNAseq data together with Seurat

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
load("../Analysis_v2_combine_5line_18line/5line_mapped_to_demuxed_barcodes_human_mapped_only_list.RData")
list_5 <- finalobj_lst
rm(finalobj_lst)
load("../Analysis_v1_normalize_by_lib_then_integrate/removed_rat_genes_list.RData")
list_18 <- finalobj_lst
rm(finalobj_lst)
load("../Analysis_v3_combine_18line_20line/20line_demuxed_removed_nonhuman_transformed_lst.RData")
list_20 <- transformed_lst
rm(transformed_lst)
gc()

DefaultAssay(list_5[[1]])
DefaultAssay(list_18[[1]])
DefaultAssay(list_20[[1]])
# integration using reciprocal PCA ####
# find anchors
features <- SelectIntegrationFeatures(object.list = c(list_5, list_18, list_20),
                                      nfeatures = 5000,
                                      fvf.nfeatures = 5000)
combined_lst <- PrepSCTIntegration(c(list_5, list_18, list_20), anchor.features = features)
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
#Found 1436 anchors
# integrate
integrated <- IntegrateData(anchorset = anchors,
                            verbose = T)

save("integrated", file = "5_18_20line_combined_after_integration_nfeat_5000.RData")
