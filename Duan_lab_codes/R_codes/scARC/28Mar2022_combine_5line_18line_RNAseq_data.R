# Chuxuan Li 02/28/2022
# Merge 5 line and 18 line data and redo analysis

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
load("../Analysis_v1_normalize_by_lib_then_integrate/removed_rat_genes_list.RData")
load("./5line_mapped_to_demuxed_barcodes_human_mapped_only_list.RData")

# integration using reciprocal PCA ####
# find anchors
features <- SelectIntegrationFeatures(object.list = c(finalobj_lst, cleanobj_lst),
                                      nfeatures = 5000,
                                      fvf.nfeatures = 5000)
combined_lst <- PrepSCTIntegration(c(finalobj_lst, cleanobj_lst), anchor.features = features)
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
#Found 9036 anchors
# integrate
integrated <- IntegrateData(anchorset = anchors,
                            verbose = T)

save("integrated", file = "5line_18line_combined_after_integration_nfeat_5000.RData")
