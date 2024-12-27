# Siwei 09 March 2022
# Use scAlign to normalise 6 libraries of RNA-seq

## init
# init
library(SingleCellExperiment)
library(scater)
library(batchelor)
library(scran)

library(BiocParallel)
library(BiocSingular)

library(Seurat)
library(Signac)

library(scAlign)
library(scMerge)

library(stringr)
library(sctransform)
library(glmGamPoi)

library(future)
library(parallel)
library(tensorflow)

plan("multisession", workers = 8)
options(future.globals.maxSize = 200 * 1024 ^ 3)

setwd("/home/zhangs3/NVME/scARC_Duan_018/R_pseudo_bulk_RNA_seq")

###
data("segList", package = "scMerge")
load(file = "scAlign_merge_6_objects.RData")

## create the combined scAlign object, add labels later
scAlign_6_objects_multi <- 
  scAlignCreateObject(sce.objects = list("g_2_0" = sce_lst[[1]],
                                         "g_2_1" = sce_lst[[2]],
                                         "g_2_6" = sce_lst[[3]],
                                         "g_8_0" = sce_lst[[4]],
                                         "g_8_1" = sce_lst[[5]],
                                         "g_8_6" = sce_lst[[6]]),
                      labels = list(sce_lst[[1]]@colData$cell.type,
                                    sce_lst[[2]]@colData$cell.type,
                                    sce_lst[[3]]@colData$cell.type,
                                    sce_lst[[4]]@colData$cell.type,
                                    sce_lst[[5]]@colData$cell.type,
                                    sce_lst[[6]]@colData$cell.type),
                      data.use = "scale.data",
                      pca.reduce = T,
                      pcs.compute = 50,
                      genes.use = common_genes,
                      cca.reduce = T,
                      ccs.compute = 15,
                      project.name = "scAlign_6_objects_multi")


save.image(file = "scAlign_merge_6_objects.RData")
# 
# load(file = "scAlign_merge_6_objects.RData")
# # run all samples in one pass
scAlign_6_objects_multi_aligned <-
  scAlignMulti(scAlign_6_objects_multi,
               options = scAlignOptions(steps = 15000,
                                        log.every = 5000,
                                        batch.size = 300,
                                        perplexity = 30,
                                        norm = TRUE,
                                        batch.norm.layer = FALSE,
                                        architecture = "large",  ## 3 layer neural network
                                        num.dim = 64),            ## Number of latent dimensions)
               encoder.data = "scale.data",
               # decoder.data = "scale.data",
               reference.data = "g_2_0",
               supervised = 'none',
               run.encoder = TRUE,
               run.decoder = F,
               device = "GPU")

save.image(file = "scAlign_merge_6_objects.RData")

# try supervised align
scAlign_6_objects_multi_aligned_supervised <-
  scAlignMulti(scAlign_6_objects_multi,
               options = scAlignOptions(steps = 15000,
                                        log.every = 5000,
                                        batch.size = 300,
                                        perplexity = 30,
                                        norm = TRUE,
                                        batch.norm.layer = FALSE,
                                        architecture = "large",  ## 3 layer neural network
                                        num.dim = 64),            ## Number of latent dimensions)
               encoder.data = "scale.data",
               # decoder.data = "scale.data",
               reference.data = "g_2_0",
               supervised = 'both',
               run.encoder = TRUE,
               run.decoder = F,
               device = "GPU")
