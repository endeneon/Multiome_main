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
# options(mc.core)

# human_only
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/human_only_traditional_normed_obj_with_cell_type_labeling.RData")
data("segList", package = "scMerge")

# load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/mapped_to_demuxed_barcodes_pure_human_obj.RData")
###
# assign ident to each library
i <- 1L


# obj_lst <- human_only
# make a sce_list object
library_list <- c("g_2_0", "g_2_1", "g_2_6",
                  "g_8_0", "g_8_1", "g_8_6")
Idents(human_only) <- "lib.ident"
obj_lst <- 
  SplitObject(human_only,
              split.by = "lib.ident")

for (i in 1:6) {
  DefaultAssay(obj_lst[[i]]) <- "RNA"
}

plan("multisession", workers = 8)

i <- 1L
for (i in 1:6) {
  obj_lst[[i]] <- 
    PercentageFeatureSet(obj_lst[[i]],
                         pattern = "^MT-",
                         col.name = "percent.mt")
  obj_lst[[i]] <- 
    SCTransform(obj_lst[[i]],
                vars.to.regress = "percent.mt",
                variable.features.n = 16000, 
                do.scale = T, 
                return.only.var.genes = T,
                method = "glmGamPoi",
                seed.use = 42,
                verbose = T)
}

head(obj_lst[[1]])
head(obj_lst[[1]]@assays$SCT@counts)

DefaultAssay(obj_lst[[1]])
common_genes <-
  Reduce(intersect,
         list(VariableFeatures(obj_lst[[1]]),
              VariableFeatures(obj_lst[[2]]),
              VariableFeatures(obj_lst[[3]]),
              VariableFeatures(obj_lst[[4]]),
              VariableFeatures(obj_lst[[5]]),
              VariableFeatures(obj_lst[[6]])))

sce_lst <- 
  lapply(obj_lst,
           function(seurat.obj) {
             SingleCellExperiment(assays = list(count = seurat.obj@assays$SCT@counts[common_genes, ],
                                                logcounts = seurat.obj@assays$SCT@counts[common_genes, ],
                                                scale.data = seurat.obj@assays$SCT@scale.data[common_genes, ]),
                                  colData = seurat.obj@meta.data)
           })

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
               reference.data = 
               supervised = 'none',
               run.encoder = TRUE,
               run.decoder = TRUE,
               device = "GPU")

save.image(file = "scAlign_merge_6_objects.RData")










# rename the cells so that when merging, weird subscripts are not added
# i <- 1L
# library_list <- c("g_2_0", "g_2_1", "g_2_6",
#                   "g_8_0", "g_8_1", "g_8_6")
# 
# for (i in 1:length(obj_lst)) {
#   print(i)
#   obj_lst[[i]] <- 
#     RenameCells(obj_lst[[i]], 
#                 add.cell.id = paste0(i, "_"))
#   obj_lst[[i]]$library.list <-
#     library_list[i]
# }
# 
# # 
# i <- 1L
# for (i in 1:length(obj_lst)) {
#   print(i)
#   if (i == 1) {
#     seurat_merged <- obj_lst[[i]]
#   } else {
#     seurat_merged <-
#       merge(seurat_merged,
#             obj_lst[[i]])
#   }
# }
# 
# colnames(seurat_merged)
# 
# # convert Seurat object to SingleCellExperiments
# sce_merged <- as.SingleCellExperiment(seurat_merged)
# #########
# 
# rm(obj_lst)
# rm(sce_merged)
# rm(seurat_merged)
# 
# DefaultAssay(human_only) <- "RNA"
# human_only@assays$integrated <- NULL
# 
# sce_merged <- as.SingleCellExperiment(human_only)
# sce_merged@assays@data@listData$counts
# 
# 
# 
# scAlign_merged <-
#   scAlignCreateObject(sce.objects = seurat_merged@assays$RNA@data,
#                       labels = list(seurat_merged$lib.ident),
#                       pca.reduce = T,
#                       pcs.compute = 50,
#                       cca.reduce = T,
#                       ccs.compute = 15,
#                       project.name = "5_lines")
# 
# scAlign_merged <-
#   scAlign(sce.object = sce_merged,
#           options = scAlignOptions(steps = 15000,
#                                    norm = T,
#                                    early.stop = T),
#           encoder.data = "counts", 
#           supervised = "none",
#           run.encoder = T,
#           run.decoder = F,
#           device = "GPU")
# 
# ################
# colnames(pure_human)
# colnames(obj_lst[[1]])
# 
# unique(str_sub(string = colnames(pure_human),
#                end = -2L))
# 
# sce_list <- vector(mode = "list", length = 6L)
# i <- 1L
# for (i in 1:length(obj_lst)) {
#   print(i)
#   temp_data_matrix <- 
#     as.matrix(obj_lst[[i]]@assays$RNA@data)
#   rownames(temp_data_matrix) <-
#              rownames(obj_lst[[i]])
#   colnames(temp_data_matrix) <-
#              colnames(obj_lst[[i]])
#   sce_list[[i]] <- 
#     SingleCellExperiment(assay = list(counts = temp_data_matrix))
# }
# 
# 
# 
# scAlign_merged <-
#   scAlignCreateObject(sce.objects = list(sce_list[[1]],
#                                          sce_list[[2]],
#                                          sce_list[[3]],
#                                          sce_list[[4]],
#                                          sce_list[[5]],
#                                          sce_list[[6]]),
#                       labels = list("g_2_0", "g_2_1", "g_2_6",
#                                     "g_8_0", "g_8_1", "g_8_6"),
#                       data.use = "counts",
#                       pca.reduce = T,
#                       pcs.compute = 50,
#                       cca.reduce = T,
#                       ccs.compute = 15,
#                       project.name = "5_lines")
# 
# scAlign_merged <-
#   scAlign(sce.object = sce_merged,
#           options = scAlignOptions(steps = 15000,
#                                    norm = T,
#                                    early.stop = T),
#           encoder.data = "counts", 
#           supervised = "none",
#           run.encoder = T,
#           run.decoder = F,
#           device = "GPU")
# 
# ########
# 
# seurat_0hr <- scM
