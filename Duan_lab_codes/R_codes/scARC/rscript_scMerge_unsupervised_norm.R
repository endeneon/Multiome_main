# Siwei 04 March 2022
# Use scMerge+scater to normalise 6 libraries of RNA-seq
# let it run from Rscript

# init
library(SingleCellExperiment)
library(scater)
library(scMerge)
library(BiocParallel)
library(BiocSingular)

library(Seurat)
library(Signac)

library(stringr)
### load data

setwd("/nvmefs/scARC_Duan_018/R_pseudo_bulk_RNA_seq")
load("sce_merged_4_cell_type_only.RData", verbose = T)
data("segList", package = "scMerge")
# load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/human_only_obj_with_cell_type_labeling.RData",
#      verbose = T)
# ###### !!! continue from here #####
# sce_batch_name <- colnames(sce_merged)
# sce_batch_name <- unlist(str_split(sce_batch_name,
#                                    pattern = '_',
#                                    simplify = T)[ ,1])
# # sce_merged$batch <- sce_batch_name
# scater::plotPCA(sce_merged,
#                 colour_by = "batch")
# scater::plotUMAP(sce_merged,
#                  colour_by = "batch")

# assign cell type by barcodes
# unique(human_only$seurat_clusters)
# unique(human_only$lib.ident)
# unique(human_only$group.ident)
# unique(human_only$time.ident)
# unique(human_only$cell.line.ident)

# sce_merged_backup <- sce_merged
# sce_merged <- sce_merged_backup
# reduce matrix size
# hist(colSums(seurat_merged), breaks = 1000, xlim = c(0, 10000))
# sce_merged <- sce_merged[i = 1:nrow(sce_merged),
#                          j = (colSums(sce_merged@assays@data@listData$counts) > 1000) &
#                            (colSums(sce_merged@assays@data@listData$counts)) < 50000]
# 
# # hist(rowSums(sce_merged@assays@data@listData$counts),
# #      breaks = 10000,
# #      xlim = c(0, 5000))
# sce_merged <- 
#   sce_merged[i = rowSums(sce_merged@assays@data@listData$counts) > 1000,
#              j = 1:ncol(sce_merged)]
# 
# 
# sce_barcodes <- colnames(sce_merged)
# sce_barcodes <- unlist(str_split(sce_barcodes,
#                                  pattern = '_',
#                                  simplify = T)[ ,3])
# 
# sce_merged$seurat_clusters <- "unknown"
# # sce_merged$seurat_clusters[str_sub(colnames(sce_merged), 
# #                                    start = 4L) %in%
# #                              unlist(str_split(colnames(human_only)[(human_only$seurat_clusters == "1") &
# #                                                                      (human_only$time.ident == "0hr")],
# #                                               pattern = "_",
# #                                               simplify = T)[, 1])] <- "type_standard_Glut"
# sce_merged$seurat_clusters[str_sub(colnames(sce_merged), 
#                                    start = 4L) %in%
#                              unlist(str_split(colnames(human_only)[(human_only$seurat_clusters == "2") &
#                                                                      (human_only$time.ident == "0hr")],
#                                               pattern = "_",
#                                               simplify = T)[, 1])] <- "type_standard_GABA"
# 
# ## remove unused objects
# rm(human_only)
# # unique(sce_merged$seurat_clusters)
# # convert counts to logcounts (will convert back later)
# # sce_merged@assays@data@listData$logcounts <-
# #   log
# # rm(human_only)
# # sce_merged <- SingleCellExperiment::counts(sce_merged)
# # run a fast and cheap calc with BiocSingular Irlba first
# try({
#   scMerge_semi_supervised_merged_fast <-
#     scMerge(sce_combine = sce_merged,
#             ctl = segList$human$human_scSEG,
#             kmeansK = c(2, 2, 2, 2, 2, 2),
#             assay_name = "scMerge_semi_supervised_Irlba",
#             cell_type = sce_merged$seurat_clusters,
#             cell_type_inc = which(sce_merged$seurat_clusters == "type_standard_GABA"),
#             BSPARAM = IrlbaParam(),
#             svd_k = 20,
#             BPPARAM = MulticoreParam(workers = 8, 
#                                      progressbar = T))
#   
#   save.image()
# })

scMerge_unsupervised_merged <-
  scMerge(sce_combine = sce_merged_4_cell_type,
          # ctl = segList$human$human_scSEG,
          kmeansK = c(4, 4, 4, 4, 4, 4),
          assay_name = "scMerge_semi_supervised",
          cell_type = sce_merged_4_cell_type$cell.type,
          cell_type_inc = NULL,
          cell_type_match = TRUE,
          # cell_type_inc = which(sce_merged$seurat_clusters == "type_standard_GABA"),
          BPPARAM = MulticoreParam(workers = 16, 
                                   progressbar = T))

save(scMerge_unsupervised_merged,
     file = "sce_merged_4_cell_type_normalised_noCtrl.RData",
     compress = T,
     compression_level = 9)
