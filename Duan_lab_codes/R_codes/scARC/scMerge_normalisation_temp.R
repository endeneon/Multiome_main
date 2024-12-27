# Siwei 01 March 2022
# Use scMerge+scater to normalise 6 libraries of RNA-seq

# init
library(SingleCellExperiment)
library(scater)
library(batchelor)
library(scran)

library(scMerge)
library(BiocParallel)
library(BiocSingular)

library(Seurat)
library(Signac)

library(stringr)

library(future)
plan("multisession", workers = 4)
options(future.globals.maxSize = 200 * 1024 ^ 3)

# plan()
### load data
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/raw_data_seurat_list.RData")
data("segList", package = "scMerge")
###

# assign ident to each library
i <- 1L
library_list <- c("g_2_0", "g_2_1", "g_2_6",
                  "g_8_0", "g_8_1", "g_8_6")

### need to merge the objects at the Seurat level,
### then convert to SingleCellExperiment
for (i in 1:length(obj_lst)) {
  obj_lst[[i]]$orig.ident <- library_list[i]
}

sce_list <- vector(mode = "list",
                   length = length(obj_lst))

# rename the cells so that when merging, weird subscripts are not added
for (i in 1:length(obj_lst)) {
  print(i)
  obj_lst[[i]] <- 
    RenameCells(obj_lst[[i]], 
                add.cell.id = paste0(i, "_"))
}

# 
i <- 1L
for (i in 1:length(obj_lst)) {
  print(i)
  if (i == 1) {
    seurat_merged <- obj_lst[[i]]
  } else {
    seurat_merged <-
      merge(seurat_merged,
            obj_lst[[i]])
  }
}

colnames(seurat_merged)

# convert Seurat object to SingleCellExperiments
sce_merged <- as.SingleCellExperiment(seurat_merged)

# for (i in 1:length(obj_lst)) {
#   print(i)
#   sce_list[[i]] <- as.SingleCellExperiment(obj_lst[[i]])
# }

# combine the sce object into one
# sce_merged <-
#   sce_cbind(sce_list = sce_list,
#             method = "intersect",
#             exprs = "counts")
colnames(sce_merged)

# change the human cell barcodes to conform changed sce barcodes
# make a df from human cells that conforms the format
# load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/human_only_obj_with_cell_type_labeling.RData",
     # verbose = T)
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/obj_list_after_SCT.RData")
i <- 1L
for (i in 1:6) {
  print(i)
  temp_seurat <-
    subset(human_only,
           subset = lib.ident == i)
  
  temp_seurat_barcode = colnames(temp_seurat)
  temp_seurat_barcode <-
    unlist(str_split(string = temp_seurat_barcode,
                     pattern = "_",
                     simplify = T)[, 1])
  temp_seurat_barcode <-
    paste0(i, "__", temp_seurat_barcode)
  
  if (i == 1) {
    human_cell_barcodes <- temp_seurat_barcode
  } else {
    human_cell_barcodes <- 
      c(human_cell_barcodes, temp_seurat_barcode)
  }
  
}


# human_cells_barcodes <-
#   unlist(str_split(string = human_cells_barcodes,
#                    pattern = "_",
#                    simplify = T)[, 1])

# take only human cells from the sce object
# sce_merged <-
#   sce_merged[, colnames(sce_merged) %in% human_cell_barcodes]
# 
# sce_merged <-
#   runPCA(sce_merged,
#          exprs_values = "counts",
#          BPPARAM = MulticoreParam(workers = 8, 
#                                   progressbar = T))
# sce_merged <-
#   runUMAP(sce_merged,
#          exprs_values = "counts",
#          BPPARAM = MulticoreParam(workers = 8, 
#                                   progressbar = T))
# 



###### !!! continue from here #####
sce_batch_name <- colnames(sce_merged)
sce_batch_name <- unlist(str_split(sce_batch_name,
                                   pattern = '_',
                                   simplify = T)[ ,1])
sce_merged$batch <- sce_batch_name
# scater::plotPCA(sce_merged,
#                 colour_by = "batch")
# scater::plotUMAP(sce_merged,
#                  colour_by = "batch")

# assign cell type by barcodes
unique(human_only$seurat_clusters)
unique(human_only$lib.ident)
unique(human_only$group.ident)
unique(human_only$time.ident)
unique(human_only$cell.line.ident)

# sce_merged_backup <- sce_merged
sce_merged <- sce_merged_backup
# reduce matrix size
hist(colSums(seurat_merged), breaks = 1000, xlim = c(0, 10000))
sce_merged <- sce_merged[i = 1:nrow(sce_merged),
                         j = (colSums(sce_merged@assays@data@listData$counts) > 1000) &
                           (colSums(sce_merged@assays@data@listData$counts)) < 50000]

hist(rowSums(sce_merged@assays@data@listData$counts),
     breaks = 10000,
     xlim = c(0, 5000))
sce_merged <- 
  sce_merged[i = rowSums(sce_merged@assays@data@listData$counts) > 1000,
             j = 1:ncol(sce_merged)]


sce_barcodes <- colnames(sce_merged)
sce_barcodes <- unlist(str_split(sce_barcodes,
                                   pattern = '_',
                                   simplify = T)[ ,3])

sce_merged$seurat_clusters <- "unknown"
# sce_merged$seurat_clusters[str_sub(colnames(sce_merged), 
#                                    start = 4L) %in%
#                              unlist(str_split(colnames(human_only)[(human_only$seurat_clusters == "1") &
#                                                                      (human_only$time.ident == "0hr")],
#                                               pattern = "_",
#                                               simplify = T)[, 1])] <- "type_standard_Glut"
sce_merged$seurat_clusters[str_sub(colnames(sce_merged), 
                                   start = 4L) %in%
                             unlist(str_split(colnames(human_only)[(human_only$seurat_clusters == "2") &
                                                                     (human_only$time.ident == "0hr")],
                                              pattern = "_",
                                              simplify = T)[, 1])] <- "type_standard_GABA"

unique(sce_merged$seurat_clusters)
# convert counts to logcounts (will convert back later)
# sce_merged@assays@data@listData$logcounts <-
#   log
rm(human_only)
# sce_merged <- SingleCellExperiment::counts(sce_merged)
# run a fast and cheap calc with BiocSingular Irlba first
scMerge_semi_supervised_merged_fast <-
  scMerge(sce_combine = sce_merged,
          ctl = segList$human$human_scSEG,
          kmeansK = c(2, 2, 2, 2, 2, 2),
          assay_name = "scMerge_semi_supervised_Irlba",
          cell_type = sce_merged$seurat_clusters,
          cell_type_inc = which(sce_merged$seurat_clusters == "type_standard_GABA"),
          BSPARAM = IrlbaParam(),
          svd_k = 20,
          BPPARAM = MulticoreParam(workers = 8, 
                                   progressbar = T))

save.image()

scMerge_semi_supervised_merged <-
  scMerge(sce_combine = sce_merged,
          ctl = segList$human$human_scSEG,
          kmeansK = c(2, 2, 2, 2, 2, 2),
          assay_name = "scMerge_semi_supervised",
          cell_type = sce_merged$seurat_clusters,
          cell_type_inc = which(sce_merged$seurat_clusters == "type_standard_GABA"),
          BPPARAM = MulticoreParam(workers = 8, 
                                   progressbar = T))

save.image()
####
# data("example_sce", package = "scMerge")
# test <- example_sce

# test sparse array conversion
library(Matrix)
library(DelayedArray)

sparse_input = sce_merged

assay(sparse_input, "counts") = as(counts(sparse_input), "dgeMatrix")
# assay(sparse_input, "logcounts") = as(logcounts(sparse_input), "dgeMatrix")


## check normalised object
sce_merged <- runPCA(sce_merged, 
                     BPPARAM = MulticoreParam(workers = 8, 
                                              progressbar = T))
scater::plotPCA(sce_merged)

scMerge_semi_supervised_merged <-
  scater::runUMAP(scMerge_semi_supervised_merged,
                  BPPARAM = MulticoreParam(workers = 8, 
                                           progressbar = T))
scater::plotUMAP(scMerge_semi_supervised_merged,
                 colour_by = "batch",
                 size = 0.2)
seurat_scMerge_semi_supervised_merged <-
  as.Seurat(scMerge_semi_supervised_merged)

library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize = 200 * 1024 ^ 3)

seurat_scMerge_semi_supervised_merged <-
  Seurat::RunUMAP(seurat_scMerge_semi_supervised_merged)

DimPlot(seurat_scMerge_semi_supervised_merged,
        reduction = "UMAP",
        group.by = "batch",
        pt.size = 0.1)


sce_merged$cell.type <- "unknown"
sce_merged$cell.type[str_sub(colnames(sce_merged), 
                             start = 4L) %in%
                       unlist(str_split(colnames(human_only)[human_only$cell.type %in% "GABA"],
                                        pattern = "_",
                                        simplify = T)[, 1])] <- "GABA"
sce_merged$cell.type[str_sub(colnames(sce_merged), 
                             start = 4L) %in%
                       unlist(str_split(colnames(human_only)[human_only$cell.type %in% "NEFM_pos_glut"],
                                        pattern = "_",
                                        simplify = T)[, 1])] <- "NEFM_pos_glut"
sce_merged$cell.type[str_sub(colnames(sce_merged), 
                             start = 4L) %in%
                       unlist(str_split(colnames(human_only)[human_only$cell.type %in% "NEFM_neg_glut"],
                                        pattern = "_",
                                        simplify = T)[, 1])] <- "NEFM_neg_glut"
sce_merged$cell.type[str_sub(colnames(sce_merged), 
                             start = 4L) %in%
                       unlist(str_split(colnames(human_only)[human_only$cell.type %in% "NPC"],
                                        pattern = "_",
                                        simplify = T)[, 1])] <- "NPC"

cell_type_list <- unique(human_only$cell.type)
sce_merged_4_cell_type <-
  sce_merged[, !(sce_merged$cell.type %in% "unknown")]

save(sce_merged_4_cell_type,
     file = "sce_merged_4_cell_type_only.RData",
     compress = T)

sce_merged_4_cell_type_seurat <- as.Seurat(scMerge_unsupervised_merged)

# sce_merged_4_cell_type_seurat <-
#   FindVariableFeatures(sce_merged_4_cell_type_seurat,
#                        selection.method = 'vst',
#                        nfeatures = 2000)

human_scSEG_filtered <-
  segList$human$human_scSEG[segList$human$human_scSEG %in% 
                              rownames(sce_merged_4_cell_type_seurat)]
sce_merged_4_cell_type_seurat[["percent.human_scSEG"]] <-
  PercentageFeatureSet(object = sce_merged_4_cell_type_seurat,
                       features = human_scSEG_filtered)


sce_merged_4_cell_type_seurat <-
  SCTransform(sce_merged_4_cell_type_seurat,
            vars.to.regress = "percent.human_scSEG",
            verbose = T)

sce_merged_4_cell_type_seurat <-
  Seurat::RunPCA(sce_merged_4_cell_type_seurat,
                 npcs = 30,
                 verbose = T)

sce_merged_4_cell_type_seurat@reductions <- NULL

sce_merged_4_cell_type_seurat <-
  Seurat::RunUMAP(sce_merged_4_cell_type_seurat,
                  dims = 1:30,
                  verbose = T)

DimPlot(sce_merged_4_cell_type_seurat,
        reduction = 'umap',
        group.by = "batch")


scMerge_unsupervised_merged <-
  scater::runUMAP(scMerge_unsupervised_merged,
                  name = 'UMAP')
scater::plotUMAP(scMerge_unsupervised_merged,
                 # name = 'umap',
                 colour_by = "batch")
