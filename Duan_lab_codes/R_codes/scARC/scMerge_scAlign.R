# Siwei 10 Mar 2022
# scMerge + scAlign

# init

library(SingleCellExperiment)

library(BiocParallel)
library(BiocSingular)

library(Seurat)
library(Signac)

library(scMerge)
library(scAlign)

library(stringr)

library(future)
plan("multisession", workers = 4)
options(future.globals.maxSize = 200 * 1024 ^ 3)

# load data
# load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/raw_data_seurat_list.RData")
# integrated
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/normalize_by_6libs_integrated.RData")

# human_only
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/human_only_traditional_normed_obj_with_cell_type_labeling.RData")
data("segList", package = "scMerge")


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
  obj_lst[[i]]@assays$integrated <- NULL
  
}

seurat_3_time_list <- 
  vector(mode = "list", length = 3L)

seurat_3_time_list[[1]] <-
  merge(obj_lst[[1]],
        obj_lst[[4]])

seurat_3_time_list[[2]] <-
  merge(obj_lst[[2]],
        obj_lst[[5]])

seurat_3_time_list[[3]] <-
  merge(obj_lst[[3]],
        obj_lst[[6]])

# sce_list <- vector(mode = "list", length = 6L)
# library_list <- c("g_0hr", "g_1hr", "g_6hr")

for (i in 1:3) {
  DefaultAssay(seurat_3_time_list[[i]]) <- "RNA"
}

# i <- 1L
# for (i in 1:length(obj_lst)) {
#   print(i)
#   temp_data_matrix <- 
#     as.matrix(obj_lst[[i]]@assays$RNA@data)
#   rownames(temp_data_matrix) <-
#     rownames(obj_lst[[i]])
#   colnames(temp_data_matrix) <-
#     colnames(obj_lst[[i]])
#   sce_list[[i]] <- 
#     SingleCellExperiment(assay = list(counts = temp_data_matrix))
#   sce_list[[i]]$batch <-
#     library_list[i]
#   sce_list[[i]]$cellTypes <- "unknown"
#   sce_list[[i]]$cellTypes[obj_lst[[i]]$cell.type %in% "NEFM_pos_glut"] <-
#     "NEFM_pos_glut"
#   sce_list[[i]]$cellTypes[obj_lst[[i]]$cell.type %in% "GABA"] <-
#     "GABA"
#   sce_list[[i]]$cellTypes[obj_lst[[i]]$cell.type %in% "NEFM_neg_glut"] <-
#     "NEFM_neg_glut"
#   sce_list[[i]]$cellTypes[obj_lst[[i]]$cell.type %in% "NPC"] <-
#     "NPC"
# }


sce_3_time_list <- vector(mode = "list",
                          length = 3L)

i <- 1L
for (i in 1:length(seurat_3_time_list)) {
  print(i)
  temp_data_matrix <- 
    as.matrix(seurat_3_time_list[[i]]@assays$RNA@data)
  rownames(temp_data_matrix) <-
    rownames(seurat_3_time_list[[i]])
  colnames(temp_data_matrix) <-
    colnames(seurat_3_time_list[[i]])
  sce_3_time_list[[i]] <- 
    SingleCellExperiment(assay = list(counts = temp_data_matrix),
                         colData = seurat_3_time_list[[i]]@meta.data)
  sce_3_time_list[[i]]$batch <-
    seurat_3_time_list[[i]]$lib.ident
  sce_3_time_list[[i]]$cellTypes <- "unknown"
  sce_3_time_list[[i]]$cellTypes[seurat_3_time_list[[i]]$cell.type %in% "NEFM_pos_glut"] <-
    "NEFM_pos_glut"
  sce_3_time_list[[i]]$cellTypes[seurat_3_time_list[[i]]$cell.type %in% "GABA"] <-
    "GABA"
  sce_3_time_list[[i]]$cellTypes[seurat_3_time_list[[i]]$cell.type %in% "NEFM_neg_glut"] <-
    "NEFM_neg_glut"
  sce_3_time_list[[i]]$cellTypes[seurat_3_time_list[[i]]$cell.type %in% "NPC"] <-
    "NPC"
}


# merge each time point by scMerge
# scMerge_semi_supervised_merged_0hr <-
#   scMerge(sce_combine = sce_3_time_list[[1]],
#           ctl = segList$human$human_scSEG,
#           kmeansK = c(4, 4),
#           assay_name = "scMerge_semi_supervised",
#           exprs = "counts",
#           cell_type = sce_3_time_list$cellTypes,
#           cell_type_inc = which(sce_3_time_list$cellTypes != "unknown"),
#           cell_type_match = T,
#           BPPARAM = MulticoreParam(workers = 8, 
#                                    progressbar = T))
# 
# scMerge_semi_supervised_merged_0hr <-
#   runPCA(scMerge_semi_supervised_merged_0hr,
#          exprs_values = "scMerge_semi_supervised",
#          BPPARAM = MulticoreParam(workers = 8, 
#                                   progressbar = T))
# 
# scater::plotPCA(scMerge_semi_supervised_merged_0hr,
#                 colour_by = "batch")

scMerge_semi_supervised_merged <- vector(mode = "list", length = 3L)
# scMerge_semi_supervised_merged[[1]] <- scMerge_semi_supervised_merged_0hr

group_list <- c("g_0hr", "g_1hr", "g_6hr")

# use scMerge to normalise each samples of each timepoint
i <- 1L
for (i in (1:3)) {
  print(i)
  scMerge_semi_supervised_merged[[i]] <-
    scMerge(sce_combine = sce_3_time_list[[i]],
            ctl = segList$human$human_scSEG,
            kmeansK = c(4, 4),
            assay_name = "scMerge_semi_supervised",
            exprs = "counts",
            cell_type = sce_3_time_list$cellTypes,
            cell_type_inc = which(sce_3_time_list$cellTypes != "unknown"),
            cell_type_match = T,
            BPPARAM = SnowParam(workers = 16, 
                                progressbar = T))
  
}

save.image(file = "scMerge_scAlign.RData")
# Use scAlign to normalise samples between time points and DE analysis

