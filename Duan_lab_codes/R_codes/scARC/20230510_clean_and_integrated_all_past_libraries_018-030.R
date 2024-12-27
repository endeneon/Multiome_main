# Chuxuan Li 05/10/2023
# Get all libraries from 018 to 030, normalize them separately, then integrate

# init ####
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(stringr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)
library(readxl)
library(readr)

library(future)

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 161061273600)

# load h5 file paths ####
h5list <- sort(list.files(path = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only", 
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T, include.dirs = T, full.names = T))
h5list
pathlist <- c(h5list)
h5list <- sort(list.files(path = "/nvmefs/scARC_Duan_022/cellranger_output", 
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T, include.dirs = T, full.names = T))
h5list
pathlist <- c(pathlist, h5list)
h5list <- sort(list.files(path = "/nvmefs/scARC_Duan_024_raw_fastq/GRCh38_output", 
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T, include.dirs = T, full.names = T))
h5list
pathlist <- c(pathlist, h5list)
h5list <- sort(list.files(path = "/nvmefs/scARC_Duan_025_GRCh38", 
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T, include.dirs = T, full.names = T))
h5list
pathlist <- c(pathlist, h5list[2:16])
h5list <- sort(list.files(path = "/data/FASTQ/Duan_Project_025_17_46", 
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T, include.dirs = T, full.names = T))
h5list
pathlist <- c(pathlist, h5list[1:6])
h5list <- sort(list.files(path = "/data/FASTQ/Duan_Project_029/GRCh38_only", 
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T, include.dirs = T, full.names = T))
h5list
pathlist <- c(pathlist, h5list[1:15])
h5list <- sort(list.files(path = "/data/FASTQ/Duan_Project_029_40201", 
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T, include.dirs = T, full.names = T))
h5list
pathlist <- c(pathlist, h5list[1:3])
h5list <- sort(list.files(path = "/data/FASTQ/Duan_Project_030", 
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T, include.dirs = T, full.names = T))
h5list
pathlist <- c(pathlist, h5list[1:12])


libnames <- pathlist
libnames[1:6] <- str_replace_all(libnames[1:6], "_", "-")
libnames <- str_extract(libnames, "[0-9]+-[0|1|6]")
names(pathlist) <- libnames
pathlist

# read demuxed barcodes ####

# 018
bc_path <- sort(list.files(path = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/common_barcodes",
                           full.names = T, pattern = ".best", recursive = T))
barcode_list <- vector(mode = "list", length = 6)
for (i in 1:length(bc_path)){
  group <- str_extract(bc_path[i], "g_[2|8]_[0|1|6]")
  if (i == 1) {
    prev.group <- "1st_group"
    list_ind <- 0
    groups <- c(group)
  } else {
    groups <- c(groups, group)
  }
  cat("group: ", group, " previous group: ", prev.group)
  if (prev.group != group) {
    list_ind <- list_ind + 1
    df <- data.frame(barcode = read.delim(bc_path[i], header = F, row.names = NULL),
                     line = str_extract(bc_path[i], "CD_[0-9][0-9]$"))
    colnames(df) <- c("barcode", "line")
  } else {
    df_to_add <- data.frame(barcode = read.delim(bc_path[i], header = F, row.names = NULL),
                            line = str_extract(bc_path[i], "CD_[0-9][0-9]$"))
    colnames(df_to_add) <- c("barcode", "line")
    df <- rbind(df, df_to_add)
    barcode_list[[list_ind]] <- df
  }
  print(list_ind)
  prev.group <- group
}
names(barcode_list) <- str_remove(unique(groups), "^g_")

# 022
bc_path <- sort(list.files(path = "/nvmefs/scARC_Duan_022/demux_barcode_output/output_barcodes",
                           full.names = T, pattern = ".tsv", recursive = T))
barcode_list_add <- vector(mode = "list", length = length(bc_path))
for (i in 1:length(bc_path)){
  if (i == 1) {
    groups <- c(str_extract(bc_path[i], "[0-9][0-9]_[0|1|6]"))
  } else {
    groups <- c(groups, str_extract(bc_path[i], "[0-9][0-9]_[0|1|6]"))
  }
  barcode_list_add[[i]] <- read.delim(bc_path[i], header = F, row.names = NULL, skip = 1)
  colnames(barcode_list_add[[i]]) <- c("barcode", "line")
}
names(barcode_list_add) <- groups
groups_022 <- groups # saving this for later removing repeated cell lines
barcode_list <- c(barcode_list, barcode_list_add)

# 024
bc_path <- list.files(path = "/nvmefs/scARC_Duan_024_raw_fastq/common_barcodes",
                      full.names = T, pattern = ".best.tsv", recursive = T)
barcode_list_add <- vector(mode = "list", length = length(bc_path))
for (i in 1:length(bc_path)){
  if (i == 1) {
    groups <- c(str_extract(bc_path[i], "[0-9][0-9]_[a-z]+\\.[0-2]"))
  } else {
    groups <- c(groups, str_extract(bc_path[i], "[0-9][0-9]_[a-z]+\\.[0-2]"))
  }
  barcode_list_add[[i]] <- read.delim(bc_path[i], header = F, row.names = NULL)
  colnames(barcode_list_add[[i]]) <- c("barcode", "line")
}
groups <- str_replace_all(groups, "\\.2", "\\.6")
groups <- str_replace_all(groups, "[a-z]+\\.", "")
names(barcode_list_add) <- groups
barcode_list <- c(barcode_list, barcode_list_add)

# 025
bc_path <- sort(list.files(path = "/nvmefs/scARC_Duan_025_GRCh38/barcodes_demuxed_by_library",
                           full.names = T, pattern = ".best.tsv", recursive = T))
barcode_list_add <- vector(mode = "list", length = length(bc_path))
for (i in 1:length(bc_path)){
  if (i == 1) {
    groups <- c(str_extract(bc_path[i], "[0-9][0-9]_[a-z]+\\.[0-2]"))
  } else {
    groups <- c(groups, str_extract(bc_path[i], "[0-9][0-9]_[a-z]+\\.[0-2]"))
  }
  barcode_list_add[[i]] <- read.delim(bc_path[i], header = F, row.names = NULL)
  colnames(barcode_list_add[[i]]) <- c("barcode", "line")
}
groups <- str_replace(groups, "\\.2", "\\.6")
groups <- str_remove(groups, "atac\\.")
names(barcode_list_add) <- groups
barcode_list <- c(barcode_list, barcode_list_add)

# 025 17&46
bc_path <- sort(list.files(path = "/data/FASTQ/Duan_Project_025_17_46/barcode_demux_output",
                           full.names = T, pattern = ".best.tsv", recursive = T))
barcode_list_add <- vector(mode = "list", length = length(bc_path))
for (i in 1:length(bc_path)){
  if (i == 1) {
    groups <- c(str_remove(str_remove(str_extract(bc_path[i], 
                                                  "pool_[0-9][0-9]_.+\\.[0-2]"), 
                                      "pool_"),
                           "_CD.+atac"))
  } else {
    groups <- c(groups, str_remove(str_remove(str_extract(bc_path[i], 
                                                          "pool_[0-9][0-9]_.+\\.[0-2]"), 
                                              "pool_"),
                                   "_CD.+atac"))
  }
  barcode_list_add[[i]] <- read.delim(bc_path[i], header = F, row.names = NULL)
  colnames(barcode_list_add[[i]]) <- c("barcode", "line")
  barcode_list_add[[i]]$line <- str_replace(barcode_list_add[[i]]$line, "-", "_")
}
groups <- str_replace(groups, "\\.2", "\\.6")
groups <- str_replace(groups, "\\.", "_")
names(barcode_list_add) <- groups
barcode_list <- c(barcode_list, barcode_list_add)


# 029
bc_path <- sort(list.files(path = "/nvmefs/scARC_Duan_018/Duan_project_029_RNA/Duan_029_human_only_demux_barcodes",
                           full.names = T, pattern = ".txt", recursive = T))
barcode_list_add <- vector(mode = "list", length = length(bc_path))
for (i in 1:length(bc_path)){
  if (i == 1) {
    groups <- c(str_extract(bc_path[i], "CW[0-9]+-[0|1|6]"))
  } else {
    groups <- c(groups, str_extract(bc_path[i], "CW[0-9]+-[0|1|6]"))
  }
  barcode_list_add[[i]] <- read_table(bc_path[i], col_names = FALSE)
  colnames(barcode_list_add[[i]]) <- c("barcode", "line")
}
names(barcode_list_add) <- str_remove(groups, "CW")
barcode_list <- c(barcode_list, barcode_list_add)

# 029 - 40201
bc_path <- sort(list.files(path = "/data/FASTQ/Duan_Project_029_40201/barcode_demux_output",
                           full.names = T, pattern = ".tsv", recursive = T))
barcode_list_add <- vector(mode = "list", length = length(bc_path))
for (i in 1:length(bc_path)){
  if (i == 1) {
    groups <- c(str_extract(bc_path[i], "[0-9]+_atac\\.[0|1|2]"))
  } else {
    groups <- c(groups, str_extract(bc_path[i], "[0-9]+_atac\\.[0|1|2]"))
  }
  barcode_list_add[[i]] <- read_table(bc_path[i], col_names = FALSE)
  colnames(barcode_list_add[[i]]) <- c("barcode", "line")
}
groups <- str_replace(groups, "\\.2", "\\.6")
groups <- str_replace(groups, "\\.", "_")
names(barcode_list_add) <- str_remove(groups, "atac_")
barcode_list <- c(barcode_list, barcode_list_add)

# 030
bc_path <- sort(list.files(path = "/data/FASTQ/Duan_Project_030/Duan_030_barcodes_by_individual/",
                           full.names = T, pattern = ".tsv", recursive = T))
barcode_list_add <- vector(mode = "list", length = length(bc_path))
for (i in 1:length(bc_path)){
  if (i == 1) {
    groups <- c(str_extract(bc_path[i], "[0-9]+_gex_[0|1|6]"))
  } else {
    groups <- c(groups, str_extract(bc_path[i], "[0-9]+_gex_[0|1|6]"))
  }
  barcode_list_add[[i]] <- read_table(bc_path[i], col_names = FALSE)
  colnames(barcode_list_add[[i]]) <- c("barcode", "line")
}
names(barcode_list_add) <- str_remove(groups, "gex_")
barcode_list <- c(barcode_list, barcode_list_add)

names(barcode_list) <- str_replace(names(barcode_list), "_", "-")

# make objects and match cells to barcodes ####
raw_obj_lst <- vector("list", length(pathlist))
transformed_lst <- vector("list", length(pathlist))
for (i in 1:length(pathlist)) {
  mat <- Read10X_h5(filename = pathlist[[i]])
  obj <- CreateSeuratObject(counts = mat$`Gene Expression`,
                            project = libnames[i])
  print(paste0(i, " number of genes: ", nrow(obj), 
               ", number of cells: ", ncol(obj)))
  obj$lib.ident <- libnames[i]
  obj <- PercentageFeatureSet(obj, pattern = c("^MT-"),
                              col.name = "percent.mt")
  raw_obj_lst[[i]] <- obj
  transformed_lst[[i]] <- SCTransform(obj, 
                                      vars.to.regress = "percent.mt",
                                      method = "glmGamPoi",
                                      variable.features.n = 8000,
                                      seed.use = 115,
                                      verbose = T)
}
names(raw_obj_lst) <- names(pathlist)
names(transformed_lst) <- names(pathlist)
save(transformed_lst, file = "transformed_lst_before_demux_bc.RData")
save(raw_obj_lst, file = "raw_obj_lst.RData")
rm(raw_obj_lst)

# remove cells not matched to human barcodes ####
transformed_lst_mapped_to_demux_bc <- vector("list", length(transformed_lst))
human_only_lst <- vector("list", length(transformed_lst))

# sep lines in each group's barcode list, then assign line by barcode
for (i in 1:length(barcode_list)){
  bcs <- barcode_list[[i]]
  gr <- names(barcode_list)[i]
  lines <- unique(bcs$line)
  cat("group: ", gr, "lines: ", lines, "\n")
  
  obj <- transformed_lst[[i]]
  obj$cell.line.ident <- "unmatched"
  for (j in 1:length(lines)){
    line_spec_barcodes <- barcode_list[[i]]$barcode[barcode_list[[i]]$line == lines[j]]
    obj$cell.line.ident[obj@assays$RNA@counts@Dimnames[[2]] %in% 
                          line_spec_barcodes] <- lines[j]
  }
  print(unique(obj$cell.line.ident))
  transformed_lst_mapped_to_demux_bc[[i]] <- obj
  human_only_lst[[i]] <- subset(obj, cell.line.ident != "unmatched")  
}

names(transformed_lst_mapped_to_demux_bc) <- names(transformed_lst)
names(human_only_lst) <- names(transformed_lst)

save(human_only_lst, file = "human_only_lst_removed_cells_unmatched_to_demuxed_bc.RData")
save(transformed_lst_mapped_to_demux_bc,
     file = "transformed_lst_matched_w_demux_bc.RData")
rm(transformed_lst, transformed_lst_mapped_to_demux_bc)

# remove repeated cell lines in 022 objects #### 
cleaned_lst <- human_only_lst
indxs_022 <- which(names(human_only_lst) %in% str_replace(groups_022, "_", "-"))
indxs_others <- which(!names(human_only_lst) %in% str_replace(groups_022, "_", "-"))

for (i in 1:length(indxs_022)) {
  obj_022 <- human_only_lst[[indxs_022[i]]]
  lines_022 <- as.character(unique(obj_022$cell.line.ident))
  for (j in 1:length(indxs_others)) { # for each library in the non-022 library list, check agains the currently library in 022
    obj_others <- human_only_lst[[indxs_others[j]]] # extract the object of a non-022 library
    lines_others <- as.character(unique(obj_others$cell.line.ident)) # extract the lines in this non-022 library
    
    obj_022$cell.line.ident[obj_022$cell.line.ident %in% lines_others] <- "redone" # redone lines in 022 will be removed later
    lines_in_022_obj <- unique(obj_022$cell.line.ident) #check the new lines after marking redone
    print(lines_in_022_obj)
    
    if (i == 1) { 
      redone_lines <- lines_022[lines_022 %in% lines_others] #the redone lines are just to keep track
      if (length(redone_lines) != 0) {
        other_libs_with_redone_lines <- unique(obj_others$lib.ident)
      }
    } else {
      redone_lines <- c(redone_lines, lines_022[lines_022 %in% lines_others])
      if (length(redone_lines) != 0) {
        other_libs_with_redone_lines <- c(other_libs_with_redone_lines,
                                          unique(obj_others$lib.ident))
      }
    }
  }
  if (length(lines_in_022_obj) == 1 & lines_in_022_obj == "redone") {
    cleaned_lst[[indxs_022[i]]] <- "removed" # delete the entire object in 022 if all the lines are redone
    print("all removed")
  } else if ("redone" %in% lines_in_022_obj) {
    print("partially removed")
    cleaned_lst[[indxs_022[i]]] <- subset(obj_022, cell.line.ident != "redone") # only keep the non-redone lines in this 022 object
  }
} 
lines_redone <- unique(redone_lines)
libs_with_redone_lines <- str_remove(unique(str_extract(unique(other_libs_with_redone_lines),
                                             "[0-9]+-")), "-")
write.table(lines_redone, file = "022_cell_lines_that_were_repeated.txt", quote = F,
            sep = "\t", row.names = F, col.names = F)
write.table(libs_with_redone_lines, file = "libraries_that_contains_repeated_lines.txt", quote = F,
            sep = "\t", row.names = F, col.names = F)

cleaned_lst$`22-0` <- NULL
cleaned_lst$`22-1` <- NULL
cleaned_lst$`22-6` <- NULL

save(cleaned_lst, file = "cleaned_lst_no_nonhuman_cells_no_redone_lines.Rdata")
rm(human_only_lst)

# do QC on each library first to reduce size ####
for (i in 1:length(cleaned_lst)) {
  obj <- cleaned_lst[[i]]
  obj <- subset(obj, nFeature_RNA < 7500 & nFeature_RNA > 300 &
                  nCount_RNA < 40000 & nCount_RNA > 500 &
                  percent.mt < 15)
  print(ncol(obj))
  cleaned_lst[[i]] <- obj
  rm(obj)
}

# prepare reciprocal PCA integration ####
# find anchors
features <- SelectIntegrationFeatures(object.list = cleaned_lst,
                                      nfeatures = 3000,
                                      fvf.nfeatures = 3000)
cleaned_lst <- PrepSCTIntegration(cleaned_lst, anchor.features = features)
cleaned_lst <- lapply(X = cleaned_lst, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(object.list = cleaned_lst,
                                  anchor.features = features,
                                  reference = c(1, 2, 3),
                                  reduction = "rpca",
                                  normalization.method = "SCT",
                                  scale = F, 
                                  dims = 1:50)
#Found 9036 anchors
#all.genes <- transformed_lst[[1]]@assays$RNA@counts@Dimnames[[1]]
save(anchors, file = "anchors_after_cleaned_lst_QC.RData")
rm(cleaned_lst)
# integrate ####
integrated <- IntegrateData(anchorset = anchors,
                            verbose = T)
#save(cleaned_lst, file = "final_lst_prepped_4_SCT_integration.RData")
#save(anchors, file = "anchors_4_integration.RData")
save(integrated, file = "integrated_018-030_RNAseq_obj_test_remove_10lines.RData")
