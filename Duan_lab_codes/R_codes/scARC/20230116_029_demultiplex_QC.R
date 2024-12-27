# Chuxuan Li 1/16/2023
# 029 data demultiplex, remove nonhuman cells, QC

# init ####
library(Seurat)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
library(readr)

# read h5, remove rat cells and genes ####
h5list <- sort(list.files(path = "/data/FASTQ/Duan_Project_029/GRCh38_only",
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T,
                          include.dirs = F, 
                          full.names = T))
objlist <- vector(mode = "list", length = length(h5list))
libs <- str_extract(h5list, "[0-9]+-[0-6]")
cell_counts <- data.frame(library = libs,
                          total = rep(0, by = length(libs)),
                          rat_or_rabbit = rep(0, by = length(libs)),
                          nonrat_or_rabbit = rep(0, by = length(libs)),
                          mapped_to_human_barcodes = rep(0, by = length(libs)),
                          not_mapped_to_human_barcodes = rep(0, by = length(libs)))
load("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/removed_rat_RNAseq_list_029_mapped_to_hybrid.RData")
noratlist <- vector(mode = "list", length = length(h5list))
for (i in 1:length(objlist)){
  h5file <- Read10X_h5(filename = h5list[i])
  cat("h5: ", str_extract(string = h5list[i],
                          pattern = "[0-9]+-[0-6]"))
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(string = h5list[i],
                                                  pattern = "[0-9]+-[0-6]"))
  print(ncol(obj))
  cell_counts$total[i] <- ncol(obj)
  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list[i],
                                                  pattern = "-[0-6]"), "-"), "hr")
  print(unique(obj$time.ident))
  # assign rat identity
  humanbc <- colnames(cleanobj_lst[[i]])
  obj$rat.ident <- "rat_rabbit"
  obj$rat.ident[colnames(obj) %in% humanbc] <- "not_rat_or_rabbit"
  cell_counts$rat_or_rabbit[i] <- sum(obj$rat.ident == "rat_rabbit")
  cell_counts$nonrat_or_rabbit[i] <- sum(obj$rat.ident == "not_rat_or_rabbit")
  objlist[[i]] <- obj
  noratobj <- subset(obj, rat.ident == "not_rat_or_rabbit")
  noratlist[[i]] <- noratobj
}
save(objlist, file = "GRCh38_mapped_raw_list.RData")

# map to demuxed barcodes ####
setwd("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/Duan_029_human_only_demux_barcodes/")
pathlist <- sort(list.files(path = ".", full.names = T, pattern = ".txt", recursive = T))
barcode_list <- vector(mode = "list", length = length(pathlist))
for (i in 1:length(pathlist)){
  barcode_list[[i]] <- read_table(pathlist[i], col_names = FALSE)
  colnames(barcode_list[[i]]) <- c("barcode", "line")
}

# sep lines in each group's barcode list, then assign line by barcode
for (i in 1:length(barcode_list)){
  lines <- unique(barcode_list[[i]]$line)
  print(lines)
  noratlist[[i]]$cell.line.ident <- "unmatched"
  for (j in 1:length(lines)){
    line_spec_barcodes <- barcode_list[[i]]$barcode[barcode_list[[i]]$line == lines[j]]
    noratlist[[i]]$cell.line.ident[noratlist[[i]]@assays$RNA@counts@Dimnames[[2]] %in% 
                                     line_spec_barcodes] <- lines[j]
  }
  cell_counts$mapped_to_human_barcodes[i] <- sum(noratlist[[i]]$cell.line.ident != "unmatched")
  cell_counts$not_mapped_to_human_barcodes[i] <- sum(noratlist[[i]]$cell.line.ident == "unmatched")
}

setwd("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/")
write.table(cell_counts, file = "cell_counts_based_solely_on_GRCh38_mapped_data.csv", 
            sep = ",", quote = F, col.names = T, row.names = F)
save(noratlist, file = "GRCh38_mapped_removed_rat_cells_assigned_demux_bc_list.RData")

# remove cells not matched to cell line barcodes
human_only_lst <- vector(mode = "list", length = length(objlist))
for (i in 1:length(noratlist)){
  print(unique(noratlist[[i]]$cell.line.ident))
  human_only_lst[[i]] <- subset(noratlist[[i]], cell.line.ident != "unmatched")
  print(unique(human_only_lst[[i]]$cell.line.ident))
}

# QC ####
# merge without correcting anything first to see distribution
for (i in 1:length(human_only_lst)) {
  human_only_lst[[i]][["percent.mt"]] <- PercentageFeatureSet(human_only_lst[[i]], pattern = "^MT-")
}
for (i in 1:(length(human_only_lst) - 1)) {
  if (i == 1) {
    rough_merged <- merge(human_only_lst[[i]], human_only_lst[[i + 1]])
  } else {
    rough_merged <- merge(rough_merged, human_only_lst[[i + 1]])
  }
} 
# check nfeature, ncount, pct mt distribution
VlnPlot(rough_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
sum(rough_merged$percent.mt > 20) #1454
sum(rough_merged$nFeature_RNA > 8500) #322
sum(rough_merged$nCount_RNA > 50000) #293
sum(rough_merged$nFeature_RNA < 500) #338
sum(rough_merged$nCount_RNA < 500) #202

QCed_lst <- vector(mode = "list", length = length(human_only_lst))
for (i in 1:length(human_only_lst)) {
  QCed_lst[[i]] <- subset(human_only_lst[[i]], 
                          subset = nFeature_RNA > 500 & 
                            nFeature_RNA < 8500 & 
                            nCount_RNA > 500 &
                            nCount_RNA < 50000 &
                            percent.mt < 20)
}
save(QCed_lst, file = "GRCh38_mapped_after_QC_list.RData")
