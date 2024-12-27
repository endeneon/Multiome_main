# Chuxuan Li 05/07/2023
# demultiplex, remove nonhuman, QC on 025 17 and 46 libraries

# init ####
library(Seurat)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
library(readr)
setwd("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/")

# read h5, remove rat cells and genes ####
h5list <- sort(list.files(path = "/data/FASTQ/Duan_Project_025_17_46",
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T,
                          include.dirs = F, 
                          full.names = T))
h5list <- h5list[str_detect(h5list, "hybrid", T)]
objlist <- vector(mode = "list", length = length(h5list))
libs <- str_extract(h5list, "[0-9]+-[0-6]")
cell_counts <- data.frame(library = libs,
                          total = rep(0, by = length(libs)),
                          rat = rep(0, by = length(libs)),
                          nonrat = rep(0, by = length(libs)),
                          human = rep(0, by = length(libs)),
                          nonhuman = rep(0, by = length(libs)))
load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part1_hybrid_genome/removed_rat_list_025_17_46_mapped_to_hybrid.RData")
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
  obj$rat.ident <- "rat"
  obj$rat.ident[colnames(obj) %in% humanbc] <- "human"
  cell_counts$rat[i] <- sum(obj$rat.ident == "rat")
  cell_counts$nonrat[i] <- sum(obj$rat.ident == "human")
  objlist[[i]] <- obj
  noratobj <- subset(obj, rat.ident == "human")
  noratlist[[i]] <- noratobj
}
save(objlist, file = "GRCh38_mapped_raw_list_17_46.RData")

# map to demuxed barcodes ####
setwd("/data/FASTQ/Duan_Project_025_17_46/barcode_demux_output/common_barcodes_17")
pathlist <- sort(list.files(path = ".", full.names = T, pattern = ".best.tsv", recursive = T))
barcode_list <- vector(mode = "list", length = 6)
for (i in 1:length(pathlist)){
  barcode_list[[i]] <- read.delim(pathlist[i], header = F, row.names = NULL)
  colnames(barcode_list[[i]]) <- c("barcode", "line")
}
names(barcode_list) <- paste0("17_", str_remove(str_extract(pathlist, "[0-2].best"), "\\.best"))
setwd("/data/FASTQ/Duan_Project_025_17_46/barcode_demux_output/common_barcodes_46")
pathlist <- sort(list.files(path = ".", full.names = T, pattern = ".best.tsv", recursive = T))
j = 1
for (i in 4:6){
  barcode_list[[i]] <- read.delim(pathlist[j], header = F, row.names = NULL)
  colnames(barcode_list[[i]]) <- c("barcode", "line")
  j = j + 1
}
names(barcode_list)[4:6] <- paste0("46_", str_remove(str_extract(pathlist, "[0-2].best"), "\\.best"))

# sep lines in each group's barcode list, then assign line by barcode
names(noratlist) <- str_extract(h5list, "[0-9]+-[0-6]")
for (i in 1:length(barcode_list)){
  lines <- unique(barcode_list[[i]]$line)
  print(lines)
  noratlist[[i]]$cell.line.ident <- "unmatched"
  for (j in 1:length(lines)){
    line_spec_barcodes <- barcode_list[[i]]$barcode[barcode_list[[i]]$line == lines[j]]
    noratlist[[i]]$cell.line.ident[noratlist[[i]]@assays$RNA@counts@Dimnames[[2]] %in% 
                                     line_spec_barcodes] <- lines[j]
  }
  cell_counts$human[i] <- sum(noratlist[[i]]$cell.line.ident != "unmatched")
  cell_counts$nonhuman[i] <- sum(noratlist[[i]]$cell.line.ident == "unmatched")
}

setwd("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/")
write.table(cell_counts, file = "cell_counts_based_solely_on_GRCh38_mapped_data_17_46.csv", 
            sep = ",", quote = F, col.names = T, row.names = F)
save(noratlist, file = "GRCh38_mapped_removed_rat_cells_assigned_demux_bc_list_17_46.RData")
# remove cells not matched to cell line barcodes
human_only_lst <- vector(mode = "list", length = length(objlist))
for (i in 1:length(objlist)){
  print(unique(noratlist[[i]]$cell.line.ident))
  human_only_lst[[i]] <- subset(noratlist[[i]], cell.line.ident != "unmatched")
  print(unique(human_only_lst[[i]]$cell.line.ident))
}
save(human_only_lst, file = "GRCh38_mapped_demux_matched_human_only_list_17_46.RData")

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
sum(rough_merged$percent.mt > 20) #2504
sum(rough_merged$nFeature_RNA > 7500) #434
sum(rough_merged$nCount_RNA > 40000) #267
sum(rough_merged$nFeature_RNA < 300) #116
sum(rough_merged$nCount_RNA < 500) #348

QCed_lst <- vector(mode = "list", length = length(human_only_lst))
for (i in 1:length(human_only_lst)) {
  QCed_lst[[i]] <- subset(human_only_lst[[i]], 
                          subset = nFeature_RNA > 300 & 
                            nFeature_RNA < 7500 & 
                            nCount_RNA > 500 &
                            nCount_RNA < 40000 &
                            percent.mt < 20)
}
save(QCed_lst, file = "GRCh38_mapped_after_QC_list_17_46.RData")
