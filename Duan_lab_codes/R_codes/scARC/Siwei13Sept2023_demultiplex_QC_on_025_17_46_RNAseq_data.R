# Siwei 13 Sept 2023
# Adapted from Chuxuan Li 05/04/2023
# map to demuxlet generated barcodes and QC for 025 17/46 RNAseq data

# init ####
library(Seurat)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
library(readr)

# need this line or R cannot get through proxy
options(download.file.method = "wget")


setwd("~/NVME/scARC_Duan_018/Duan_project_025_17_46_RNA/Analysis_pt2_mapped_to_Hg38_only_siwei")
# setwd("/nvmefs/scARC_Duan_018/Duan_project_030_RNA/Analysis_part2_mapped_to_Hg38_only/")

# read h5, remove rat cells and genes ####
h5list <- sort(list.files(path = "/data/FASTQ/Duan_Project_025_17_46",
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T,
                          include.dirs = F,
                          full.names = T))
h5list <- h5list[str_detect(h5list, "hg38_only", negate = F)]
objlist <- vector(mode = "list", length = length(h5list))
libs <- str_extract(h5list, "[0-9]+-[0-6]")
names(objlist) <- libs
cell_counts <- data.frame(library = libs,
                          total = rep(0, by = length(libs)),
                          mouse = rep(0, by = length(libs)),
                          nonmouse = rep(0, by = length(libs)),
                          human = rep(0, by = length(libs)),
                          nonhuman = rep(0, by = length(libs)))
load("../Analysis_pt1_hybrid_genome_Siwei/removed_mouse_list_025_17_46_mapped_to_hybrid.RData")
nomouselist <- vector(mode = "list", length = length(h5list))
names(nomouselist) <- libs

for (i in 1:length(objlist)) {
  h5file <- Read10X_h5(filename = h5list[i])
  cat("h5: ", str_extract(string = h5list[i],
                          pattern = "[0-9]+-[0-6]"))
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(string = h5list[i],
                                                  pattern = "[0-9]+-[0-6]"))
  print(ncol(obj))
  cell_counts$total[i] <- ncol(obj)
  # assign library identity
  obj$time.ident <-
    paste0(str_remove(str_extract(string = h5list[i],
                                  pattern = "-[0-6]"), "-"),
           "hr")
  print(unique(obj$time.ident))
  # assign mouse identity
  humanbc <- colnames(cleanobj_lst[[i]])
  obj$mouse.ident <- "mouse"
  obj$mouse.ident[colnames(obj) %in% humanbc] <- "human"
  cell_counts$mouse[i] <- sum(obj$mouse.ident == "mouse")
  cell_counts$nonmouse[i] <- sum(obj$mouse.ident == "human")
  objlist[[i]] <- obj
  nomouseobj <- subset(obj, mouse.ident == "human")
  nomouselist[[i]] <- nomouseobj
}
save(objlist, file = "GRCh38_mapped_raw_list_with_mouse_ident.RData")

# map to demuxed barcodes ####
# setwd("/data/FASTQ/Duan_Project_030/Duan_030_barcodes_by_individual")
pathlist <-
  sort(list.files(path = "/home/zhangs3/Data/FASTQ/Duan_Project_025_17_46/barcode_demux_output/common_barcodes",
                  full.names = T,
                  pattern = ".tsv",
                  recursive = T))
barcode_list <-
  vector(mode = "list",
         length = length(pathlist))

for (i in 1:length(pathlist)) {
  barcode_list[[i]] <-
    read.delim(pathlist[i],
               header = F,
               row.names = NULL)
  colnames(barcode_list[[i]]) <-
    c("barcode", "line")
}

names(barcode_list) <-
  c("Duan-025-17-0hr",
    "Duan-025-17-1hr",
    "Duan-025-17-6hr",
    "Duan-025-46-0hr",
    "Duan-025-46-1hr",
    "Duan-025-46-6hr")
#   str_remove(str_extract(pathlist, "[0-9]+_gex_[0|1|6]hr"), "gex_")

# sep lines in each group's barcode list, then assign line by barcode
for (i in 1:length(barcode_list)) {
  lines <- unique(barcode_list[[i]]$line)
  #print(lines)
  print(unique(nomouselist[[i]]$orig.ident))
  print(names(barcode_list)[i])
  nomouselist[[i]]$cell.line.ident <- "unmatched"
  for (j in 1:length(lines)) {
    line_spec_barcodes <- barcode_list[[i]]$barcode[barcode_list[[i]]$line == lines[j]]
    nomouselist[[i]]$cell.line.ident[nomouselist[[i]]@assays$RNA@counts@Dimnames[[2]] %in%
                                     line_spec_barcodes] <- lines[j]
  }
  cell_counts$human[i] <- sum(nomouselist[[i]]$cell.line.ident != "unmatched")
  cell_counts$nonhuman[i] <- sum(nomouselist[[i]]$cell.line.ident == "unmatched")
}

setwd("/nvmefs/scARC_Duan_018/Duan_project_030_RNA/Analysis_part2_mapped_to_Hg38_only/")
write.table(cell_counts,
            file = "cell_counts_based_solely_on_GRCh38_mapped_data.csv",
            sep = ",",
            quote = F,
            col.names = T,
            row.names = F)

save(nomouselist,
     file = "GRCh38_mapped_removed_mouse_cells_assigned_demux_bc_list.RData")

# remove cells not matched to cell line barcodes
human_only_lst <- vector(mode = "list",
                         length = length(nomouselist))
for (i in 1:length(nomouselist)) {
  print(unique(nomouselist[[i]]$cell.line.ident))
  human_only_lst[[i]] <- subset(nomouselist[[i]], cell.line.ident != "unmatched")
  print(unique(human_only_lst[[i]]$cell.line.ident))
}
saveRDS(human_only_lst, file = "GRCh38_mapped_demux_matched_human_only_list.RData")

# QC ####
# merge without correcting anything first to see distribution
for (i in 1:length(human_only_lst)) {
  human_only_lst[[i]][["percent.mt"]] <-
    PercentageFeatureSet(human_only_lst[[i]],
                         pattern = "^MT-")
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
sum(rough_merged$percent.mt > 15) # 7997
sum(rough_merged$nFeature_RNA > 7500) # 434
sum(rough_merged$nCount_RNA > 40000) # 269
sum(rough_merged$nFeature_RNA < 300) # 116
sum(rough_merged$nCount_RNA < 500) # 348

QCed_lst <-
  vector(mode = "list",
         length = length(human_only_lst))
for (i in 1:length(human_only_lst)) {
  QCed_lst[[i]] <-
    subset(human_only_lst[[i]],
           subset = nFeature_RNA > 300 &
             nFeature_RNA < 7500 &
             nCount_RNA > 500 &
             nCount_RNA < 40000 &
             percent.mt < 20)
}
saveRDS(QCed_lst, file = "GRCh38_mapped_after_QC_list.RData")


# test convert count to tpm
# install.packages("tiledbsoma",
#                  repos = c('https://tiledb-inc.r-universe.dev',
#                            'https://cloud.r-project.org'))
install.packages("cellxgene.census",
                 repos = c('https://chanzuckerberg.r-universe.dev',
                           'https://cloud.r-project.org'))



QCed_lst <-
  readRDS(file = "GRCh38_mapped_after_QC_list.RData")

test_seurat_obj <-
  QCed_lst[[1]]
