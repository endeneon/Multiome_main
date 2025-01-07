library(tidyverse)
library(ArchR)

set.seed(1)
setwd("/scratch/midway3/zichengwang/arrow_files")
addArchRThreads(threads = 12) 
addArchRGenome("hg38")

# List all input files from 10X 
batch_folders <- list.files("/project/xinhe/zicheng/scATAC_seq")
batch_df <- strsplit(batch_folders, "-|_") %>% do.call(rbind, .) %>% as.data.frame()
colnames(batch_df) <- c("Batch", "Lib", "Time")
saveRDS(batch_df, "/project/xinhe/zicheng/neuron_stimulation/data/batch_df.rds")

# Read barcodes from each batch
Generate_barcodes <- function(each_row) {
  Batch <- each_row[1]
  Lib <- each_row[2]
  Time <- each_row[3]
  
  barcode_dir <- paste0("/project/xinhe/zicheng/neuron_stimulation/data/barcodes/Barcodes_Cleaned_Oct23/Duan_",
                        Batch, "/lib_", Lib, "/")
  barcode_txts <- list.files(barcode_dir, pattern = paste0("*", Time, "hr.*\\.txt$"))
  stopifnot(length(barcode_txts) > 0)
  barcodes <- sapply(barcode_txts, function(x) scan(paste0(barcode_dir, x), character())) %>% unlist()
  names(barcodes) <- NULL
  
  return(barcodes)
}

barcode_list <- apply(batch_df, 1, Generate_barcodes)
names(barcode_list) <- paste0("Batch_", batch_df$Batch, "_", batch_df$Lib, "_", batch_df$Time)

# List all input files from 10X
inputFiles <- paste0("/project/xinhe/zicheng/scATAC_seq/", batch_df$Batch, "_", batch_df$Lib, "-", batch_df$Time, "/atac_fragments.tsv.gz")
names(inputFiles) <- paste0("Batch_", batch_df$Batch, "_", batch_df$Lib, "_", batch_df$Time)

# Create Arrow files
ArrowFiles <- createArrowFiles(inputFiles = inputFiles, validBarcodes = barcode_list)
ArrowFiles


set.seed(1)
setwd("/scratch/midway3/zichengwang/arrow_files")
addArchRThreads(threads = 16) 

arrow_dir <- "/scratch/midway3/zichengwang/arrow_files"
batch_df <- readRDS("/project/xinhe/zicheng/neuron_stimulation/data/batch_df.rds")

# Extract meta data
extract_meta_per_sample <- function(sample_name, Batch, Lib, Time, barcodes) {
  split_name <- unlist(strsplit(sample_name, "_(?=[^_]*$)", perl = TRUE))
  sample_id <- split_name[1]
  cell_type <- split_name[2]
  meta_per_sample <- data.frame(batch_barcodes = paste0("Batch_", Batch, "_", Lib, "_", Time, "#", barcodes[[sample_name]]),
                                Time = paste0(Time, "hr"), sample_id = sample_id, cell_type = cell_type)
  return(meta_per_sample)
}

extract_meta <- function(each_row) {
  names(each_row) <- NULL
  Batch <- each_row[1]
  Lib <- each_row[2]
  Time <- each_row[3]
  
  barcode_dir <- paste0("/project/xinhe/zicheng/neuron_stimulation/data/barcodes/Barcodes_Cleaned_Oct23/Duan_",
                        Batch, "/lib_", Lib, "/")
  barcode_txts <- list.files(barcode_dir, pattern = paste0("*", Time, "hr.*\\.txt$"))
  stopifnot(length(barcode_txts) > 0)
  barcodes <- sapply(barcode_txts, function(x) scan(paste0(barcode_dir, x), character()))
  names(barcodes) <- gsub("_barcodes.txt", "", names(barcodes)) %>% gsub(paste0(Time, "hr_"), "", .) %>% 
    gsub("Duan030_", "", .) %>% gsub("029_", "", .) %>% gsub("NEFM_pos_glut", "npglut", .) %>% gsub("NEFM_neg_glut", "nmglut", .)
  
  meta_per_batch <- lapply(names(barcodes), function(x) extract_meta_per_sample(x, Batch, Lib, Time, barcodes)) %>% bind_rows()
  
  return(meta_per_batch)
}

barcode_by_time <- lapply(c("0", "1", "6"), function(x) apply(batch_df[batch_df$Time == x, ], 1, extract_meta) %>% bind_rows()) %>% bind_rows()

arrow_files <- paste0("Batch_", batch_df$Batch, "_", batch_df$Lib, "_", batch_df$Time, ".arrow")

addArchRGenome("hg38")

# Integrate all data into one ArchR project
ArchR_AllCells <- ArchRProject(
  ArrowFiles = arrow_files, 
  outputDirectory = "../ArchR_AllCells",
  copyArrows = TRUE
)

barcode_by_time <- barcode_by_time %>% filter(batch_barcodes %in% ArchR_AllCells$cellNames)

# Add metadata to ArchR project
ArchR_AllCells <- addCellColData(ArchRProj = ArchR_AllCells, data = barcode_by_time$Time,
                                 cells = barcode_by_time$batch_barcodes, name = "Time")

ArchR_AllCells <- addCellColData(ArchRProj = ArchR_AllCells, data = barcode_by_time$sample_id,
                                 cells = barcode_by_time$batch_barcodes, name = "sample_id")

ArchR_AllCells <- addCellColData(ArchRProj = ArchR_AllCells, data = barcode_by_time$cell_type,
                                 cells = barcode_by_time$batch_barcodes, name = "cell_type")

ArchR_AllCells$cell_cluster <- paste0(ArchR_AllCells$cell_type, "__", ArchR_AllCells$Time)
ArchR_AllCells <- ArchR_AllCells[!is.na(ArchR_AllCells$sample_id), ]

ArchR_AllCells <- addGroupCoverages(ArchRProj = ArchR_AllCells,
                                    sampleLabels = "sample_id",
                                    groupBy = "cell_cluster",
                                    maxReplicates = 100)

ArchR_AllCells <- addReproduciblePeakSet(ArchRProj = ArchR_AllCells,
                                         groupBy = "cell_cluster",
                                         excludeChr = c("chrX", "chrY"),
                                         maxPeaks = 250000,
                                         cutOff = 0.05)

ArchR_AllCells <- addPeakMatrix(ArchR_AllCells)

saveArchRProject(ArchRProj = ArchR_AllCells, load = FALSE)

################### Prepare RNA Matrix for ArchR
library(ArchR)
library(parallel)

source("/project/xinhe/zicheng/bin/R_functions.R")

batch_folders <- list.files("/project/xinhe/zicheng/scATAC_seq")
batch_df <- strsplit(batch_folders, "-|_") %>% do.call(rbind, .) %>% as.data.frame()
colnames(batch_df) <- c("Batch", "Lib", "Time")
inputFiles <- paste0("/project/xinhe/zicheng/scATAC_seq/", batch_df$Batch, "_", batch_df$Lib, "-", batch_df$Time, "/filtered_feature_bc_matrix.h5")
names(inputFiles) <- paste0("Batch_", batch_df$Batch, "_", batch_df$Lib, "_", batch_df$Time)

ArchR_w_scRNA <- loadArchRProject("/scratch/midway3/zichengwang/ArchR_w_scRNA/")

import_10Xh5 <- function(file_name, batch_name) {
  batch_barcodes <- ArchR_w_scRNA$cellNames[ArchR_w_scRNA$Sample == batch_name]
  
  batch_SE <- import10xFeatureMatrix(
    input = c(file_name),
    names = c(batch_name),
    featureType = "Gene Expression"
  )
  
  shared_barcodes <- intersect(batch_barcodes, colnames(batch_SE))
  batch_SE <- batch_SE[, shared_barcodes]
  
  return(batch_SE)
}

SE_list <- lapply(1:72, function(x) import_10Xh5(inputFiles[x], names(inputFiles)[x]))
names(SE_list) <- names(inputFiles)

scRNA_SE <- read10x_remove_rows(SE_list, names = names(SE_list))
saveRDS.gz(scRNA_SE, "/scratch/midway3/zichengwang/h5seurat/SE_scRNA_filtered.rds", 8)

############# Add expression matrix to ArchR
library(ArchR)
source("/project/xinhe/zicheng/bin/R_functions.R")

setwd("/scratch/midway3/zichengwang")
addArchRThreads(threads = 8) 

scRNA_SE <- readRDS.gz("/scratch/midway3/zichengwang/h5seurat/SE_scRNA_filtered.rds", 8)

ArchR_w_scRNA <- loadArchRProject("/scratch/midway3/zichengwang/ArchR_w_scRNA/")
ArchR_w_scRNA <- ArchR_w_scRNA[ArchR_w_scRNA$cellNames %in% colnames(scRNA_SE), ]

ArchR_w_scRNA <- addGeneExpressionMatrix(
  input = ArchR_w_scRNA,
  seRNA = scRNA_SE,
  excludeChr = c("chrM", "chrX", "chrY"),
  strictMatch = TRUE,
  force = TRUE,
  threads = 8
)

saveArchRProject(ArchR_w_scRNA)