# Chuxuan Li Sept 27, 2021
# build a pipeline to separate data by cell line and do QC

# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)
library(future)
library(stringr)
# set threads and parallelization
plan("multisession", workers = 6)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

# 1. read in all the .h5 files, use only the RNA part
file_names <- list.files("")  # get a list of .h5 file names from the path that needs specified
n_line <- 0 # total number of lines analyzed in one pass, by group makes more sense 
n_group <- 0 # total number of groups analyzed in one pass

n_obj <- n_group * 3 # number of lines times 3 time points is the final number of files to work with

gex_lst <- vector(mode = "list", length = n_obj) # initialize a list to store all the gene expression objects

for (i in 1:length(file_names)){
  f <- file_lst[i]
  raw <- Read10X_h5(filename = f)
  raw_gex <- CreateSeuratObject(counts = raw$`Gene Expression`,
                            project = "group____")
  gex <- raw_gex
  gex_lst[[i]] <- gex
}

# get percent mitochondrial genes as a new column, then filter out cells with percent.mt > 10
for (i in 1:length(gex_lst)){
  g <- gex_lst[[i]]
  g <- PercentageFeatureSet(g, 
                            pattern = c("^MT-"),
                            col.name = "percent.mt")
  g <- subset(g, subset = percent.mt < 10)
  gex_lst[[i]] <- g
}

# SCTransform (slow step)
for (i in 1:length(gex_lst)){
  g <- gex_lst[[i]]
  g <- SCTransform(g,
                   vars.to.regress = "percent.mt", 
                   variable.features.n = 8000,
                   method = "glmGamPoi",
                   verbose = T)
  
  g <- RunPCA(g, 
              verbose = T)
  g <- RunUMAP(g,
               dims = 1:30,
               verbose = T)
  g <- FindNeighbors(g,
                     dims = 1:30,
                     verbose = T)
  g <- FindClusters(g,
                    verbose = T)
  
  g@meta.data$cell.line.ident <- "unknown"
  
  gex_lst[[i]] <- g
}

# read in cell line barcodes files
barcodes_names <- list.files("")  # get a list of barcode file names from the path that needs specified
n_sep <- n_line * 3 # number of lines times 3 time points is the total number of objects to separate into

barcodes_lst <- vector(mode = "list", length = n_sep)
# read in the barcode files, unlist them, then store in the list
for (i in 1:n_sep){
  f <- barcodes_names[[i]]
  bc <- read_csv(i, col_names = FALSE)
  bc <- unlist(bc)
  barcodes_lst[[i]] <- bc
}

# get cell line names in a list
cell_line_names <- vector(mode = "list", length = n_line)
for (i in 1:n_line){
  f <- barcodes_names[[i]]
  cl <- str_sub(f, 1:4) # depends on the name of the barcode file, extract the cell line names
}

# assign cell line identities 
for (i in 1:length(gex_lst)){
  g <- gex_lst[[i]]
  for (j in 1:length(barcodes_lst)){
    l <- barcodes_lst[j]
    g@meta.data$cell.line.ident[g@assays$SCT@counts@Dimnames[[2]] %in% l <- cell_line_names[[j]]]
    g@meta.data$cell.line.ident <- factor(g@meta.data$cell.line.ident)
  }
  gex_lst[[i]] <- g
}



