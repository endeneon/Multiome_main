# Chuxuan Li 05/15/2023
# make 018-029 combined Seurat object into scanpy annData class

# init ####
library(reticulate)
library(scater)
library(Seurat)
library(cowplot)
library(SeuratDisk)
#Sys.setenv(RETICULATE_PYTHON = "r-reticulate/bin/python")

load("./018-029_RNA_integrated_labeled_with_harmony.RData")
ad <- import("anndata", convert = F)
integrated_ad <- Convert(from = integrated_labeled, to = "anndata", 
                         filename = "integrated_labeled.h5ad")
