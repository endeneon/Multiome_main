# Siwei 17 Oct 2024

# init####
library(Seurat)
library(Signac)
library(sctransform)
library(glmGamPoi)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

library(patchwork)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)
library(stringr)
library(future)

plan("multisession", workers = 4)
set.seed(2001)
options(future.globals.maxSize = 229496729600)

# set working env
setwd("/nvmefs/scARC_Duan_018/018-029_combined_analysis")

# load
load("../Duan_project_024_ATAC/macs2_called_new_peaks.RData")

# load("multiomic_obj_new_plot_umap_w_gact.RData")