# Chuxuan Li 10/28/2021
# Subset the Integrated and labeled RNAseq dataset to get barcodes by sub-cell
#types, by cell line, and by time, for SNP calling

# init
library(Seurat)
library(future)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)
library(patchwork)
library(cowplot)
library(gplots)

# set threads and parallelization
plan("multisession", workers = 6)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

# load data

