#Chuxuan Li 1/17/2022
#Using the original Seurat-integrated dataset in Analysis_RNAseq_v2, using monocle
#to re-count and analyze the differential expression

library(monocle)
library(Seurat)
library(future)
library(sctransform)
library(readr)
library(ggplot2)
library(future)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)
library(patchwork)
library(cowplot)
library(gplots)


# set threads and parallelization
plan("multisession", workers = 2)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 11474836480)

