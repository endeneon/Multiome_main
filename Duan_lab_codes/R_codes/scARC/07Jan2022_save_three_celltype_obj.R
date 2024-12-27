setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes")
load("02Dec2021_calc_da_peaks_motif.RData")

# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)

library(readr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(graphics)
library(stringr)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)


set.seed(1105)

three_celltype_obj <- subset(obj_complete, subset = broad.cell.type %in% 
                               c("GABA", "NPC", "glut"))

save("three_celltype_obj", file = "three_celltype_obj.RData")
