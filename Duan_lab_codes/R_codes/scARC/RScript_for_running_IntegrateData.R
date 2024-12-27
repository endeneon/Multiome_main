# Chuxuan Li 05/14/2023
# script to run integrateData

# init ####
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(stringr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)
library(readxl)
library(readr)

library(future)

plan("multicore", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 161061273600)

load("./anchors_after_cleaned_lst_QC.RData")
integrated <- IntegrateData(anchorset = anchors,
                            verbose = T)
save(integrated, file = "integrated_018-030_RNAseq_obj_test_QC.RData")
