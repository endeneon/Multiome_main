# Chuxuan Li 04/18/2023
# test if running 20230324_ATACseq_afer_aggr_map_to_RNAseq.R with bash can help

library(Seurat)
library(Signac)
library(sctransform)+
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

plan("multisession", workers = 2)
set.seed(2001)
options(future.globals.maxSize = 5368709120)

# load data ####
# load part 1
library(Matrix)
counts1 <- Read10X("./018_to_029_combined/outs/filtered_feature_bc_matrix/part1/gzfiles/")
counts2 <- Read10X("./018_to_029_combined/outs/filtered_feature_bc_matrix/part2/gzfiles/")
counts3 <- Read10X("./018_to_029_combined/outs/filtered_feature_bc_matrix/part3/gzfiles/")
frag <- "./018_to_029_combined/outs/atac_fragments.tsv.gz"
counts <- counts1$Peaks + counts2$Peaks + counts3$Peaks

# load annotation
load("../R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")

ATAC <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = frag,
  annotation = ens_use
)
ATAC # 225934 features for 795150 cells
sum(rowSums(ATAC@data))
# 4185464619
save(ATAC, file = "018-029_ATAC_chromatin_assay.RData")

