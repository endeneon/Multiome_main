# MACS2 call peak from Rscript
# Siwei 07 Apr 2022
# use 18 lines, Glut cells at 0 hr to check 
# rs2027349-c1orf54 connections

# need chr1 only

#init

library(Seurat)
library(SeuratWrappers)
library(Signac)
library(future)

library(remotes)
library(monocle3)
# library(monocle)
library(cicero)

library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(stringr)

library(ggplot2)
library(patchwork)

plan("multisession", workers = 4)
options(future.globals.maxSize = 100 * 1024^3)


# setwd
setwd("/nvmefs/scARC_Duan_018/R_peak_permutation")

# load data
load("chr1_18_samples.RData")

# call peaks use MACS2
MACS2_peaks <-
  CallPeaks(object = chr1_glut_ATAC_new,
            macs2.path = "~/Data/Anaconda3-envs/signac/bin/macs2",
            extsize = 150,
            shift = -73,
            additional.args = '-q 0.001 -g hs')
# save data
save.image()
