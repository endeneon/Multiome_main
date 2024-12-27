# Siwei 02 Sept 2024
# Use Lexi's multiomic data of 28 Apr 2023 (new peakset)
# and replot ATACseq data w new peak set projected by RNAseq cell types


# init ####
{
  library(Seurat)
  library(Signac)
  library(sctransform)
  library(glmGamPoi)
  # library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
  # library(GenomicFeatures)
  library(AnnotationDbi)
  # library(BSgenome.Hsapiens.UCSC.hg38)
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
  
}

plan("multisession", workers = 4)
set.seed(2001)
options(future.globals.maxSize = 229496729600)

setwd("/nvmefs/scARC_Duan_018/018-029_combined_analysis")

# load
multiomic_obj_new_plot_umap <-
	readRDS(file = "multiomic_obj_w_gact_TSS.RData")

DefaultAssay(multiomic_obj_new_plot_umap) <- "ATAC"
multiomic_obj_new_plot_umap <-
  NucleosomeSignal(multiomic_obj_new_plot_umap)
saveRDS(multiomic_obj_new_plot_umap,
        file = "multiomic_obj_w_gact_TSS_NS.RDs")
q(save = "no")
