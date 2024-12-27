# Siwei 13 Feb 2024
# Rscript for counting npglut vs GABA 1=0hr SNP-flanking peaks

# init ####
{
  library(Seurat)
  library(Signac)
  
  # library(EnsDb.Hsapiens.v86)
  library(GenomicFeatures)
  # library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  
  library(stringr)
  library(future)
  
  # library(ggplot2)
  
  library(parallel)
  # library(doParallel)
  library(future)
  
  library(MASS)
}

# param #####data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCC
plan("multisession", workers = 8)
set.seed(42)
options(future.globals.maxSize = 229496729600)

setwd(".")
# load data ####

load("multiomic_obj_new_470K_cells_frag_file.RData")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/df_0hr_npglut_GABA_cells_50K.RData")
# run
load("peaks_inside_npglut_0hr_vs_GABA_SNPs_1000bp.RData")

# confirm all data have been properly loaded
print(head(ATAC_fragment))
print(head(npglut_GABA_0hr_barcodes))
print(head(peaks_flanking_npglut_0hr_SNPs))


# run
system.time({
  npglut_vs_GABA_hr0_peaks_count <-
    FeatureMatrix(fragments = ATAC_fragment, 
                  features = peaks_flanking_npglut_0hr_SNPs,
                  process_n = 2000, 
                  cells = npglut_GABA_0hr_barcodes,
                  verbose = T)
})
print("counting done")
save(npglut_vs_GABA_hr0_peaks_count,
     file = "npglut_vs_GABA_hr0_peaks_count_1000bp_FeatureMatrix_76_lines_12Feb2024.RData")
