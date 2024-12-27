# Siwei 09 Feb 2024
# Rscript for counting npglut 1hr SNP-flanking peaks

# init ####
{
  library(Seurat)
  library(Signac)
  
  # library(EnsDb.Hsapiens.v86)
  library(GenomicFeatures)
  # library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  
  # library(edgeR)
  
  # library(RColorBrewer)
  
  library(stringr)
  library(future)
  
  # library(ggplot2)
  
  library(parallel)
  # library(doParallel)
  library(future)
  # library(foreach)
  # library()
  
  library(MASS)
}

# param #####
plan("multisession", workers = 8)
set.seed(42)
options(future.globals.maxSize = 229496729600)

setwd(".")
# load data ####

load("multiomic_obj_new_470K_cells_frag_file.RData")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/df_npglut_cells_62K.RData")
# run
load("peaks_inside_npglut_1hr_SNPs_500bp.RData")
load("peaks_inside_npglut_6hr_SNPs_500bp.RData")

# confirm all data have been properly loaded
print(length(ATAC_fragment))
print(length(npglut_barcodes))
print(length(peaks_flanking_npglut_1hr_SNPs))
print(length(peaks_flanking_npglut_6hr_SNPs))


# run
system.time({
  npglut_hr1_peaks_count <-
    FeatureMatrix(fragments = ATAC_fragment, 
                  features = peaks_flanking_npglut_1hr_SNPs,
                  process_n = 2000, 
                  cells = npglut_barcodes,
                  verbose = T)
})
print("1 hr counting done")
save(npglut_hr1_peaks_count,
     file = "npglut_hr1_peaks_count_FeatureMatrix_76_lines_12Feb2024.RData")
rm(npglut_hr1_peaks_count)
rm(peaks_flanking_npglut_1hr_SNPs)

# run
system.time({
  npglut_hr6_peaks_count <-
    FeatureMatrix(fragments = ATAC_fragment, 
                  features = peaks_flanking_npglut_6hr_SNPs,
                  process_n = 2000, 
                  cells = npglut_barcodes,
                  verbose = T)
})
print("6 hr counting done")
save(npglut_hr6_peaks_count,
     file = "npglut_hr6_peaks_count_FeatureMatrix_76_lines_12Feb2024.RData")
