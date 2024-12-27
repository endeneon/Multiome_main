# Siwei 14 Feb 2022
# Run permutation on the peaks for Xin's lab
# to see the peak quality

# init
library(regioneR)
library(readr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)

# test if the masked version of hg38 is available
getGenomeAndMask("hg38")
library(BiocManager)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38") # no mask
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
# get human genome and apply hg38 blacklist mask
blacklisted_region <- toGRanges("~/Data/Databases/Genomes/hg38/hg38_blacklisted_regions.bed",
                                genome = "hg38")
hg38_genome <- getGenomeAndMask("hg38",
                                mask = blacklisted_region)

chrom_ranges <- read_delim("~/Data/Databases/Genomes/hg38/hg38.chrom.sizes", 
                           delim = "\t", escape_double = FALSE, 
                           col_names = FALSE, trim_ws = TRUE)
chrom_ranges$X3 <- 1
chrom_ranges <- chrom_ranges[1:24, c(1, 3, 2)]
chrom_ranges <- as.data.frame(chrom_ranges)
chrom_ranges <- toGRanges(chrom_ranges, genome = "hg38")

peaks_original <- toGRanges("peaks_for_He_Lab/GABA_0hr_specific_peaks.bed",
                            genome = "hg38")
set.seed(42)

# run permutation test
rand_region_hg38 <- function(x) {
  randomizeRegions(x, 
                   genome = "hg38",
                   mask = blacklisted_region)
}

pt_results <- permTest(A = peaks_original,
                       randomize.function = randomizeRegions,
                       genome = "hg38",
                       B = chrom_ranges,
                       universe = chrom_ranges,
                       evaluate.function = numOverlaps,
                       alternative = "less",
                       n = 50,
                       verbose = T)

pt_results <- permTest(A = peaks_original,
                       randomize.function = randomizeRegions,
                       genome = "hg38",
                       per.chorosome = T,
                       B = randomizeRegions(A = peaks_original,
                                            genome = "hg38",
                                            per.chromosome = T),
                       universe = chrom_ranges,
                       evaluate.function = numOverlaps,
                       mask = blacklisted_region,
                       # alternative = "greater",
                       n = 500,
                       force.parallel = T,
                       mc.set.seed = F,
                       mc.cores = 8,
                       count.once = T,
                       verbose = T)

pt_results <- permTest(A = randomizeRegions(A = peaks_original,
                                            genome = "hg38",
                                            per.chromosome = T),
                       randomize.function = randomizeRegions,
                       genome = "hg38",
                       per.chorosome = T,
                       B = peaks_original,
                       universe = chrom_ranges,
                       evaluate.function = meanDistance,
                       mask = blacklisted_region,
                       allow.overlaps = F,
                       # alternative = "greater",
                       n = 50,
                       force.parallel = T,
                       mc.set.seed = F,
                       mc.cores = 8,
                       count.once = T,
                       verbose = T)

pt_results <- permTest(A = sample(peaks_original, 
                                  size = 2000, replace = F),
                       randomize.function = randomizeRegions,
                       genome = "hg38",
                       per.chorosome = T,
                       B = peaks_original,
                       universe = chrom_ranges,
                       evaluate.function = numOverlaps,
                       mask = blacklisted_region,
                       allow.overlaps = F,
                       # alternative = "greater",
                       n = 50,
                       force.parallel = T,
                       mc.set.seed = F,
                       mc.cores = 8,
                       count.once = T,
                       verbose = T)
pt_results
summary(pt_results)
plot(pt_results)
# B = toGRanges("/home/zhangs3/Data/Databases/Genomes/hg38/gencode.v35.genes.bed",
#               genome = "hg38"),




######


load_file_list <- list.files(path = "./peaks_for_He_Lab/",
                             pattern = "*.bed")
save_file_list <- str_replace_all(string = load_file_list,
                                  pattern = "bed",
                                  replacement = "png")

dev.off()
i <- 1
for (i in 1:length(load_file_list)) {
  try({
  print(load_file_list[i])
  peaks_original <- toGRanges(paste("peaks_for_He_Lab/", load_file_list[i], sep = ""),
                                    genome = "hg38")
                              
  pt_overlap <- overlapPermTest(sample(peaks_original, 
                                       size = 400, replace = F),
                                randomizeRegions(A = peaks_original,
                                                 genome = "hg38"),
                                ntimes = 2000,
                                genome = "hg38",
                                count.once = T,
                                verbose = T,
                                force.parallel = T,
                                mc.set.seed = F,
                                mc.cores = 8)
  
    
  })
  
  png(filename = paste("permutation_PNG/", save_file_list[i], sep = ""),
      bg = "white", 
      width = 640, height = 480)
  plot(pt_overlap)
  
  dev.off()
}


# pt_overlap <- overlapPermTest(peaks_original, 
#                               randomizeRegions(A = peaks_original,
#                                                genome = "hg38"),
#                               ntimes = 100,
#                               genome = "hg38",
#                               count.once = T,
#                               verbose = T,
#                               force.parallel = T,
#                               mc.set.seed = F,
#                               mc.cores = 8)
pt_overlap
plot(pt_overlap)
