# Siwei 09 Feb 2024
# Calc chromatin accessibility correlation using Yifan's code
# use peaks from npglut 1hr and 6 hr to compare with 0 hr

# init ####
{
  library(Seurat)
  library(Signac)

  library(EnsDb.Hsapiens.v86)
  library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  
  library(edgeR)
  
  library(RColorBrewer)

  library(stringr)
  library(future)
  
  library(ggplot2)

  library(parallel)
  library(doParallel)
  library(future)
  library(foreach)
  # library()

  library(MASS)
}

# param #####
plan("multisession", workers = 8)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# load lexi's pre-made 500 bp peak data set ####
load("/nvmefs/scARC_Duan_018/018-029_combined_analysis/df_Lexi_new_peak_set_500bp.RData")

## use npglut_6hr only
peaks_500_width <-
  peaks_500_width[str_detect(string = names(peaks_500_width),
                             pattern = "npglut_6hr",
                             negate = F)]

## merge as one df
df_master_peak_set <-
  do.call(rbind,
          args = peaks_500_width)
sum(is.na(df_master_peak_set$chr))
unique(df_master_peak_set$chr)
df_master_peak_set <-
  df_master_peak_set[str_detect(string = df_master_peak_set$chr,
                                pattern = "^chr[1-9].*"), ]
length(sort(unique(df_master_peak_set$chr)))
## re-sort by chr names
df_master_peak_set$chr <-
  factor(df_master_peak_set$chr,
         levels = sort(unique(df_master_peak_set$chr)))

## split by chr name
all_peaks_by_chr <-
  split(x = df_master_peak_set,
        f = df_master_peak_set$chr)
names(all_peaks_by_chr)

discontinuous_peaks_list_by_chr <-
  vector(mode = "list",
         length = length(all_peaks_by_chr))
names(discontinuous_peaks_list_by_chr) <-
  names(all_peaks_by_chr)

## register a parallel backend for foreach #####
cluster_parallel <-
  makeCluster(8)
registerDoParallel(cluster_parallel)

## run the foreach for var i ####
discontinuous_peaks_list_by_chr <-
  foreach(i = 1:length(all_peaks_by_chr)) %dopar% {
    peaks_by_chr <-
      all_peaks_by_chr[[i]]
    peaks_by_chr <-
      peaks_by_chr[order(peaks_by_chr$start), ]
    for (j in 1:(nrow(peaks_by_chr))) {
      print(paste0("chr",
                   i,
                   ", peak ",
                   j))
      current_peak <-
        peaks_by_chr[j, ]
      if (j == 1) {
        previous_peak <- current_peak
        df_combined_peak <- current_peak
      } else {
        if (previous_peak$end <= current_peak$start) {
          df_combined_peak <-
            as.data.frame(rbind(df_combined_peak,
                                current_peak))
          previous_peak <- current_peak
        } else {
          if (previous_peak$score < current_peak$score) {
            df_combined_peak[nrow(df_combined_peak), ] <-
              unlist(current_peak)
            previous_peak <- current_peak
          }
        }
      }
    }
    return(df_combined_peak)
  }
stopCluster(cluster_parallel)
# 
# for (i in 1:length(all_peaks_by_chr)) {
#   peaks_by_chr <-
#     all_peaks_by_chr[[i]]
#   peaks_by_chr <-
#     peaks_by_chr[order(peaks_by_chr$start), ]
#   for (j in 1:(nrow(peaks_by_chr))) {
#     print(paste0("chr",
#                  i,
#                  ", peak ",
#                  j))
#     current_peak <-
#       peaks_by_chr[j, ]
#     if (j == 1) {
#       previous_peak <- current_peak
#       df_combined_peak <- current_peak
#     } else {
#       if (previous_peak$end <= current_peak$start) {
#         df_combined_peak <-
#           as.data.frame(rbind(df_combined_peak,
#                               current_peak))
#         previous_peak <- current_peak
#       } else {
#         if (previous_peak$score < current_peak$score) {
#           df_combined_peak[nrow(df_combined_peak), ] <-
#             unlist(current_peak)
#           previous_peak <- current_peak
#         }
#       }
#     }
#   }
#   discontinuous_peaks_list_by_chr[[i]] <-
#     df_combined_peak
# }
names(discontinuous_peaks_list_by_chr) <-
  names(all_peaks_by_chr)

full_discontinuous_peaks_list <-
  do.call(rbind,
          args = discontinuous_peaks_list_by_chr)
sum(is.na(full_discontinuous_peaks_list$chr))

df_2_write <-
  full_discontinuous_peaks_list
df_2_write$chr <-
  str_split(string = rownames(df_2_write),
            pattern = "\\.",
            simplify = T)[, 1]
unique(df_2_write$chr)
# df_2_write <-
#   df_2_write[, c(7, 2:6)]
df_2_write$start <-
  as.integer(df_2_write$start)
df_2_write$end <-
  as.integer(df_2_write$end)
df_2_write$strand <- '*'

GRanges_all_peaks <-
  makeGRangesFromDataFrame(df = df_2_write,
                           keep.extra.columns = T, 
                           na.rm = T)

save(GRanges_all_peaks,
     file = "Granges_npglut_6hr_peaks_non_overlapping_500bp.RData")

# load all SNPs ####
load("~/NVME/scARC_Duan_018/R_ASoC/Sum_100_lines_21Dec2023.RData")

## extract the lists w ASoC SNPs only ([[3]])
vcf_lst_ASoC_only <-
  vector(mode = "list",
         length = length(master_vcf_list))

for (i in 1:length(vcf_lst_ASoC_only)) {
  print(i)
  vcf_lst_ASoC_only[[i]] <-
    master_vcf_list[[i]][[3]]
}
names(vcf_lst_ASoC_only) <-
  names(master_vcf_list)

## extract the lists w all SNPs ([[1]])
vcf_lst_all_SNPs <-
  vector(mode = "list",
         length = length(master_vcf_list))

for (i in 1:length(vcf_lst_all_SNPs)) {
  print(i)
  vcf_lst_all_SNPs[[i]] <-
    master_vcf_list[[i]][[1]]
}
names(vcf_lst_all_SNPs) <-
  names(master_vcf_list)

## 1/6 hr specific criteria: FDR < 0.05 in 1/6 hr and #####
## binomP > 0.05 in 0hr
names(vcf_lst_all_SNPs)

## get npglut 6 hr specific SNPs ####
npglut_6hr_ASoC_SNPs <-
  vcf_lst_all_SNPs[['6hr_npglut']][vcf_lst_all_SNPs[['6hr_npglut']]$FDR < 0.05, ]
npglut_0hr_pBinom_lessthan_005_SNPs <-
  vcf_lst_all_SNPs[['0hr_npglut']][vcf_lst_all_SNPs[['0hr_npglut']]$pVal < 0.05, ]

npglut_6hr_specific_snps <-
  npglut_6hr_ASoC_SNPs[!(npglut_6hr_ASoC_SNPs$ID %in% npglut_0hr_pBinom_lessthan_005_SNPs$ID), ]
npglut_6hr_specific_snps <-
  npglut_6hr_specific_snps[!duplicated(npglut_6hr_specific_snps$ID), ]

npglut_6hr_specific_interval <-
  npglut_6hr_specific_snps[, c(1, 2, 2, 3, 15)]
npglut_6hr_specific_interval$POS.1 <-
  as.numeric(npglut_6hr_specific_interval$POS)
npglut_6hr_specific_interval$POS <-
  as.numeric(npglut_6hr_specific_interval$POS) - 1
npglut_6hr_specific_interval$strand <- "*"

colnames(npglut_6hr_specific_interval)[1:3] <-
  c("chr", "start", "end")
npglut_6hr_specific_interval <-
  makeGRangesFromDataFrame(df = npglut_6hr_specific_interval,
                           keep.extra.columns = T,
                           na.rm = T)

peaks_flanking_npglut_6hr_SNPs <-
  GRanges_all_peaks[GenomicRanges::nearest(x = npglut_6hr_specific_interval,
                                           subject = GRanges_all_peaks,
                                           ignore.strand = T,
                                           select = "arbitrary")]
peaks_flanking_npglut_6hr_SNPs_distance <-
  GenomicRanges::distanceToNearest(x = npglut_6hr_specific_interval,
                                   subject = peaks_flanking_npglut_6hr_SNPs,
                                   ignore.strand = T)
sum(peaks_flanking_npglut_6hr_SNPs_distance@elementMetadata$distance == 0)
# hist(peaks_flanking_npglut_1hr_SNPs_distance@elementMetadata$distance)

peaks_flanking_npglut_6hr_SNPs <-
  peaks_flanking_npglut_6hr_SNPs[peaks_flanking_npglut_6hr_SNPs_distance@elementMetadata$distance == 0]
peaks_flanking_npglut_6hr_SNPs <-
  subsetByOverlaps(x = peaks_flanking_npglut_6hr_SNPs,
                   ranges = blacklist_hg38_unified,
                   invert = T)
peaks_flanking_npglut_6hr_SNPs <-
  peaks_flanking_npglut_6hr_SNPs[!duplicated(peaks_flanking_npglut_6hr_SNPs@elementMetadata@listData$peak_ID)]
save(peaks_flanking_npglut_6hr_SNPs,
     file = "peaks_inside_npglut_6hr_SNPs_500bp.RData")

# load("multiomic_obj_with_new_peaks_labeled.RData")
# 
# multiomic_obj_new <-
#   UpdateSeuratObject(multiomic_obj_new)
# DefaultAssay(multiomic_obj_new) <- "ATAC"
# ATAC_fragment <-
#   Fragments(multiomic_obj_new)
# save(ATAC_fragment,
#      file = "multiomic_obj_new_470K_cells_frag_file.RData")
# load lexi's pre-made 500 bp peak data set ####
load("/nvmefs/scARC_Duan_018/018-029_combined_analysis/df_Lexi_new_peak_set_500bp.RData")

## use npglut_1hr only #####
peaks_500_width <-
  peaks_500_width[str_detect(string = names(peaks_500_width),
                             pattern = "npglut_1hr",
                             negate = F)]

## merge as one df
df_master_peak_set <-
  do.call(rbind,
          args = peaks_500_width)
sum(is.na(df_master_peak_set$chr))
unique(df_master_peak_set$chr)
df_master_peak_set <-
  df_master_peak_set[str_detect(string = df_master_peak_set$chr,
                                pattern = "^chr[1-9].*"), ]
length(sort(unique(df_master_peak_set$chr)))
## re-sort by chr names
df_master_peak_set$chr <-
  factor(df_master_peak_set$chr,
         levels = sort(unique(df_master_peak_set$chr)))

## split by chr name
all_peaks_by_chr <-
  split(x = df_master_peak_set,
        f = df_master_peak_set$chr)
names(all_peaks_by_chr)

discontinuous_peaks_list_by_chr <-
  vector(mode = "list",
         length = length(all_peaks_by_chr))
names(discontinuous_peaks_list_by_chr) <-
  names(all_peaks_by_chr)

## register a parallel backend for foreach #####
cluster_parallel <-
  makeCluster(8)
registerDoParallel(cluster_parallel)

## run the foreach for var i ####
discontinuous_peaks_list_by_chr <-
  foreach(i = 1:length(all_peaks_by_chr)) %dopar% {
    peaks_by_chr <-
      all_peaks_by_chr[[i]]
    peaks_by_chr <-
      peaks_by_chr[order(peaks_by_chr$start), ]
    for (j in 1:(nrow(peaks_by_chr))) {
      print(paste0("chr",
                   i,
                   ", peak ",
                   j))
      current_peak <-
        peaks_by_chr[j, ]
      if (j == 1) {
        previous_peak <- current_peak
        df_combined_peak <- current_peak
      } else {
        if (previous_peak$end <= current_peak$start) {
          df_combined_peak <-
            as.data.frame(rbind(df_combined_peak,
                                current_peak))
          previous_peak <- current_peak
        } else {
          if (previous_peak$score < current_peak$score) {
            df_combined_peak[nrow(df_combined_peak), ] <-
              unlist(current_peak)
            previous_peak <- current_peak
          }
        }
      }
    }
    return(df_combined_peak)
  }
stopCluster(cluster_parallel)
# 

names(discontinuous_peaks_list_by_chr) <-
  names(all_peaks_by_chr)

full_discontinuous_peaks_list <-
  do.call(rbind,
          args = discontinuous_peaks_list_by_chr)
sum(is.na(full_discontinuous_peaks_list$chr))

df_2_write <-
  full_discontinuous_peaks_list
df_2_write$chr <-
  str_split(string = rownames(df_2_write),
            pattern = "\\.",
            simplify = T)[, 1]
unique(df_2_write$chr)

df_2_write$start <-
  as.integer(df_2_write$start)
df_2_write$end <-
  as.integer(df_2_write$end)
df_2_write$strand <- '*'

GRanges_all_peaks <-
  makeGRangesFromDataFrame(df = df_2_write,
                           keep.extra.columns = T, 
                           na.rm = T)


## get npglut 1 hr specific SNPs
npglut_1hr_ASoC_SNPs <-
  vcf_lst_all_SNPs[['1hr_npglut']][vcf_lst_all_SNPs[['1hr_npglut']]$FDR < 0.05, ]
npglut_0hr_pBinom_lessthan_005_SNPs <-
  vcf_lst_all_SNPs[['0hr_npglut']][vcf_lst_all_SNPs[['0hr_npglut']]$pVal < 0.05, ]

npglut_1hr_specific_snps <-
  npglut_1hr_ASoC_SNPs[!(npglut_1hr_ASoC_SNPs$ID %in% npglut_0hr_pBinom_lessthan_005_SNPs$ID), ]
npglut_1hr_specific_snps <-
  npglut_1hr_specific_snps[!duplicated(npglut_1hr_specific_snps$ID), ]

npglut_1hr_specific_interval <-
  npglut_1hr_specific_snps[, c(1, 2, 2, 3, 15)]
npglut_1hr_specific_interval$POS.1 <-
  as.numeric(npglut_1hr_specific_interval$POS)
npglut_1hr_specific_interval$POS <-
  as.numeric(npglut_1hr_specific_interval$POS) - 1
npglut_1hr_specific_interval$strand <- "*"

colnames(npglut_1hr_specific_interval)[1:3] <-
  c("chr", "start", "end")
npglut_1hr_specific_interval <-
  makeGRangesFromDataFrame(df = npglut_1hr_specific_interval,
                           keep.extra.columns = T,
                           na.rm = T)

peaks_flanking_npglut_1hr_SNPs <-
  GRanges_all_peaks[GenomicRanges::nearest(x = npglut_1hr_specific_interval,
                                           subject = GRanges_all_peaks,
                                           ignore.strand = T,
                                           select = "arbitrary")]
peaks_flanking_npglut_1hr_SNPs_distance <-
  GenomicRanges::distanceToNearest(x = npglut_1hr_specific_interval,
                                   subject = peaks_flanking_npglut_1hr_SNPs,
                                   ignore.strand = T)
sum(peaks_flanking_npglut_1hr_SNPs_distance@elementMetadata$distance == 0)
hist(peaks_flanking_npglut_1hr_SNPs_distance@elementMetadata$distance)

peaks_flanking_npglut_1hr_SNPs <-
  peaks_flanking_npglut_1hr_SNPs[peaks_flanking_npglut_1hr_SNPs_distance@elementMetadata$distance == 0]
peaks_flanking_npglut_1hr_SNPs <-
  subsetByOverlaps(x = peaks_flanking_npglut_1hr_SNPs,
                   ranges = blacklist_hg38_unified,
                   invert = T)
peaks_flanking_npglut_1hr_SNPs <-
  peaks_flanking_npglut_1hr_SNPs[!duplicated(peaks_flanking_npglut_1hr_SNPs@elementMetadata@listData$peak_ID)]
save(peaks_flanking_npglut_1hr_SNPs,
     file = "peaks_inside_npglut_1hr_SNPs_500bp.RData")
