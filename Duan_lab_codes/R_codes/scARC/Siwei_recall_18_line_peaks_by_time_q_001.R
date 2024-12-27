# Siwei 06 Feb 2024
# re-call 18 line peaks at q=0.01
# to compare w/ PsychEncode peaks

# init####
library(Seurat)
library(Signac)
# library(sctransform)
# library(glmGamPoi)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
# library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

library(stringr)
library(future)

plan("multisession", workers = 6)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# load data
load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/ATAC_obj_added_MACS2_peaks.RData")

ATAC@assays$peaks <- NULL

DefaultAssay(ATAC) <- "ATAC"
unique(ATAC$cell.type)
ATAC <-
  ATAC[, !(ATAC$cell.type %in% "unknown")]

unique(ATAC$time.ident)
ATAC$time.ident <-
  factor(ATAC$time.ident,
         levels = c("0hr", "1hr", "6hr"))

peaks_by_time <-
  CallPeaks(object = ATAC,
            assay = "ATAC",
            group.by = "time.ident",
            macs2.path = "~/Data/Anaconda3-envs/aligners/bin/macs2",
            combine.peaks = F,
            cleanup = T,
            verbose = T,
            additional.args = '-q 0.01 --seed 42')
# save(peaks_by_time,
#      file = "peaks_by_time_3_cell_types_raw.RData")
save(peaks_by_time,
     file = "peaks_by_time_3_cell_types_w_lambda_raw.RData")
names(peaks_by_time) <-
  c("0hr", "1hr", "6hr")
# ATAC[['ATAC']] <-
#   as(ATAC[['ATAC']],
#      Class = "Assay5")
for (i in 1:length(peaks_by_time)) {
  peaks_by_time[[i]] <-
    keepStandardChromosomes(peaks_by_time[[i]], 
                            pruning.mode = "coarse")
  peaks_by_time[[i]] <-
    subsetByOverlaps(x = peaks_by_time[[i]],
                     ranges = blacklist_hg38_unified,
                     invert = T)
}

all_peaks_list <-
  as.data.frame(rbind(as.data.frame(peaks_by_time[[1]]),
                      as.data.frame(peaks_by_time[[2]]),
                      as.data.frame(peaks_by_time[[3]])))
all_peaks_list <-
  all_peaks_list[!(all_peaks_list$seqnames %in% c("chrX", "chrY", "chrM")), ]

all_peaks_list$start_500 <-
  all_peaks_list$start + all_peaks_list$relative_summit_position - 250
all_peaks_list$end_500 <-
  all_peaks_list$start + all_peaks_list$relative_summit_position + 250

all_peaks_by_chr <-
  split(x = all_peaks_list,
        f = all_peaks_list$seqnames)
all_peaks_by_chr <-
  all_peaks_by_chr[1:22]

all_peaks_by_time <-
  split(x = all_peaks_list,
        f = all_peaks_list$ident)

discontinuous_peaks_list_by_chr <-
  vector(mode = "list",
         length = length(all_peaks_by_chr))
names(discontinuous_peaks_list_by_chr) <-
  names(all_peaks_by_chr)

for (i in 1:length(all_peaks_by_chr)) {
  peaks_by_chr <-
    all_peaks_by_chr[[i]]
  peaks_by_chr <-
    peaks_by_chr[order(peaks_by_chr$start_500), ]
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
      if (previous_peak$end_500 <= current_peak$start_500) {
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
  discontinuous_peaks_list_by_chr[[i]] <-
    df_combined_peak
}

full_discontinuous_peaks_list <-
  do.call(rbind,
          args = discontinuous_peaks_list_by_chr)

full_discontinuous_peaks_by_time <-
  split(x = full_discontinuous_peaks_list,
        f = full_discontinuous_peaks_list$ident)

# write all peaks by time
for (i in 1:length(all_peaks_by_time)) {
  print(names(all_peaks_by_time)[i])
  df_2_write <-
    all_peaks_by_time[[i]]
  df_2_write <-
    df_2_write[, c(1, 13, 14, 6, 7, 5)]
  df_2_write$start_500 <-
    as.integer(df_2_write$start_500)
  df_2_write$end_500 <-
    as.integer(df_2_write$end_500)
  write.table(df_2_write,
              file = paste0("recalled_peaks_18_line_by_time_",
                            names(all_peaks_by_time)[i],
                            "_500bp.bed"),
              quote = F, sep = "\t",
              row.names = F, col.names = F)
}

# write discontinuous linear peaks
df_2_write <-
  full_discontinuous_peaks_list[, c(1, 13, 14, 6, 7, 5)]
df_2_write$chr <-
  str_split(string = rownames(df_2_write),
            pattern = "\\.",
            simplify = T)[, 1]
unique(df_2_write$chr)
df_2_write <-
  df_2_write[, c(7, 2:6)]
df_2_write$start_500 <-
  as.integer(df_2_write$start_500)
df_2_write$end_500 <-
  as.integer(df_2_write$end_500)
df_2_write$strand <- '*'

write.table(df_2_write,
            file = "recalled_peaks_18_line_by_time_non_overlapping_500bp.bed",
            quote = F, sep = "\t",
            row.names = F, col.names = F)
# make_501bp_peak_by_summit <-
#   function()