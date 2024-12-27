# Siwei 05 Aug 2022
# calculate histogram of BAM file sizes to determine resampling 

# init
library(ggplot2)
library(readr)

library(parallel)
library(stringr)
library(Rsamtools)

library(RColorBrewer)
library(grDevices)

## get all bam files (un de-deduped, will deal with deduped later)

setwd("~/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/bam_dump_4_fasta_05Jul2022_hg38/subsetted_raw_bams/")

# get all bam files
input_file_list <-
  list.files(path = ".", 
             pattern = "0hr_CD_63_GABA_bwa_barcode_WASPed.bam$",
             full.names = T, 
             recursive = T)

countBam(input_file_list[1], 
         index = paste(input_file_list[1], 
                       '.bai', 
                       sep = ""))

## remove file names contain bwa, etc.
input_file_list <-
  input_file_list[str_detect(string = input_file_list,
                             pattern = "_bwa_",
                             negate = T)]

input_file_list <-
  input_file_list[str_detect(string = input_file_list,
                             pattern = "_GABA",
                             negate = F)]

df_file_info <-
  file.info(input_file_list,
            extra_cols = T)

hist(df_file_info$size, breaks = 50)
median(df_file_info$size) # NEFM_pos: 2.5e9; NEFM_neg: 1.8e9; GABA: 1.8e9
mean(df_file_info$size) # 2290912524
sum(df_file_info$size < 1e9)

sum(df_file_info$size) # 373418741464
sum(df_file_info$size[df_file_info$size > 2290912520]) # 2.39033e+11
min(df_file_info$size) * nrow(df_file_info) # 72148771663


df_file_info <-
  df_file_info[order(df_file_info$size), ]
df_file_info$order <- 1:nrow(df_file_info)


ggplot(df_file_info,
       aes(x = order,
           y = size)) +
  geom_col() +
  theme_classic()
