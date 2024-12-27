# Siwei 22 Feb 2022
# concatenate peaks and make unique, non-overlapping ones
# with parallel foreach package

# init
library(doParallel)
library(foreach)
library(readr)
library(stringr)

# check core numbers
detectCores() # 64

# register 10 cores fore doParallel
registerDoParallel(cores = 10)

# load data

# change here and load as many as you want to merge

# df_1 <- read_delim("scARC_5_lines_summit/500bp/GABA_bam_0hr_summits_cleaned.bed.500bp.bed", 
#                    delim = "\t", escape_double = FALSE, 
#                    col_names = FALSE, trim_ws = TRUE)
# df_2 <- read_delim("scARC_5_lines_summit/500bp/GABA_bam_1hr_summits_cleaned.bed.500bp.bed", 
#                    delim = "\t", escape_double = FALSE, 
#                    col_names = FALSE, trim_ws = TRUE)
# df_3 <- read_delim("scARC_5_lines_summit/500bp/GABA_bam_6hr_summits_cleaned.bed.500bp.bed", 
#                    delim = "\t", escape_double = FALSE, 
#                    col_names = FALSE, trim_ws = TRUE)

df_file_readin <-
  list.files(path = "scARC_5_lines_summit/500bp/",
             pattern = "\\.bed$",
             full.names = T)

df_list_by_file <- vector(mode = "list", length = length(df_file_readin))

i <- 1L

for (i in 1:length(df_file_readin)) {
  df_list_by_file[[i]] <- 
    read_delim(file = df_file_readin[i], 
               delim = "\t", escape_double = FALSE, 
               col_names = FALSE, trim_ws = TRUE)
  
}

# concatenate all peaks
df_master <-
  do.call("rbind", df_list_by_file)

# assign colnames to df_master
colnames(df_master) <- c("CHR", "START", "END", "NAME", "logQ", "STRAND")
# sort to make sure the concatenated input bed was not properly sorted
df_master <- df_master[order(df_master$CHR, df_master$START), ]

# split df_master by CHR
df_chr_list <- unique(df_master$CHR)
df_list <- vector(mode = "list",
                  length = length(df_chr_list))
names(df_list) <- df_chr_list

i <- 1L
for (i in 1:length(df_chr_list)) {
  df_list[[i]] <- df_master[df_master$CHR %in% df_chr_list[i], ]
}

##### define the ordering variable
sort_bed_per_chr <- function(df_raw) {

  # sort again in case the input bed was not properly sorted
  df_raw <- df_raw[order(df_raw$CHR, df_raw$START), ]
  i <- 1
  for (i in 1:nrow(df_raw)) {
    # print(i)
    if (i == 1) {
      df_output <- data.frame(df_raw[i, ],
                              stringsAsFactors = F)
      previous_peak <- df_raw[i, ]
    } else {
      if ((previous_peak[[3]][1] < df_raw[[2]][i]) ) { # bed file is semi-closed at the end
        df_temp <- df_raw[i, ]
        # print(colnames(df_temp))
        df_output <- data.frame(rbind(df_output,
                                      df_temp),
                                stringsAsFactors = F)
        previous_peak <- df_raw[i, ]
      }
    }
  }
  return(df_output)
}
#####

##### write the foreach() loop
k <- 1L

df_final_output <-
  foreach(k = 1:length(df_chr_list),
          .combine = 'rbind',
          .multicombine = T) %dopar%
  sort_bed_per_chr(df_list[[k]])

df_final_output <- 
  df_final_output[order(df_final_output$CHR, df_final_output$START), ]

peak_summarise_set <- 
  data.frame(str_split(string = df_final_output$NAME, 
                       pattern = "_",
                       simplify = T))
# factor(peak_summarise_set$X1)
# summary(peak_summarise_set$X1)
# levels(peak_summarise_set$X1)
sum(peak_summarise_set$X1 %in% "GABA") # 108013
sum(peak_summarise_set$X1 %in% "NEFMm") # 96116
sum(peak_summarise_set$X1 %in% "NEFMp") # 187782
sum(peak_summarise_set$X1 %in% "NPC") # 153295


write.table(df_final_output,
            file = "scARC_5_lines_summit/unique_non_overlap_peaks/sum_5_lines_GABA_Glut_unique_non_overlap_peaks.bed",
            quote = F, sep = "\t",
            row.names = F, col.names = F)

