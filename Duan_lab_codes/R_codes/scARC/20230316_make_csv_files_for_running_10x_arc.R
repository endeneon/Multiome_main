# Chuxuan Li 03/16/2023
# make .csv files for each 10x multiomic library for running 10x arc pipeline

# init ####
library(stringr)
library(readr)

# make csv in loop ####
ATAC_samples <- unique(str_extract(sort(list.files("/nvmefs/scARC_Duan_022/ATAC_FASTQs/",
                                                   recursive = F, full.names = F)), 
                         "Duan_021_[0-9][0-9]-[0|1|6]_ATAC_FL"))
GEX_samples <- unique(str_extract(sort(list.files("/nvmefs/scARC_Duan_022/GEX_FASTQs/",
                                                   recursive = F, full.names = F)), 
                                  "Duan_021_[0-9][0-9]-[0|1|6]_GEX_FL"))
for (i in 1:length(ATAC_samples)) {
  df <-data.frame(fastqs = c("/home/cli/NVME/scARC_Duan_022/ATAC_FASTQs",
                             "/home/cli/NVME/scARC_Duan_022/GEX_FASTQs"),
                  sample = c(ATAC_samples[i], GEX_samples[i]),
                  library_type = c("Chromatin Accessibility",
                                   "Gene Expression"))
  write.csv(df, file = paste0("/nvmefs/scARC_Duan_022/csv_by_lib/Duan_022_", 
                              str_extract(ATAC_samples[i], "[0-9][0-9]-[0|1|6]"),
                              "_input_lib.csv"), quote = F, row.names = F, col.names = T)
}

