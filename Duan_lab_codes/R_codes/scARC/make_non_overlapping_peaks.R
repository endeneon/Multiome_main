# Siwei 21 Feb 2022
# make unique non-overlapping peaks from bed file
# peak size can be customised, but usually +/- 250 bp from summit
# (Greenberg, 2016)

# init
library(readr)

half_peak_size = 250

# this part will "flatten out" any overlapping peaks based on their
# position and value (column 5). If two peaks overlap, use the one
# with larger -logQ value (column 5)

# check the peak no. that passes different -logQ threshold later


bed_list <- list.files(path = "scARC_5_lines_summit/cleaned/", 
                       pattern = "\\.bed$",
                       full.names = F)

k <- 1
for (k in 1:length(bed_list)) {
  writeout_df <- make_500_bp_bed(bed_list[k])
  write.table(writeout_df,
              file = paste('scARC_5_lines_summit/500bp/', bed_list[k], '.500bp.bed', sep = ""),
              quote = F, sep = "\t", 
              row.names = F, col.names = F)
}



make_500_bp_bed <- function(x) {
  df_raw <- read_delim(paste("scARC_5_lines_summit/cleaned/", x, sep = ""), 
                       delim = "\t", escape_double = FALSE, 
                       col_names = FALSE, trim_ws = TRUE)
  colnames(df_raw) <- c("CHR", "START", "END", "NAME", "logQ")
  df_raw$STRAND <- "*"
  
  # sort again in case the input bed was not properly sorted
  df_raw <- df_raw[order(df_raw$CHR, df_raw$START), ]
  df_raw$START <- df_raw$START - 250
  df_raw$END <- df_raw$END + 250
  
  i <- 1
  
  for (i in 1:nrow(df_raw)) {
    # print(i)
    if (i == 1) {
      df_output <- data.frame(df_raw[i, ],
                              stringsAsFactors = F)
      previous_peak <- df_raw[i, ]
    } else {
      if ((previous_peak[[3]][1] < df_raw[[2]][i]) | 
          (previous_peak[[1]][1] != df_raw[[1]][i])) {
        df_temp <- df_raw[i, ]
        # print(colnames(df_temp))
        df_output <- data.frame(rbind(df_output,
                                      df_temp),
                                stringsAsFactors = F)
        if ((previous_peak[[1]][1] != df_raw[[1]][i])) {
          print(df_raw[[1]][i])
        }
        previous_peak <- df_raw[i, ]
        
      }
    }
  }
  
  return(df_output)
}

