# Siwei 19 Oct 2022
# calc 15 Aug 2022 results of DE ASoC SNP numbers

library(stringr)
library(readr)


cell_type <- "NEFM_pos"

SNPFile_list <-
  list.files(path = "NotDuplicateReadFilter_15Aug2022/",
             pattern = paste(".*",
                             cell_type,
                             ".*\\.tsv",
                             sep = ""),
             full.names = T)




hr_0 <- 
  read_delim(SNPFile_list[1], 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
hr_1 <-
  read_delim(SNPFile_list[2], 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
hr_6 <-
  read_delim(SNPFile_list[3], 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

hr_0 <- 
  hr_0[hr_0$FDR < 0.05, ]
hr_1 <- 
  hr_1[hr_1$FDR < 0.05, ]
hr_6 <- 
  hr_6[hr_6$FDR < 0.05, ]


union_list <-
  union(union(hr_0$ID, hr_1$ID), hr_6$ID)
intersect_list <-
  
