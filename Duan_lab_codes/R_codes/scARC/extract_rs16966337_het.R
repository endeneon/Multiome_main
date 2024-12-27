# Siwei 24 Oct 2023

library(readxl)
library(stringr)

df_raw <-
  read_excel("rs16966337_het_excel_list.xlsx", 
             col_names = FALSE)

sample_names <-
  unique(unlist(df_raw[1, ]))

dir.create("ASoC_het_samples")

i <- 1

for (i in 1:length(sample_names)) {
  df_sub <-
    df_raw[ ,
            unlist(df_raw[1, ] %in% sample_names[i])]
  colnames(df_sub) <-
    c('...1',
      '...2')
  df_sub <-
    df_sub[11:nrow(df_sub), ]
  df_sub <-
    df_sub[!is.na(df_sub[[1]]), ]
  df_sub$sample_names <-
    str_split(string = df_sub$...1,
              pattern = "\\.var",
              simplify = T)[, 1]
  df_sub$sample_names <-
    str_replace(string = df_sub$sample_names,
                pattern = "-",
                replacement = "_")
  # df_sub$sample_names <-
  #   str_replace(string = df_sub$sample_names,
  #               pattern = "CW",
  #               replacement = "CW_")
  list_2_writeout <-
    data.frame(samples = df_sub$sample_names[str_detect(string = df_sub$...2,
                                                 pattern = "0\\/1")])
  writeout_file_name <-
    str_split(string = sample_names[i],
              pattern = "VQSR_",
              simplify = T)[, 2]
  writeout_file_name <-
    str_split(string = writeout_file_name,
              pattern = "_25Sept",
              simplify = T)[, 1]
  writeout_file_name <-
    str_c(writeout_file_name,
          '_het_sample_names.tsv')
  write.table(list_2_writeout,
              file = paste0("ASoC_het_samples/",
                            writeout_file_name),
              quote = F,
              sep = "\t",
              row.names = F,
              col.names = F)
}
