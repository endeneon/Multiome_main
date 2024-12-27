# Chuxuan Li 05/17/2023
# get cell count on before and after QC by library


libs <- names(raw_obj_lst)
cell_count_df <- data.frame(library = c(libs, "total"),
                            raw = rep_len(0, length(libs) + 1),
                            after_demux = rep_len(0, length(libs) + 1),
                            after_QC = rep_len(0, length(libs) + 1))

load("raw_obj_lst.RData")
for (i in 1:length(libs)) {
  cell_count_df[i, 2] <- ncol(raw_obj_lst[[i]])
}
cell_count_df[nrow(cell_count_df), 2] <- sum(cell_count_df$raw)
rm(raw_obj_lst)

load("./cleaned_lst_no_nonhuman_cells_no_redone_lines.Rdata")
for (i in 1:length(cleaned_lst)) {
  cell_count_df[cell_count_df$library == names(cleaned_lst)[i], 3] <- ncol(cleaned_lst[[i]])
}
cell_count_df[nrow(cell_count_df), 3] <- sum(cell_count_df$after_demux)
rm(cleaned_lst)

load("./integrated_018-030_RNAseq_obj_test_QC.RData")
libs <- unique(integrated$lib.ident)
for (i in 1:length(libs)) {
  cell_count_df[cell_count_df$library == libs[i], 4] <- sum(integrated$lib.ident == libs[i])
}
cell_count_df[nrow(cell_count_df), 4] <- sum(cell_count_df$after_QC)

write.table(cell_count_df, file = "cell_counts_per_lib_raw_after_filtering_repeated_and_nonhuman_after_QC.csv",
            quote = F, sep = ",", row.names = F)
