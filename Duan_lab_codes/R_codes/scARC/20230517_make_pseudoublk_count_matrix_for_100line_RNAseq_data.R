# Chuxuan Li 05/17/2023
# make pseudobulk count matrices for Lifan

load("/018-030_RNA_integrated_labeled_with_harmony.RData")
lines <- sort(unique(integrated_labeled$cell.line.ident))
times <- sort(unique(integrated_labeled$time.ident))
types <- sort(unique(integrated_labeled$cell.type))
for (i in 1:length(lines)) {
  for (j in 1:length(times)) {
    for (k in 1:length(types)) {
      df <-
        rowSums(integrated_labeled@assays$RNA@counts[, integrated_labeled$cell.line.ident == lines[i] &
                                                       integrated_labeled$time.ident == times[j] &
                                                       integrated_labeled$cell.type == types[k]])
      print(paste(lines[i], times[j], types[k]))
      if (i == 1 & j == 1 & k == 1) {
        pseudobulk_df <- df
        cnames <- paste(lines[i], times[j], types[k])
      } else {
        pseudobulk_df <- cbind(pseudobulk_df, df)
        cnames <- c(cnames, paste(lines[i], times[j], types[k]))
      }
    }
  }
}
colnames(pseudobulk_df) <- cnames
save(pseudobulk_df, file = "./018-030_100line_pseudobulk_df_by_celltype_time_cellline.RData")
