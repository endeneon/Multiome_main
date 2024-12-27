# Chuxuan Li 04/20/2023
# recount the number of significantly expressed genes from MAST scDE results

deg_counts <- array(dim = c(3, 10), 
                    dimnames = list(c("GABA", "nmglut", "npglut"),
                                    c("total genes", "passed filter - 0hr", 
                                      "upregulated - 0hr", "downregulated - 0hr",
                                      "passed filter - 1hr", 
                                      "upregulated - 1hr", "downregulated - 1hr",
                                      "passed filter - 6hr", 
                                      "upregulated - 6hr", "downregulated - 6hr")))
results_list <- list.files("./MAST_scDE", pattern = "^results", full.names = T)
results_list <- results_list[c(1:5, 7, 9:11)] # use the redo ones for counting
times_list <- str_extract(results_list, "[0|1|6]hr")
types_list <- str_sub(str_extract(results_list, "[A-Za-z]+_[0|1|6]"), end = -3L)
for (i in 1:length(results_list)) {
  load(results_list[i])
  fcHurdle$coef <- (-1) * fcHurdle$coef
  if (types_list[i] == "GABA") {
    rowind <- 1
  } else if (types_list[i] == "nmglut") {
    rowind <- 2
  } else {
    rowind <- 3
  }
  if (times_list[i] == "0hr") {
    deg_counts[rowind, 2] <- nrow(fcHurdle)
    deg_counts[rowind, 3] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef > 0)
    deg_counts[rowind, 4] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef < 0)
  } else if (times_list[i] == "1hr") {
    deg_counts[rowind, 5] <- nrow(fcHurdle)
    deg_counts[rowind, 6] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef > 0)
    deg_counts[rowind, 7] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef < 0)
  } else {
    deg_counts[rowind, 8] <- nrow(fcHurdle)
    deg_counts[rowind, 9] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef > 0)
    deg_counts[rowind, 10] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef < 0)
  }
}

deg_counts[, 1] <- rep_len(32823, 3)
write.table(deg_counts, file = "./MAST_scDE/deg_counts_after_redo.csv", quote = F, sep = ",")
