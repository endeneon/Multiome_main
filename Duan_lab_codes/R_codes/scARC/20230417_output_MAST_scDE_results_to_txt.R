# Chuxuan Li 04/17/2023
# output MAST case vs control DE analysis results into txt

library(readr)
library(stringr)
library(readxl)

rdatas <- list.files("./MAST_scDE/", pattern = "^results_", full.names = T, recursive = F)
rdatas <- rdatas[c(1:5, 7, 9:11)] # use the redo ones for counting
names <- str_extract(rdatas, "[A-Za-z]+_[0|1|6]hr")
for (i in 1:length(rdatas)) {
  load(rdatas[i])
  fcHurdle$coef <- (-1) * fcHurdle$coef
  write.table(fcHurdle, file = paste0("./MAST_scDE/DE_results_after_redo/unfiltered/", 
                                      names[i], "_unfiltered_full_results.csv"), 
              quote = F, sep = ",", row.names = F, col.names = T)
  write.table(fcHurdle[fcHurdle$fdr < 0.05 & fcHurdle$coef > 0, ],
              file = paste0("./MAST_scDE/DE_results_after_redo/upregulated_in_case/", 
                            names[i], "_upregulated_significant_results.csv"), 
              quote = F, sep = ",", row.names = F, col.names = T)
  write.table(fcHurdle[fcHurdle$fdr < 0.05 & fcHurdle$coef < 0, ],
              file = paste0("./MAST_scDE/DE_results_after_redo/downregulated_in_case/", 
                            names[i], "_downregulated_significant_results.csv"), 
              quote = F, sep = ",", row.names = F, col.names = T)
  write.table(fcHurdle$primerid[fcHurdle$fdr < 0.05 & fcHurdle$coef > 0],
              file = paste0("./MAST_scDE/DE_results_after_redo/gene_only/upregulated/", 
                            names[i], "_upregulated_significant_genes.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(fcHurdle$primerid[fcHurdle$fdr < 0.05 & fcHurdle$coef < 0],
              file = paste0("./MAST_scDE/DE_results_after_redo/gene_only/downregulated/", 
                            names[i], "_downregulated_significant_genes.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
}
