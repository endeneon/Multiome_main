# Chuxuan Li 05/05/2023
# output MAST DE analysis results into txt

library(readr)
library(stringr)
library(readxl)
setwd("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/oc2_mm10_analysis")
rdatas <- sort(list.files("./MAST_scDE/", pattern = "^results_", full.names = T, recursive = F))
names <- c("GABA_glut_combined", "GABA", "glut", "oligodendrocyte")

for (i in 1:length(rdatas)) {
  load(rdatas[i])
  write.table(fcHurdle, file = paste0("./MAST_scDE/DE_results/unfiltered/", 
                                      names[i], "_unfiltered_full_results.csv"), 
              quote = F, sep = ",", row.names = F, col.names = T)
  write.table(fcHurdle[fcHurdle$fdr < 0.05 & fcHurdle$coef > 0, ],
              file = paste0("./MAST_scDE/DE_results/upregulated_in_anesthesia/", 
                            names[i], "_upregulated_significant_results.csv"), 
              quote = F, sep = ",", row.names = F, col.names = T)
  write.table(fcHurdle[fcHurdle$fdr < 0.05 & fcHurdle$coef < 0, ],
              file = paste0("./MAST_scDE/DE_results/downregulated_in_anesthesia/", 
                            names[i], "_downregulated_significant_results.csv"), 
              quote = F, sep = ",", row.names = F, col.names = T)
  write.table(fcHurdle$primerid[fcHurdle$fdr < 0.05 & fcHurdle$coef > 0],
              file = paste0("./MAST_scDE/DE_results/gene_only/upregulated/", 
                            names[i], "_upregulated_significant_genes.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(fcHurdle$primerid[fcHurdle$fdr < 0.05 & fcHurdle$coef < 0],
              file = paste0("./MAST_scDE/DE_results/gene_only/downregulated/", 
                            names[i], "_downregulated_significant_genes.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
}
