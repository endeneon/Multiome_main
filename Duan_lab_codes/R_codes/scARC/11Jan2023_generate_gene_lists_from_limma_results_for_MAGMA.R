# Chuxuan Li 1/11/2023
# generate list of genes from limma DE results for MAGMA analysis

# init ####
library(readr)
library(stringr)

pathlist <- list.files("./Analysis_part2_GRCh38/limma_DE_results/unfiltered_by_padj",
                       pattern = "[6|1]v0", full.names = T, include.dirs = F)
namelist <- str_extract(pathlist, "[A-Za-z]+_res_[1|6]v0")
dflist <- vector("list", length(pathlist))
for (i in 1:length(dflist)) {
  df <- read_csv(pathlist[i])
  df <- df[, c("genes", "logFC", "adj.P.Val")]
  df$set <- rep_len("set1", nrow(df))
  dflist[[i]] <- df
}

# filter by padj and logFC ####
for (i in 1:length(dflist)) {
  write.table(dflist[[i]][, c("set", "genes")][dflist[[i]]$logFC < 0 & dflist[[i]]$adj.P.Val < 0.05, ],
              file = paste0("./Analysis_part2_GRCh38/magma/limma_DEG_lists/", 
                            namelist[i], "_downregulated.tsv"), 
              sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(dflist[[i]][, c("set", "genes")][dflist[[i]]$logFC > 0 & dflist[[i]]$adj.P.Val < 0.05, ],
              file = paste0("./Analysis_part2_GRCh38/magma/limma_DEG_lists/", 
                            namelist[i], "_upregulated.tsv"), 
              sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(dflist[[i]][, c("set", "genes")][dflist[[i]]$adj.P.Val > 0.05, ],
              file = paste0("./Analysis_part2_GRCh38/magma/limma_DEG_lists/", 
                            namelist[i], "_nonsignificant.tsv"), 
              sep = "\t", quote = F, row.names = F, col.names = F)
}
