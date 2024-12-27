# Chuxuan Li 1/12/2022
# read and change .csv files into .bed files

# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)

library(readr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(graphics)
library(stringr)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)


set.seed(1105)

# four cell type files ####
setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/Output/four_cell_types_use_MAST/")
path_lst <- list.files(path = ".", pattern = "*.csv")
df_lst <- vector(mode = "list", length = length(path_lst))

for (i in 1:length(path_lst)){
  f <- read_csv(file = path_lst[[i]], col_names = F)
  #print("1")
  df_lst[[i]] <- data.frame(chr = str_split(string = f$X6, pattern = "-", simplify = T)[, 1],
                            start = gsub(pattern = "-",
                                         replacement = "",
                                         x = str_extract_all(string = f$X6,
                                                             pattern = "-[0-9]+-"
                                         )),
                            end = gsub(pattern = "-",
                                       replacement = "",
                                       x = str_extract_all(string = f$X6,
                                                           pattern = "-[0-9]+$"
                                       )),
                            id = paste0(f$X2, '^', f$X1, '^', f$X5), # log2FC^pval^pval_adj
                            p_val = f$X1,
                            strand = rep("+", length(f$X1)))
  #print("2")
  file_name <- paste0("./bed/", 
                      gsub(pattern = ".csv", replacement = ".bed", x = path_lst[[i]]))
  print(file_name)
  write.table(df_lst[[i]], 
              file = file_name, 
              quote = F, 
              sep = "\t", 
              row.names = F, 
              col.names = F)
}


# three cell type files ####
setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/Output/three_cell_types_da_peaks")
path_lst <- list.files(path = ".", pattern = "*.csv")
df_lst <- vector(mode = "list", length = length(path_lst))

for (i in 1:length(path_lst)){
  f <- read_csv(file = path_lst[[i]], col_names = F)
  df_lst[[i]] <- data.frame(chr = str_split(string = f$X6, pattern = "-", simplify = T)[, 1],
                            start = gsub(pattern = "-",
                                         replacement = "",
                                         x = str_extract_all(string = f$X6,
                                                             pattern = "-[0-9]+-"
                                         )),
                            end = gsub(pattern = "-",
                                       replacement = "",
                                       x = str_extract_all(string = f$X6,
                                                           pattern = "-[0-9]+$"
                                       )),
                            id = paste0(f$X2, '^', f$X1, '^', f$X5),
                            p_val = f$X1,
                            strand = rep("+", length(f$X1)))
  file_name <- paste0("./bed/", 
                      gsub(pattern = ".csv", replacement = ".bed", x = path_lst[[i]]))
  print(file_name)
  write.table(df_lst[[i]], 
              file = file_name, 
              quote = F, 
              sep = "\t", 
              row.names = F, 
              col.names = F)
}
