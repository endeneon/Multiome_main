# Siwei 05 Feb 2024
# Test if there is any log2 diff changes caused by sex factor

# Need to use sva/Combat to get normalised counts
# Adapt to 100 lines

# init ####
library(Seurat)
library(Signac)

library(limma)
library(edgeR)
library(ggplot2)
library(stringr)

library(sva)
# library(harmony)

library(reshape2)

library(readr)
library(readxl)

plan("multisession", workers = 4)
set.seed(42)
options(future.globals.maxSize = 229496729600)
options(Seurat.object.assay.version = "v5")

# Aux functions ####
createDGE <- 
  function(count_matrix) {
    y <- 
      DGEList(counts = count_matrix, 
              genes = rownames(count_matrix), 
              group = cell_lines)
    A <- rowSums(y$counts)
    hasant <- rowSums(is.na(y$genes)) == 0
    y <- y[hasant, , keep.lib.size = F]
    print(dim(y))
    return(y)
  }

filterByCpm <- 
  function(df1, df2, cutoff, proportion = 1) {
    cpm1 <- cpm(df1)
    cpm2 <- cpm(df2)
    passfilter <- (rowSums(cpm1 >= cutoff) >= ncol(cpm1) * proportion |
                     # rowSums(cpm2 >= cutoff) >= ncol(cpm2) * proportion |
                     (rowSums(cpm1 > 0) & (rowSums(cpm2 >= cutoff) >= ncol(cpm2) * proportion))) # number of samples from either group > ns
    return(passfilter)
  }

cnfV <- 
  function(y, design) {
    y <- calcNormFactors(y)
    v <- voom(y, design, plot = F)
    return(v)  
  }
# deal three cell types one by one, 
# extract them from adj_mat_lst

# load data #####
load("~/NVME/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/pseudo_bulk_100_lines_data_typexlinextime_mat_lst_05Feb2024.RData")
load("~/NVME/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/multiomic_obj_4_plotting_gene_expression_30Jan2024.RData")

## assemble adj_mat_list
# adj_mat_list <-
#   vector(mode = "list",
#          length = 3L)
# names(adj_mat_list) <-
#   c("GABA",
#     "nmglut",
#     "npglut")

npglut_030 <-
  adj_mat_lst[[3]]
colnames(npglut_030)
rownames(npglut_030)
save(npglut_030,
     file = "npglut_030_matrix.RData")
