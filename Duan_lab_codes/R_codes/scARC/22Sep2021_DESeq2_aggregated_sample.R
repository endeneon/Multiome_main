# Chuxuan 9/22/2021
# Load data from each group, separate into different cell lines, separate by the cell types,
# take the count matrices, aggregate the counts within one cell line for each gene, merge into a count matrix,
# write design matrix with time identity and cell line identity
# create deseq2 object

rm(list=ls())
# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(stringr)

# library(data.table)
library(ggplot2)

# library(Matrix)
# library(Rfast) # This one is incompatible with Seurat SCTransform(), load separately!
# library(plyr)
# library(dplyr)
# library(stringr)


library(future)
# set threads and parallelization

plan("multisession", workers = 8)
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

# load data
load("pre-SCT_QCed_g2_g8.RData")

# # load data use read10x_h5
# # note this h5 file contains both atac-seq and gex information
# g_2_0_raw <- 
#   Read10X_h5(filename = "//nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/libraries_2_0/outs/filtered_feature_bc_matrix.h5")
# # load gex information first
# g_2_0_raw_gex <- 
#   CreateSeuratObject(counts = g_2_0_raw$`Gene Expression`,
#                      project = "g_2_0_raw_gex")
# 
# g_2_1_raw <- 
#   Read10X_h5(filename = "//nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/libraries_2_1/outs/filtered_feature_bc_matrix.h5")
# # load gex information first
# g_2_1_raw_gex <- 
#   CreateSeuratObject(counts = g_2_1_raw$`Gene Expression`,
#                      project = "g_2_1_raw_gex")
# 
# g_2_6_raw <- 
#   Read10X_h5(filename = "//nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/libraries_2_6/outs/filtered_feature_bc_matrix.h5")
# # load gex information first
# g_2_6_raw_gex <- 
#   CreateSeuratObject(counts = g_2_6_raw$`Gene Expression`,
#                      project = "g_2_6_raw_gex")
# 
# 
# g_8_0_raw <- 
#   Read10X_h5(filename = "//nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/libraries_8_0/outs/filtered_feature_bc_matrix.h5")
# # load gex information first
# g_8_0_raw_gex <- 
#   CreateSeuratObject(counts = g_8_0_raw$`Gene Expression`,
#                      project = "g_8_0_raw_gex")
# 
# g_8_1_raw <- 
#   Read10X_h5(filename = "//nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/libraries_8_1/outs/filtered_feature_bc_matrix.h5")
# # load gex information first
# g_8_1_raw_gex <- 
#   CreateSeuratObject(counts = g_8_1_raw$`Gene Expression`,
#                      project = "g_8_1_raw_gex")
# 
# g_8_6_raw <- 
#   Read10X_h5(filename = "//nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/libraries_8_6/outs/filtered_feature_bc_matrix.h5")
# # load gex information first
# g_8_6_raw_gex <- 
#   CreateSeuratObject(counts = g_8_6_raw$`Gene Expression`,
#                      project = "g_8_6_raw_gex")
# 
# 
# # make a seurat object for subsequent operation
# # Note: all genes with total counts = 0 have been pre-removed
# g_2_0_gex <- g_2_0_raw_gex
# g_2_1_gex <- g_2_1_raw_gex
# g_2_6_gex <- g_2_6_raw_gex
# g_8_0_gex <- g_8_0_raw_gex
# g_8_1_gex <- g_8_1_raw_gex
# g_8_6_gex <- g_8_6_raw_gex
# 
# # create a list of Seurat objects for loop use
# obj_lst <- c(g_2_0_gex, g_2_1_gex, g_2_6_gex, g_8_0_gex, g_8_1_gex, g_8_6_gex)
# 
# store mitochondrial percentage in object meta data
# human starts with "MT-", rat starts with "Mt-"

g_2_0@meta.data$cell.line.ident <- "unknown"
g_2_1@meta.data$cell.line.ident <- "unknown"
g_2_6@meta.data$cell.line.ident <- "unknown"
g_8_0@meta.data$cell.line.ident <- "unknown"
g_8_1@meta.data$cell.line.ident <- "unknown"
g_8_6@meta.data$cell.line.ident <- "unknown"

## import cell ident barcodes
g20_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/backup/g_2_0_CD27_barcodes.txt", 
           col_names = FALSE)
g20_CD_27_barcodes <- unlist(g20_CD_27_barcodes)

g20_CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/backup/g_2_0_CD54_barcodes.txt", 
           col_names = FALSE)
g20_CD_54_barcodes <- unlist(g20_CD_54_barcodes)


g21_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_1_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
g21_CD_27_barcodes <- unlist(g21_CD_27_barcodes)

g21_CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_1_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
g21_CD_54_barcodes <- unlist(g21_CD_54_barcodes)


g26_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_6_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
g26_CD_27_barcodes <- unlist(g26_CD_27_barcodes)

g26_CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_6_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
g26_CD_54_barcodes <- unlist(g26_CD_54_barcodes)


g80_CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_0_common_CD08_CD25_CD26.best.CD_08", 
           col_names = FALSE)
g80_CD_08_barcodes <- unlist(g80_CD_08_barcodes)

g80_CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_0_common_CD08_CD25_CD26.best.CD_25", 
           col_names = FALSE)
g80_CD_25_barcodes <- unlist(g80_CD_25_barcodes)

g80_CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_0_common_CD08_CD25_CD26.best.CD_26", 
           col_names = FALSE)
g80_CD_26_barcodes <- unlist(g80_CD_26_barcodes)


g81_CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_1_common_CD08_CD25_CD26.best.CD_08", 
           col_names = FALSE)
g81_CD_08_barcodes <- unlist(g81_CD_08_barcodes)

g81_CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_1_common_CD08_CD25_CD26.best.CD_25", 
           col_names = FALSE)
g81_CD_25_barcodes <- unlist(g81_CD_25_barcodes)

g81_CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_1_common_CD08_CD25_CD26.best.CD_26", 
           col_names = FALSE)
g81_CD_26_barcodes <- unlist(g81_CD_26_barcodes)


g86_CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_6_common_CD08_CD25_CD26.best.CD_08", 
           col_names = FALSE)
g86_CD_08_barcodes <- unlist(g86_CD_08_barcodes)

g86_CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_6_common_CD08_CD25_CD26.best.CD_25", 
           col_names = FALSE)
g86_CD_25_barcodes <- unlist(g86_CD_25_barcodes)

g86_CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_6_common_CD08_CD25_CD26.best.CD_26", 
           col_names = FALSE)
g86_CD_26_barcodes <- unlist(g86_CD_26_barcodes)

## assign cell line identity

g_2_0@meta.data$cell.line.ident[g_2_0@assays$RNA@counts@Dimnames[[2]] %in% g20_CD_27_barcodes] <- "CD_27"
g_2_0@meta.data$cell.line.ident[g_2_0@assays$RNA@counts@Dimnames[[2]] %in% g20_CD_54_barcodes] <- "CD_54"
g_2_1@meta.data$cell.line.ident[g_2_1@assays$RNA@counts@Dimnames[[2]] %in% g21_CD_27_barcodes] <- "CD_27"
g_2_1@meta.data$cell.line.ident[g_2_1@assays$RNA@counts@Dimnames[[2]] %in% g21_CD_54_barcodes] <- "CD_54"
g_2_6@meta.data$cell.line.ident[g_2_6@assays$RNA@counts@Dimnames[[2]] %in% g26_CD_27_barcodes] <- "CD_27"
g_2_6@meta.data$cell.line.ident[g_2_6@assays$RNA@counts@Dimnames[[2]] %in% g26_CD_54_barcodes] <- "CD_54"

g_8_0@meta.data$cell.line.ident[g_8_0@assays$RNA@counts@Dimnames[[2]] %in% g80_CD_08_barcodes] <- "CD_08"
g_8_0@meta.data$cell.line.ident[g_8_0@assays$RNA@counts@Dimnames[[2]] %in% g80_CD_25_barcodes] <- "CD_25"
g_8_0@meta.data$cell.line.ident[g_8_0@assays$RNA@counts@Dimnames[[2]] %in% g80_CD_26_barcodes] <- "CD_26"
g_8_1@meta.data$cell.line.ident[g_8_1@assays$RNA@counts@Dimnames[[2]] %in% g81_CD_08_barcodes] <- "CD_08"
g_8_1@meta.data$cell.line.ident[g_8_1@assays$RNA@counts@Dimnames[[2]] %in% g81_CD_25_barcodes] <- "CD_25"
g_8_1@meta.data$cell.line.ident[g_8_1@assays$RNA@counts@Dimnames[[2]] %in% g81_CD_26_barcodes] <- "CD_26"
g_8_6@meta.data$cell.line.ident[g_8_6@assays$RNA@counts@Dimnames[[2]] %in% g86_CD_08_barcodes] <- "CD_08"
g_8_6@meta.data$cell.line.ident[g_8_6@assays$RNA@counts@Dimnames[[2]] %in% g86_CD_25_barcodes] <- "CD_25"
g_8_6@meta.data$cell.line.ident[g_8_6@assays$RNA@counts@Dimnames[[2]] %in% g86_CD_26_barcodes] <- "CD_26"

g_2_0@meta.data$cell.line.ident <-
  factor(g_2_0@meta.data$cell.line.ident)
g_2_1@meta.data$cell.line.ident <-
  factor(g_2_1@meta.data$cell.line.ident)
g_2_6@meta.data$cell.line.ident <-
  factor(g_2_6@meta.data$cell.line.ident)
g_8_0@meta.data$cell.line.ident <-
  factor(g_8_0@meta.data$cell.line.ident)
g_8_1@meta.data$cell.line.ident <-
  factor(g_8_1@meta.data$cell.line.ident)
g_8_6@meta.data$cell.line.ident <-
  factor(g_8_6@meta.data$cell.line.ident)

summary(g_8_6@meta.data$cell.line.ident)

# assign time identity
g_2_0@meta.data$time.ident <- "0hr"
g_8_0@meta.data$time.ident <- "0hr"
g_2_1@meta.data$time.ident <- "1hr"
g_8_1@meta.data$time.ident <- "1hr"
g_2_6@meta.data$time.ident <- "6hr"
g_8_6@meta.data$time.ident <- "6hr"

# assign group identity
g_2_0@meta.data$group.ident <- "group_2"
g_2_6@meta.data$group.ident <- "group_2"
g_2_1@meta.data$group.ident <- "group_2"
g_8_0@meta.data$group.ident <- "group_8"
g_8_1@meta.data$group.ident <- "group_8"
g_8_6@meta.data$group.ident <- "group_8"

glut_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/combined_barcodes/group_2_Glut_barcodes.txt",
                            header = F)
glut_barcodes <- unlist(glut_barcodes)
glut_barcodes <- str_sub(glut_barcodes, end = -3)

GABA_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/combined_barcodes/group_2_GABA_barcodes.txt",
                            header = F)
GABA_barcodes <- unlist(GABA_barcodes)
GABA_barcodes <- str_sub(GABA_barcodes, end = -3)

# separate cell lines from the object of each time point, combine them, then divide them by cell type
# eventually get 2*5 = 10 samples for pseudo bulk
CD_27 <- merge(subset(g_2_0, subset = cell.line.ident == "CD_27"), subset(g_2_1, subset = cell.line.ident == "CD_27"))

CD_27 <- merge(CD_27, subset(g_2_6, subset = cell.line.ident == "CD_27"))



## Siwei  subset CD_27 use cb_glut (since cb_glut is the glut subset of combined_gex)
# CD_27 <- subset(cb_glut, 
#                 subset = cell.line.ident == "CD_27")
# 
# sum((CD_27$time.ident %in% "0hr")) # 2017
# sum((CD_27$time.ident %in% "1hr")) # 2237
# sum((CD_27$time.ident %in% "6hr")) # 2066

### Cause of the issue:
### Barcodes merged from g_2_1 and g_2_6 are ended as 
### xxxxxxxxxxxxx-1_1, xxxxxxxxxxxxx-1_2, xxxxxxxxxxxxx-1_3
### e.g. execute the following line
unique(str_sub(CD_54@assays$RNA@counts@Dimnames[[2]],
               start = 17L))
### hence the barcodes could not match,
### need to trim the last 3 bytes of the Dimnames[[2]] as well
### use the following code solves the issue
###
CD_27 <- merge(subset(g_2_0, subset = cell.line.ident == "CD_27"), 
               c(subset(g_2_1, subset = cell.line.ident == "CD_27"),
                 subset(g_2_6, subset = cell.line.ident == "CD_27")))

CD_27$cell.type.ident <- "unknown"

CD_27$cell.type.ident[str_sub(CD_27@assays$RNA@counts@Dimnames[[2]],
                              end = -3) %in% 
                        str_sub(test_gb, 
                                end = -5)] <- "glut"
sum(str_sub(CD_27@assays$RNA@counts@Dimnames[[2]],
            end = -3) %in% 
      str_sub(test_gb, 
              end = -5))

sum((CD_27$cell.type.ident %in% "glut") & (CD_27$time.ident %in% "0hr")) # 2138
sum((CD_27$cell.type.ident %in% "glut") & (CD_27$time.ident %in% "1hr")) # 2381
sum((CD_27$cell.type.ident %in% "glut") & (CD_27$time.ident %in% "6hr")) # 2176
## Siwei ##


# use Siwei's solution to assign cell type to CD_27
CD_27$cell.type.ident <- "unknown"
CD_27$cell.type.ident[str_sub(CD_27@assays$RNA@counts@Dimnames[[2]], end = -3) %in% str_sub(glut_barcodes, end = -3)] <- "glut"
CD_27$cell.type.ident[str_sub(CD_27@assays$RNA@counts@Dimnames[[2]], end = -3) %in% str_sub(GABA_barcodes, end = -3)] <- "GABA"

CD_27_glut <- subset(CD_27, subset = cell.type.ident == "glut")
CD_27_glut_0hr <- subset(CD_27_glut, subset = time.ident == "0hr")
CD_27_glut_1hr <- subset(CD_27_glut, subset = time.ident == "1hr")
CD_27_glut_6hr <- subset(CD_27_glut, subset = time.ident == "6hr")

sum((CD_27$cell.type.ident %in% "glut") & (CD_27$time.ident %in% "0hr")) # 2138
sum((CD_27$cell.type.ident %in% "glut") & (CD_27$time.ident %in% "1hr")) # 2381
sum((CD_27$cell.type.ident %in% "glut") & (CD_27$time.ident %in% "6hr")) # 2176

CD_27_GABA <- subset(CD_27, subset = cell.type.ident == "GABA")
CD_27_GABA_0hr <- subset(CD_27_GABA, subset = time.ident == "0hr")
CD_27_GABA_1hr <- subset(CD_27_GABA, subset = time.ident == "1hr")
CD_27_GABA_6hr <- subset(CD_27_GABA, subset = time.ident == "6hr")

sum((CD_27$cell.type.ident %in% "GABA") & (CD_27$time.ident %in% "0hr")) # 1733
sum((CD_27$cell.type.ident %in% "GABA") & (CD_27$time.ident %in% "1hr")) # 1775
sum((CD_27$cell.type.ident %in% "GABA") & (CD_27$time.ident %in% "6hr")) # 1510

# assign cell type to CD_54, split into 0, 1, 6hrs
CD_54 <- merge(subset(g_2_0, subset = cell.line.ident == "CD_54"), 
               c(subset(g_2_1, subset = cell.line.ident == "CD_54"),
                 subset(g_2_6, subset = cell.line.ident == "CD_54")))

CD_54$cell.type.ident <- "unknown"
CD_54@meta.data$cell.type.ident[str_sub(CD_54@assays$RNA@counts@Dimnames[[2]], end = -3) %in% str_sub(glut_barcodes, end = -3)] <- "glut"
CD_54@meta.data$cell.type.ident[str_sub(CD_54@assays$RNA@counts@Dimnames[[2]], end = -3) %in% str_sub(GABA_barcodes, end = -3)] <- "GABA"

CD_54_glut <- subset(CD_54, subset = cell.type.ident == "glut")
CD_54_glut_0hr <- subset(CD_54_glut, subset = time.ident == "0hr")
CD_54_glut_1hr <- subset(CD_54_glut, subset = time.ident == "1hr")
CD_54_glut_6hr <- subset(CD_54_glut, subset = time.ident == "6hr")

CD_54_GABA <- subset(CD_54, subset = cell.type.ident == "GABA")
CD_54_GABA_0hr <- subset(CD_54_GABA, subset = time.ident == "0hr")
CD_54_GABA_1hr <- subset(CD_54_GABA, subset = time.ident == "1hr")
CD_54_GABA_6hr <- subset(CD_54_GABA, subset = time.ident == "6hr")

# CD_08
CD_08 <- merge(subset(g_8_0, subset = cell.line.ident == "CD_08"), 
               c(subset(g_8_1, subset = cell.line.ident == "CD_08"),
                 subset(g_8_6, subset = cell.line.ident == "CD_08")))

CD_08@meta.data$cell.type.ident[str_sub(CD_08@assays$RNA@counts@Dimnames[[2]], end = -3) %in% str_sub(glut_barcodes, end = -3)] <- "glut"
CD_08@meta.data$cell.type.ident[str_sub(CD_08@assays$RNA@counts@Dimnames[[2]], end = -3) %in% str_sub(GABA_barcodes, end = -3)] <- "GABA"
CD_08_glut <- subset(CD_08, subset = cell.type.ident == "glut")
CD_08_glut_0hr <- subset(CD_08_glut, subset = time.ident == "0hr")
CD_08_glut_1hr <- subset(CD_08_glut, subset = time.ident == "1hr")
CD_08_glut_6hr <- subset(CD_08_glut, subset = time.ident == "6hr")

CD_08_GABA <- subset(CD_08, subset = cell.type.ident == "GABA")
CD_08_GABA_0hr <- subset(CD_08_GABA, subset = time.ident == "0hr")
CD_08_GABA_1hr <- subset(CD_08_GABA, subset = time.ident == "1hr")
CD_08_GABA_6hr <- subset(CD_08_GABA, subset = time.ident == "6hr")

# CD_26
CD_26 <- merge(subset(g_8_0, subset = cell.line.ident == "CD_26"), 
               c(subset(g_8_1, subset = cell.line.ident == "CD_26"),
                 subset(g_8_6, subset = cell.line.ident == "CD_26")))
CD_26@meta.data$cell.type.ident[str_sub(CD_26@assays$RNA@counts@Dimnames[[2]], end = -3) %in% str_sub(glut_barcodes, end = -3)] <- "glut"
CD_26@meta.data$cell.type.ident[str_sub(CD_26@assays$RNA@counts@Dimnames[[2]], end = -3) %in% str_sub(GABA_barcodes, end = -3)] <- "GABA"

CD_26_glut <- subset(CD_26, subset = cell.type.ident == "glut")
CD_26_glut_0hr <- subset(CD_26_glut, subset = time.ident == "0hr")
CD_26_glut_1hr <- subset(CD_26_glut, subset = time.ident == "1hr")
CD_26_glut_6hr <- subset(CD_26_glut, subset = time.ident == "6hr")

CD_26_GABA <- subset(CD_26, subset = cell.type.ident == "GABA")
CD_26_GABA_0hr <- subset(CD_26_GABA, subset = time.ident == "0hr")
CD_26_GABA_1hr <- subset(CD_26_GABA, subset = time.ident == "1hr")
CD_26_GABA_6hr <- subset(CD_26_GABA, subset = time.ident == "6hr")

# CD_25
CD_25 <- merge(subset(g_8_0, subset = cell.line.ident == "CD_25"), 
               c(subset(g_8_1, subset = cell.line.ident == "CD_25"),
                 subset(g_8_6, subset = cell.line.ident == "CD_25")))

CD_25@meta.data$cell.type.ident[str_sub(CD_25@assays$RNA@counts@Dimnames[[2]], end = -3) %in% str_sub(glut_barcodes, end = -3)] <- "glut"
CD_25@meta.data$cell.type.ident[str_sub(CD_25@assays$RNA@counts@Dimnames[[2]], end = -3) %in% str_sub(GABA_barcodes, end = -3)] <- "GABA"

CD_25_glut <- subset(CD_25, subset = cell.type.ident == "glut")
CD_25_glut_0hr <- subset(CD_25_glut, subset = time.ident == "0hr")
CD_25_glut_1hr <- subset(CD_25_glut, subset = time.ident == "1hr")
CD_25_glut_6hr <- subset(CD_25_glut, subset = time.ident == "6hr")

CD_25_GABA <- subset(CD_25, subset = cell.type.ident == "GABA")
CD_25_GABA_0hr <- subset(CD_25_GABA, subset = time.ident == "0hr")
CD_25_GABA_1hr <- subset(CD_25_GABA, subset = time.ident == "1hr")
CD_25_GABA_6hr <- subset(CD_25_GABA, subset = time.ident == "6hr")

# finished subsetting the data into cell line/cell type/time objects

# merging glut/GABA count matrices, prepare for DESeq2
library(Rcpp)

library(DESeq2)
library(pasilla)

library(BiocParallel)
# register Biocparallel multithread
register(MulticoreParam(4))

# Rcpp function to convert dgCMatrix into normal matrix
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){

  int k = z.size() ;

  IntegerMatrix  mat(nrows, ncols);

  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }

  return mat;
}
' )


as_matrix <- function(mat){
  
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

# create a list of samples for future loop use, separate two cell types
glut_sample_lst <- c(CD_08_glut_0hr, CD_08_glut_1hr, CD_08_glut_6hr, 
                     CD_25_glut_0hr, CD_25_glut_1hr, CD_25_glut_6hr, 
                     CD_26_glut_0hr, CD_26_glut_1hr, CD_26_glut_6hr, 
                     CD_27_glut_0hr, CD_27_glut_1hr, CD_27_glut_6hr, 
                     CD_54_glut_0hr, CD_54_glut_1hr, CD_54_glut_6hr)
GABA_sample_lst <- c(CD_08_GABA_0hr, CD_08_GABA_1hr, CD_08_GABA_6hr, 
                     CD_25_GABA_0hr, CD_25_GABA_1hr, CD_25_GABA_6hr, 
                     CD_26_GABA_0hr, CD_26_GABA_1hr, CD_26_GABA_6hr, 
                     CD_27_GABA_0hr, CD_27_GABA_1hr, CD_27_GABA_6hr, 
                     CD_54_GABA_0hr, CD_54_GABA_1hr, CD_54_GABA_6hr)

# initialize list to store count matrices
glut_count_mat_lst <- vector(mode = "list", length = length(glut_sample_lst))
GABA_count_mat_lst <- vector(mode = "list", length = length(GABA_sample_lst))

# for loop to generate count matrices

for (i in 1:length(glut_sample_lst)){
  s <- glut_sample_lst[[i]]
  s_mat <- as_matrix(s@assays$RNA@counts)
  glut_count_mat_lst[[i]] <- s_mat
}

for (i in 1:length(GABA_sample_lst)){
  s <- GABA_sample_lst[[i]]
  s_mat <- as_matrix(s@assays$RNA@counts)
  GABA_count_mat_lst[[i]] <- s_mat
}

# initialize list to store row sum matrices
glut_sum_mat_lst <- vector(mode = "list", length = length(glut_sample_lst))
GABA_sum_mat_lst <- vector(mode = "list", length = length(GABA_sample_lst))

for (i in 1:length(glut_sample_lst)){
  sum_mat <- rowSums(glut_count_mat_lst[[i]])
  glut_sum_mat_lst[[i]] <- sum_mat
}
for (i in 1:length(GABA_sample_lst)){
  sum_mat <- rowSums(glut_count_mat_lst[[i]])
  GABA_sum_mat_lst[[i]] <- sum_mat
}

# merge all the samples into one matrix
glut_merged_count_matrix <- do.call(cbind, glut_sum_mat_lst)
GABA_merged_count_matrix <- do.call(cbind, GABA_sum_mat_lst)

# assign column names by sample names
names_glut <- c("CD_08_glut_0hr", "CD_08_glut_1hr", "CD_08_glut_6hr", 
                "CD_25_glut_0hr", "CD_25_glut_1hr", "CD_25_glut_6hr", 
                "CD_26_glut_0hr", "CD_26_glut_1hr", "CD_26_glut_6hr", 
                "CD_27_glut_0hr", "CD_27_glut_1hr", "CD_27_glut_6hr", 
                "CD_54_glut_0hr", "CD_54_glut_1hr", "CD_54_glut_6hr")
names_GABA <- c("CD_08_GABA_0hr", "CD_08_GABA_1hr", "CD_08_GABA_6hr", 
                "CD_25_GABA_0hr", "CD_25_GABA_1hr", "CD_25_GABA_6hr", 
                "CD_26_GABA_0hr", "CD_26_GABA_1hr", "CD_26_GABA_6hr", 
                "CD_27_GABA_0hr", "CD_27_GABA_1hr", "CD_27_GABA_6hr", 
                "CD_54_GABA_0hr", "CD_54_GABA_1hr", "CD_54_GABA_6hr")
colnames(glut_merged_count_matrix) <- names_glut
colnames(GABA_merged_count_matrix) <- names_GABA

# write manually the group and time columns, merge with samples names to get design matrices
group_col <- c(rep("group_8", 9), rep("group_2", 6))
time_col <- rep(c("0hr", "1hr", "6hr"), 5)
line_col <- str_sub(names_glut, 1, 5)
glut_design_matrix <- cbind(group_col, time_col, line_col)
GABA_design_matrix <- cbind(group_col, time_col, line_col)
# assign column names
rownames(glut_design_matrix) <- names_glut
rownames(GABA_design_matrix) <- names_GABA
# check we have the cells in the same order
all(rownames(GABA_design_matrix) == colnames(GABA_merged_count_matrix)) # TRUE
all(rownames(glut_design_matrix) == colnames(glut_merged_count_matrix)) # TRUE

# construct DESeqDataSet objects:
dds_glut <- DESeqDataSetFromMatrix(countData = glut_merged_count_matrix, 
                                   colData = glut_design_matrix,
                                   design = ~ time_col + line_col)

dds_GABA <- DESeqDataSetFromMatrix(countData = GABA_merged_count_matrix, 
                                   colData = GABA_design_matrix,
                                   design = ~ time_col + line_col)
# take a look
dds_glut
# class: DESeqDataSet 
# dim: 36601 15 
# metadata(1): version
# assays(1): counts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(0):
#   colnames(15): CD_08_glut_0hr CD_08_glut_1hr ... CD_54_glut_1hr CD_54_glut_6hr
# colData names(2): group_col time_col
dds_GABA
# class: DESeqDataSet 
# dim: 36601 15 
# metadata(1): version
# assays(1): counts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(0):
#   colnames(15): CD_08_GABA_0hr CD_08_GABA_1hr ... CD_54_GABA_1hr CD_54_GABA_6hr
# colData names(2): group_col time_col

# remove all-zero rows
dds_glut <- dds_glut[rowSums(counts(dds_glut)) > 0, ]
dds_GABA <- dds_GABA[rowSums(counts(dds_GABA)) > 0, ]

saveRDS(dds_glut, file = "dds_glut.rds")
saveRDS(dds_GABA, file = "dds_GABA.rds")


# differential expression analysis
dds_glut_dea <- DESeq(dds_glut, 
                      parallel = T)
res_glut <- results(dds_glut_dea)

dds_GABA_dea <- DESeq(dds_GABA, 
                      parallel = T)
res_GABA <- results(dds_GABA_dea)

# compare 0 and 1 hr, 0 and 6 hr
res_0_1_glut <- results(dds_glut_dea, 
                        contrast = c("time_col", "1hr", "0hr"), 
                        parallel = T)
res_0_6_glut <- results(dds_glut_dea, 
                        contrast = c("time_col", "6hr", "0hr"), 
                        parallel = T)
res_1_6_glut <- results(dds_glut_dea, 
                        contrast = c("time_col", "6hr", "1hr"), 
                        parallel = T)

res_0_1_GABA <- results(dds_GABA_dea, 
                        contrast = c("time_col", "1hr", "0hr"), 
                        parallel = T)
res_0_6_GABA <- results(dds_GABA_dea, 
                        contrast = c("time_col", "6hr", "0hr"), 
                        parallel = T)

# check specific genes

res_0_1_glut@listData$log2FoldChange[res_0_1_glut@rownames == "GAPDH"]
res_0_1_glut@listData$pvalue[res_0_1_glut@rownames == "GAPDH"]

rowSums(glut_merged_count_matrix)[rownames(glut_merged_count_matrix) == "BDNF"]
rowSums(glut_merged_count_matrix)[rownames(glut_merged_count_matrix) == "FOS"]

res_0_6_glut@listData$log2FoldChange[res_0_6_glut@rownames == "BDNF"]
res_0_6_glut@listData$pvalue[res_0_6_glut@rownames == "BDNF"]

res_0_1_glut@listData$log2FoldChange[res_0_1_glut@rownames == "BDNF"]
res_0_1_glut@listData$pvalue[res_0_1_glut@rownames == "BDNF"]

# PNOC
res_0_1_GABA@listData$log2FoldChange[res_0_1_GABA@rownames == "PNOC"]
res_0_1_GABA@listData$pvalue[res_0_1_GABA@rownames == "PNOC"]

res_0_6_GABA@listData$log2FoldChange[res_0_6_GABA@rownames == "PNOC"]
res_0_6_GABA@listData$pvalue[res_0_6_GABA@rownames == "PNOC"]
