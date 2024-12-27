# Chuxuan 9/20/2021
# export count matrix and design matrix from Seurat object,
# use DESeq2 to do pseudo bulk analysis

#rm(list=ls())
# init
library(Seurat)
library(sctransform)
library(glmGamPoi)

library(future)

library(readr)
library(ggplot2)

library(Rcpp)

library(DESeq2)
library(pasilla)

library(BiocParallel)

# set threads and parallelization
plan("multisession", workers = 6)
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

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

# read in labeled Seurat object from previous Seurat analysis
combined_mat <- readRDS("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Seurat_objects/combined_g2_g8.rds")
glut_mat <- subset(combined_mat, subset = cell.type.ident == "Glut")
GABA_mat <- subset(combined_mat, subset = cell.type.ident == "GABA")


# use the RNA (preprocessed) count matrix, use Rcpp function to convert to large matrix
glut_count_matrix <- as_matrix(glut_mat@assays$RNA@counts)
glut_count_matrix <- as.matrix(glut_count_matrix)
GABA_count_matrix <- as_matrix(GABA_mat@assays$RNA@counts)
GABA_count_matrix <- as.matrix(GABA_count_matrix)

save(list = c("glut_count_matrix", "glut_design_matrix"), 
     file = "glut_matrices.RData", 
     compression_level = 9, 
     precheck = T)
save(list = c("GABA_count_matrix", "GABA_design_matrix"), 
     file = "GABA_matrices.RData", 
     compression_level = 9, 
     precheck = T)


# read in time, cell type, group from Seurat object
glut_design_matrix <- cbind(glut_mat@assays$RNA@counts@Dimnames[[2]], 
                            glut_mat@meta.data$time.ident, 
                            glut_mat@meta.data$group.ident)
GABA_design_matrix <- cbind(GABA_mat@assays$RNA@counts@Dimnames[[2]], 
                            GABA_mat@meta.data$time.ident, 
                            GABA_mat@meta.data$group.ident)


# assign row names of the design matrix to be barcodes
rownames(glut_design_matrix) <- glut_design_matrix[, 1]
glut_design_matrix <- glut_design_matrix[, 2:3]
colnames(glut_design_matrix) <- c("time.ident", "group.ident")

rownames(GABA_design_matrix) <- GABA_design_matrix[, 1]
GABA_design_matrix <- GABA_design_matrix[, 2:3]
colnames(GABA_design_matrix) <- c("time.ident", "group.ident")


# check we have the cells in the same order
all(rownames(GABA_design_matrix) == colnames(GABA_count_matrix)) # TRUE
all(rownames(glut_design_matrix) == colnames(glut_count_matrix)) # TRUE

# construct a DESeqDataSet object:
dds_glut <- DESeqDataSetFromMatrix(countData = glut_count_matrix, 
                                   colData = glut_design_matrix,
                                   design = ~ time.ident + group.ident)

dds_GABA <- DESeqDataSetFromMatrix(countData = GABA_count_matrix, 
                                   colData = GABA_design_matrix,
                                   design = ~ time.ident + group.ident)
dds_glut
# class: DESeqDataSet 
# dim: 36601 37545 
# metadata(1): version
# assays(1): counts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(0):
#   colnames(37545): AAACAGCCAGTAAAGC-1_1_1 AAACAGCCATCCAGGT-1_1_1 ... TTTGTTGGTTCCTCCT-1_2_3 TTTGTTGGTTGGCCGA-1_2_3
# colData names(2): time.ident group.ident
dds_GABA
# class: DESeqDataSet 
# dim: 36601 16936 
# metadata(1): version
# assays(1): counts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(0):
#   colnames(16936): AAACAGCCAAGGTATA-1_1_1 AAACAGCCAGCATGAG-1_1_1 ... TTTGTGGCAAGACTCC-1_2_3 TTTGTGTTCGTTAAGC-1_2_3
# colData names(2): time.ident group.ident

# remove all-zero rows
dds_glut <- dds_glut[rowSums(counts(dds_glut)) > 0, ]
dds_GABA <- dds_GABA[rowSums(counts(dds_GABA)) > 0, ]


# sum(is.na(counts(dds_glut)))
# sum(is.infinite(counts(dds_glut)))
# test <- counts(dds_glut)
# sum(is.integer(dds_glut@assays@data$counts))

# save into a different object
dds_glut_nonzero <- dds_glut
# add one to remove zero values
dds_glut_nonzero@assays@data$counts <- dds_glut_nonzero@assays@data$counts + 1L

dds_GABA_nonzero <- dds_GABA
dds_GABA_nonzero@assays@data$counts <- dds_GABA_nonzero@assays@data$counts + 1L

save.image("DESeq2_pre_analysis.RData")

# differential expression analysis
# !slow!
dds_glut_nonzero <- DESeq(dds_glut_nonzero, 
                          parallel = T)
res_glut <- results(dds_glut_nonzero)

dds_GABA_nonzero <- DESeq(dds_GABA_nonzero, 
                          parallel = T)
res_GABA <- results(dds_GABA_nonzero)

# check histogram of counts in each cell
hist(colSums(counts(dds_glut)),
     main = "glut count/gene distribution",
     xlim = c(0, 6e4),
     col = "darkred",
     breaks = 1000
)

hist(colSums(counts(dds_GABA)),
     main = "GABA cell count/gene distribution",
     xlim = c(0, 6e4),
     col = "darkred",
     breaks = 100
)

# check histogram of counts for each gene
hist(rowSums(counts(dds_glut)),
     main = "glut count/cell distribution",
     xlim = c(0, 100),
     col = "darkgreen",
     breaks = 1000000
)
sum(rowSums(counts(dds_glut)) == 3)
hist(rowSums(counts(dds_GABA)),
     main = "GABA count/cell distribution",
     xlim = c(0, 2e4),
     col = "darkgreen",
     breaks = 10000
)
# check the total counts of low-expression genes
# BDNF for glut
# PNOC for GABA
sum(dds_glut@assays@data@listData$counts[res_glut@rownames == "BDNF"])
sum(dds_glut_nonzero@assays@data@listData$counts[res_glut@rownames == "BDNF"])

sum(dds_GABA@assays@data@listData$counts[res_glut@rownames == "PNOC"])
sum(dds_GABA_nonzero@assays@data@listData$counts[res_glut@rownames == "PNOC"])


# compare 0hr and 1hr, 0hr and 6hr
res_0_1_glut <- results(dds_glut_nonzero, 
                   contrast = c("time.ident", "1hr", "0hr"), 
                   parallel = T)
res_0_6_glut <- results(dds_glut_nonzero, 
                   contrast = c("time.ident", "6hr", "0hr"), 
                   parallel = T)

res_0_1_GABA <- results(dds_GABA_nonzero, 
                   contrast = c("time.ident", "1hr", "0hr"), 
                   parallel = T)
res_0_6_GABA <- results(dds_GABA_nonzero, 
                   contrast = c("time.ident", "6hr", "0hr"), 
                   parallel = T)
# res_1_6 <- results(dds_nonzero, contrast = c("time.ident", "1hr", "6hr"))

res_0_1_glut@listData$log2FoldChange[res_0_1_glut@rownames == "BDNF"]
res_0_1_glut@listData$pvalue[res_0_1_glut@rownames == "BDNF"]

res_0_6_glut@listData$log2FoldChange[res_0_6_glut@rownames == "BDNF"]
res_0_6_glut@listData$pvalue[res_0_6_glut@rownames == "BDNF"]
