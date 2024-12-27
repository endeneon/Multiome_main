# Chuxuan 9/14/2021
# Already obtained barcodes for cells to exclude, now subset human only data with the barcodes, 
#merge cell lines, and preprocess the data; finally, check marker gene expression

rm(list=ls())
# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)
library(future)
# set threads and parallelization
plan("multisession", workers = 6)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)


# load human only data, create Seurat objects
# load data use read10x_h5
# note this h5 file contains both atac-seq and gex information
g_2_0_read <- 
  Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_2_0/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_2_0_read_gex <- 
  CreateSeuratObject(counts = g_2_0_read$`Gene Expression`,
                     project = "g_2_0_read_gex")

g_2_1_read <- 
  Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_2_1/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_2_1_read_gex <- 
  CreateSeuratObject(counts = g_2_1_read$`Gene Expression`,
                     project = "g_2_1_read_gex")

g_2_6_read <- 
  Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_2_6/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_2_6_read_gex <- 
  CreateSeuratObject(counts = g_2_6_read$`Gene Expression`,
                     project = "g_2_6_read_gex")


g_8_0_read <- 
  Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_8_0/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_8_0_read_gex <- 
  CreateSeuratObject(counts = g_8_0_read$`Gene Expression`,
                     project = "g_8_0_read_gex")

g_8_1_read <- 
  Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_8_1/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_8_1_read_gex <- 
  CreateSeuratObject(counts = g_8_1_read$`Gene Expression`,
                     project = "g_8_1_read_gex")

g_8_6_read <- 
  Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_8_6/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_8_6_read_gex <- 
  CreateSeuratObject(counts = g_8_6_read$`Gene Expression`,
                     project = "g_8_6_read_gex")

# make a separate seurat object for subsequent operation
# Note: all genes with total counts = 0 have been pre-removed
g_2_0 <- g_2_0_read_gex
g_2_1 <- g_2_1_read_gex
g_2_6 <- g_2_6_read_gex
g_8_0 <- g_8_0_read_gex
g_8_1 <- g_8_1_read_gex
g_8_6 <- g_8_6_read_gex

# read barcode files
g_2_0_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_2_0_exclude_barcodes.txt",
                             header = F)
g_2_1_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_2_1_exclude_barcodes.txt",
                             header = F)
g_2_6_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_2_6_exclude_barcodes.txt",
                             header = F)
g_8_0_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_8_0_exclude_barcodes.txt",
                             header = F)
g_8_1_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_8_1_exclude_barcodes.txt",
                             header = F)
g_8_6_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_8_6_exclude_barcodes.txt",
                             header = F)

# test6 <- c(unlist(g_2_6_barcodes, use.names = F), unlist(g_8_6_barcodes, use.names = F))

# select the data with barcodes not in the barcode files
# the subset() function does not take counts@dimnames, has to give cells an ident first
g_2_0@meta.data$exclude.status <- "include"
g_2_0@meta.data$exclude.status[colnames(g_2_0@assays$RNA@counts) %in% g_2_0_barcodes] <- "exclude"
g_2_0 <- subset(x = g_2_0, subset = exclude.status %in% "include")

g_2_1@meta.data$exclude.status <- "include"
g_2_1@meta.data$exclude.status[colnames(g_2_1@assays$RNA@counts) %in% g_2_1_barcodes] <- "exclude"
g_2_1 <- subset(x = g_2_1, subset = exclude.status %in% "include")

g_2_6@meta.data$exclude.status <- "include"
g_2_6@meta.data$exclude.status[colnames(g_2_6@assays$RNA@counts) %in% g_2_0_barcodes] <- "exclude"
g_2_6 <- subset(x = g_2_6, subset = exclude.status %in% "include")

g_8_0@meta.data$exclude.status <- "include"
g_8_0@meta.data$exclude.status[colnames(g_8_0@assays$RNA@counts) %in% g_8_0_barcodes] <- "exclude"
g_8_0 <- subset(x = g_8_0, subset = exclude.status %in% "include")

g_8_1@meta.data$exclude.status <- "include"
g_8_1@meta.data$exclude.status[colnames(g_8_1@assays$RNA@counts) %in% g_8_1_barcodes] <- "exclude"
g_8_1 <- subset(x = g_8_1, subset = exclude.status %in% "include")

g_8_6@meta.data$exclude.status <- "include"
g_8_6@meta.data$exclude.status[colnames(g_8_6@assays$RNA@counts) %in% g_8_0_barcodes] <- "exclude"
g_8_6 <- subset(x = g_8_6, subset = exclude.status %in% "include")

#save.image("pre-SCT_QCed_g2_g8.RData")

# merge the cell lines at 0, 1, 6
time_0 = merge(g_2_0, g_8_0)
time_1 = merge(g_2_1, g_8_1)
time_6 = merge(g_2_6, g_8_6)

# assign time point ident
time_0@meta.data$time.ident <- "0hr"
time_1@meta.data$time.ident <- "1hr"
time_6@meta.data$time.ident <- "6hr"

# merge 0, 1, 6 time points into one object
combined_gex_read <- merge(time_0, c(time_1, time_6))

# analyze gene expression of marker genes
# pre-treat data 

# store mitochondrial percentage in object meta data
# human starts with "MT-", rat starts with "Mt-"
combined_gex <- 
  PercentageFeatureSet(combined_gex_read,
                       pattern = c("^MT-"),
                       col.name = "percent.mt")

#use SCTransform()
combined_gex <-
  SCTransform(combined_gex,
              vars.to.regress = "percent.mt",
              method = "glmGamPoi",
              variable.features.n = 10000, #do not go above, will crash if you do
              verbose = T,
              seed.use = 42)

# combined_gex <- NormalizeData(combined_gex, 
#                               normalization.method = "LogNormalize", 
#                               scale.factor = 10000)
# combined_gex <- FindVariableFeatures(combined_gex, 
#                                      selection.method = "vst", 
#                                      nfeatures = 2000)
# all.genes <- rownames(combined_gex)
# combined_gex <- ScaleData(combined_gex, 
#                           features = all.genes)

combined_gex <- RunPCA(combined_gex, 
                    verbose = T)
combined_gex <- RunUMAP(combined_gex,
                     dims = 1:30,
                     verbose = T)
combined_gex <- FindNeighbors(combined_gex,
                           dims = 1:30,
                           verbose = T)
combined_gex <- FindClusters(combined_gex,
                          verbose = T, 
                          resolution = 0.7)

saveRDS(combined_gex, file = "combined_seurat_object.rds")

# make dimension plot
DimPlot(combined_gex, 
        label = T) +
  NoLegend()

library(stringr)
combined_gex@meta.data$cell.line.ident[str_sub(combined_gex@assays$SCT@counts@Dimnames[[2]], end = -5L) %in% 
                                         c(g_2_0_CD_27_barcodes, g_2_1_CD_27_barcodes, g_2_6_CD_27_barcodes)] <- "CD_27"
combined_gex@meta.data$cell.line.ident[str_sub(combined_gex@assays$SCT@counts@Dimnames[[2]], end = -5L) %in% 
                                         c(g_8_0_CD_25_barcodes, g_8_1_CD_25_barcodes, g_8_6_CD_25_barcodes)] <- "CD_25"
combined_gex@meta.data$cell.line.ident[str_sub(combined_gex@assays$SCT@counts@Dimnames[[2]], end = -5L) %in% 
                                         c(g_8_0_CD_26_barcodes, g_8_1_CD_26_barcodes, g_8_6_CD_26_barcodes)] <- "CD_26"
combined_gex@meta.data$cell.line.ident[str_sub(combined_gex@assays$SCT@counts@Dimnames[[2]], end = -5L) %in% 
                                         c(g_8_0_CD_08_barcodes, g_8_1_CD_08_barcodes, g_8_6_CD_08_barcodes)] <- "CD_08"
combined_gex@meta.data$cell.line.ident[str_sub(combined_gex@assays$SCT@counts@Dimnames[[2]], end = -5L) %in% 
                                         c(g_2_0_CD_54_barcodes, g_2_1_CD_54_barcodes, g_2_6_CD_54_barcodes)] <- "CD_54"

DimPlot(combined_gex,
        cols = c("red", "green", "yellow", "cyan", "pink", "black"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "cell.line.ident",
        label = F) +
  ggtitle("clustering by cell line")

DimPlot(combined_gex,
        cols = c("white", "white", "white", "cyan", "white", "white"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "cell.line.ident",
        label = F) +
  ggtitle("clustering by cell line")

DimPlot(combined_gex,
        cols = c("red", "green", "black"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("clustering by time")

# check specific genes
FeaturePlot(combined_gex, 
            features = c("GAD1", "SLC17A6", "SLC17A7", 
                         "SOX2"), 
            pt.size = 0.2,
            ncol = 2)

StackedVlnPlot(obj = combined_gex, 
               features = c("GAD1", "GAD2", "PNOC", 
                            "SLC17A6", "SLC17A7", "DLG4", "GLS", 
                            "VIM", "NES", "SOX2")) +
  coord_flip()


# confirm some clusters have low counts and low features 
combined_gex@meta.data$log_nCount_RNA <- log(combined_gex@meta.data$nCount_RNA)
combined_gex@meta.data$log_nFeature_RNA <- log(combined_gex@meta.data$nFeature_RNA)

FeaturePlot(combined_gex, 
            features = "log_nCount_RNA", 
            pt.size = 0.2)

FeaturePlot(combined_gex, 
            features = "log_nFeature_RNA", 
            pt.size = 0.2)

FeaturePlot(combined_gex, 
            features = "percent.mt", 
            pt.size = 0.2)

# plot time as the feature, see if times are well-separated
DimPlot(combined_gex,
        cols = c("red", "green", "black"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("clustering timed by 0, 1, 6hr")


# assign group identity to the cells 
combined_gex@meta.data$group.ident <- "unassigned"
combined_gex@meta.data$group.ident[combined_gex@meta.data$orig.ident %in% c("g_2_0_read_gex", "g_2_1_read_gex", "g_2_6_read_gex")] <- "group_2"
combined_gex@meta.data$group.ident[combined_gex@meta.data$orig.ident %in% c("g_8_0_read_gex", "g_8_1_read_gex", "g_8_6_read_gex")] <- "group_8"
#check if there is any cell left out
sum(combined_gex@meta.data$group.ident == "unassigned") # 0 

# already SCTransformed, subset into 3 time points again
resplit_0 <- subset(combined_gex, subset = time.ident == "0hr")
resplit_1 <- subset(combined_gex, subset = time.ident == "1hr")
resplit_6 <- subset(combined_gex, subset = time.ident == "6hr")

# plot at 3 time points the behaviors of two groups to compare and see if they are similar
DimPlot(resplit_0,
        cols = c("cyan", "magenta"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "group.ident",
        label = F) +
  ggtitle("group_2 and group_8 compared, 0hr")

DimPlot(resplit_1,
        cols = c("cyan", "magenta"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "group.ident",
        label = F) +
  ggtitle("group_2 and group_8 compared, 1hr")

DimPlot(resplit_6,
        cols = c("cyan", "magenta"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "group.ident",
        label = F) +
  ggtitle("group_2 and group_8 compared, 6hr")

# assign cell type based on the stacked vlnplot
# Glut_barcodes <- 
#   combined_gex@assays$SCT@counts@Dimnames[[2]][combined_gex@meta.data$seurat_clusters %in% 
#                                               c("2", "4", "5", "6", "7", "8", "10", "12", "14", "17", "19", "22", "23", "24", "25")]
# GABA_barcodes <- 
#   combined_gex@assays$SCT@counts@Dimnames[[2]][combined_gex@meta.data$seurat_clusters %in% 
#                                               c("0", "1", "3", "20")]
# Glut_barcodes <- 
#   combined_gex@assays$SCT@counts@Dimnames[[2]][combined_gex@meta.data$seurat_clusters %in% 
#                                                  c("1", "3", "4", "15", 
#                                                    "19", "5", "6", "8", 
#                                                    "9", "10", "20", "21",
#                                                    "11", "24", "16")]
# GABA_barcodes <- 
#   combined_gex@assays$SCT@counts@Dimnames[[2]][combined_gex@meta.data$seurat_clusters %in% 
#                                                  c("0", "2", "14", "7")]

Glut_barcodes <- 
  combined_gex@assays$SCT@counts@Dimnames[[2]][combined_gex@meta.data$seurat_clusters %in% 
                                                 c("2", "3", "5", "19", 
                                                   "16", "17", "8", "7",
                                                   "9", "12", "13", "18", "22")]
GABA_barcodes <- 
  combined_gex@assays$SCT@counts@Dimnames[[2]][combined_gex@meta.data$seurat_clusters %in% 
                                                 c("1", "11", "6", "4")]

cb_glut <- subset(combined_gex, subset = seurat_clusters %in% 
         c("2", "3", "5", "19", 
           "16", "17", "8", "7",
           "9", "12", "13", "18", "22"))
sum((cb_glut@meta.data$cell.line.ident == "CD_27") & (cb_glut@meta.data$time.ident == "6hr"))
test_gb <- cb_glut@assays$RNA@counts@Dimnames[[2]]

## write out line-type barcodes
write.table(Glut_barcodes,
            file = "../combined_barcodes/group_2_Glut_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)

write.table(GABA_barcodes,
            file = "../combined_barcodes/group_2_GABA_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)



combined_gex@meta.data$cell.type.ident <- "not assigned"
combined_gex@meta.data$cell.type.ident[combined_gex@assays$SCT@counts@Dimnames[[2]] %in% Glut_barcodes] <- "Glut"
combined_gex@meta.data$cell.type.ident[combined_gex@assays$SCT@counts@Dimnames[[2]] %in% GABA_barcodes] <- "GABA"

# already SCTransformed, subset data
combined_glut <- subset(combined_gex, subset = cell.type.ident == "Glut")
combined_GABA <- subset(combined_gex, subset = cell.type.ident == "GABA")

saveRDS(combined_gex, file = "./combined_g2_g8_USE_THIS.rds")


# write a function to plot repeatedly
fgplot_draw <- function(sobj, gene_lst){
  count = 0
  for(gene in gene_lst){
    count = count + 1
    p <- VlnPlot(sobj, 
                 features = gene,
                 group.by = "time.ident",
                 slot = "scale.data",
                 pt.size = 0) + 
      geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1) +
      ylim(c(-1,2.5))
    print(p)
  }
}
# marker genes are checked against homo sapien genome; referenced from Spiegel et al., 2014
markergenes <- c("EGR1", "FOS", "NPAS4", "VGF" )
#markergenes <- c("GAPDH", "RPL8", "ACTB", "BDNF")
lowexpgenes <- c("FOSB", "RGS2", "ADCYAP1", "NRN1")
lowexpgenes <- c("FOSB", "RGS2", "IGF1", "PNOC", "PTHLH")
fgplot_draw(combined_glut, lowexpgenes)
fgplot_draw(combined_GABA, lowexpgenes)

excitatory_markers <- c("EGR1", "FOS", "FOSB", "NPAS4", "RGS2", "VGF", "ADCYAP1", "NRN1")
inhibitory_markers <- c("EGR1", "FOS", "FOSB", "NPAS4", "RGS2", "VGF", "IGF1", "PNOC", "PTHLH")

# create mean and median matrices to store values
glut_mean_mat = matrix(data = NA, nrow = 3, ncol=8)
glut_med_mat = matrix(data = NA, nrow = 3, ncol=8)
GABA_mean_mat = matrix(data = NA, nrow = 3, ncol=9)
GABA_med_mat = matrix(data = NA, nrow = 3, ncol=9)
timelst <- c("0hr", "1hr", "6hr")

# calculate median and mean
for(i in 1:3){
  for(j in 1:8){
    p <- data.frame(val = combined_glut@assays$SCT@scale.data[rownames(combined_gex@assays$SCT@scale.data) %in% excitatory_markers[j], ],
                    time = combined_glut$time.ident,
                    stringsAsFactors = F)
    glut_mean_mat[i, j] <- mean(p$val[p$time %in% timelst[i]])
    glut_med_mat[i, j] <- median(p$val[p$time %in% timelst[i]])
  }
}

for(i in 1:3){
  for(j in 1:9){
    p <- data.frame(val = combined_glut@assays$SCT@scale.data[rownames(combined_gex@assays$SCT@scale.data) %in% inhibitory_markers[j], ],
                    time = combined_glut$time.ident,
                    stringsAsFactors = F)
    GABA_mean_mat[i, j] <- mean(p$val[p$time %in% timelst[i]])
    GABA_med_mat[i, j] <- median(p$val[p$time %in% timelst[i]])
  }
}

rownames(glut_mean_mat) <- c(timelst)
rownames(glut_med_mat) <- c(timelst)
rownames(GABA_mean_mat) <- c(timelst)
rownames(GABA_med_mat) <- c(timelst)
colnames(glut_mean_mat) <- excitatory_markers
colnames(glut_med_mat) <- excitatory_markers
colnames(GABA_mean_mat) <- inhibitory_markers
colnames(GABA_med_mat) <- inhibitory_markers

# do t test
mat0 <- combined_glut@assays$SCT@scale.data[rownames(combined_glut@assays$SCT@scale.data) %in% "VGF", 
                                            combined_glut@meta.data$triple.ident == "0hr"]

for(i in excitatory_markers){
  print("marker gene is: ")
  print(i)
  mat0 <- combined_glut@assays$SCT@scale.data[rownames(combined_glut@assays$SCT@scale.data) %in% i, 
                                              combined_glut@meta.data$time.ident == "0hr"]
  mat1 <- combined_glut@assays$SCT@scale.data[rownames(combined_glut@assays$SCT@scale.data) %in% i, 
                                              combined_glut@meta.data$time.ident == "1hr"]
  mat6 <- combined_glut@assays$SCT@scale.data[rownames(combined_glut@assays$SCT@scale.data) %in% i, 
                                              combined_glut@meta.data$time.ident == "6hr"]
  
  # t test comparing 0 and 1 hrs
  print(t.test(mat0, mat1, 
         paired = F, 
         var.equal = F, 
         alternative = "t"))
  # t test comparing 0 and 6 hrs
  print(t.test(mat0, mat6, 
         paired = F, 
         var.equal = F, 
         alternative = "t"))
  # t test comparing 1 and 6 hrs
  print(t.test(mat1, mat6, 
         paired = F, 
         var.equal = F, 
         alternative = "t"))
}

for(i in inhibitory_markers){
  print("marker gene is: ")
  print(i)
  mat0 <- combined_GABA@assays$SCT@scale.data[rownames(combined_GABA@assays$SCT@scale.data) %in% i, 
                                              combined_GABA@meta.data$time.ident == "0hr"]
  mat1 <- combined_GABA@assays$SCT@scale.data[rownames(combined_GABA@assays$SCT@scale.data) %in% i, 
                                              combined_GABA@meta.data$time.ident == "1hr"]
  mat6 <- combined_GABA@assays$SCT@scale.data[rownames(combined_GABA@assays$SCT@scale.data) %in% i, 
                                              combined_GABA@meta.data$time.ident == "6hr"]
  
  # t test comparing 0 and 1 hrs
  print(t.test(mat0, mat1, 
               paired = F, 
               var.equal = F, 
               alternative = "t"))
  # t test comparing 0 and 6 hrs
  print(t.test(mat0, mat6, 
               paired = F, 
               var.equal = F, 
               alternative = "t"))
  # t test comparing 1 and 6 hrs
  print(t.test(mat1, mat6, 
               paired = F, 
               var.equal = F, 
               alternative = "t"))
}
# combined_0hr_SCTmat <- combined_glut@assays$SCT@scale.data[rownames(combined_glut@assays$SCT@scale.data) %in% "PNOC", 
#                                                           combined_glut@meta.data$triple.ident == "0hr"]
# combined_1hr_SCTmat <- combined_glut@assays$SCT@scale.data[rownames(combined_glut@assays$SCT@scale.data) %in% "PNOC",
#                                                           combined_glut@meta.data$triple.ident == "1hr"]
# combined_6hr_SCTmat <- combined_glut@assays$SCT@scale.data[rownames(combined_glut@assays$SCT@scale.data) %in% "PNOC",
#                                                           combined_glut@meta.data$triple.ident == "6hr"]
# 
# # t test comparing 0 and 1 hrs
# t.test(combined_0hr_SCTmat, combined_0hr_SCTmat, 
#        paired = F, 
#        var.equal = F, 
#        alternative = "t")
# # t test comparing 0 and 6 hrs
# t.test(combined_0hr_SCTmat, combined_6hr_SCTmat, 
#        paired = F, 
#        var.equal = F, 
#        alternative = "t")
# # t test comparing 1 and 6 hrs
# t.test(combined_1hr_SCTmat, combined_6hr_SCTmat, 
#        paired = F, 
#        var.equal = F, 
#        alternative = "t")
