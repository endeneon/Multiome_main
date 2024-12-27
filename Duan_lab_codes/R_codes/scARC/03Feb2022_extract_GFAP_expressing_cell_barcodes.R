# Chuxuan Li 02/03/2022
# use libraries 36 and 49 that mapped onto mixed human and mouse genome,
#check which cells express GFAP and extract their barcodes


library(Seurat)
library(stringr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)

library(future)

# set threads and parallelization
plan("multisession", workers = 1)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 107374182400)

# load data ####
setwd("/data/FASTQ/Duan_Project_022_Reseq/")
# note this h5 file contains both atac-seq and gex information
# Read files, separate the ATACseq data from the .h5 matrix
dirlist <- list.files(path = ".", 
                      pattern = "_hybrid$",
                      recursive = F,
                      include.dirs = T)
h5list <- list.files(path = dirlist, 
                     pattern = "filtered_feature_bc_matrix.h5",
                     recursive = T,
                     include.dirs = T)
pathlist <- rep_len("", length(dirlist))
for (i in 1:length(dirlist)){
  pathlist[i] <- paste(".", dirlist[i], h5list[i], sep = "/")
}
objlist <- vector(mode = "list", length = length(h5list))
transformed_lst <- vector(mode = "list", length(h5list))

# normalization
for (i in 1:length(objlist)){
  h5file <- Read10X_h5(filename = pathlist[i])
  print(str_extract(string = pathlist[i], 
                    pattern = "[0-9][0-9]_[0-6]"))
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(string = h5list[i], 
                                                  pattern = "[0-9][0-9]_[0-6]"))
  obj$lib.ident <- str_extract(string = h5list[i], 
                               pattern = "[0-9][0-9]_[0-6]")
  print(paste0(i, " number of genes: ", nrow(obj), 
               ", number of cells: ", ncol(obj)))
  obj <- 
    PercentageFeatureSet(obj,
                         pattern = c("^MT-"),
                         col.name = "percent.mt")
  objlist[[i]] <- FindVariableFeatures(obj, nfeatures = 3000)
  transformed_lst[[i]] <- SCTransform(objlist[[i]], 
                                      vars.to.regress = "percent.mt",
                                      method = "glmGamPoi",
                                      variable.features.n = 8000,
                                      seed.use = 115,
                                      verbose = T)
}

setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat")
features <- SelectIntegrationFeatures(object.list = transformed_lst, 
                                      nfeatures = 3000, 
                                      fvf.nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = transformed_lst, 
                                  anchor.features = features,
                                  dims = 1:50)
hybrid_combined <- IntegrateData(anchorset = anchors)

hybrid_combined <- ScaleData(hybrid_combined)
hybrid_combined <- RunPCA(hybrid_combined, 
                     verbose = T, 
                     seed.use = 11)
hybrid_combined <- RunUMAP(hybrid_combined, 
                      reduction = "pca", 
                      dims = 1:30,
                      seed.use = 11)
hybrid_combined <- FindNeighbors(hybrid_combined, 
                            reduction = "pca", 
                            dims = 1:30)
hybrid_combined <- FindClusters(hybrid_combined, 
                           resolution = 0.5,
                           random.seed = 11)

save("hybrid_combined", file = "human_rat_hybrid_integrated_clustered_obj.RData")
DimPlot(hybrid_combined, 
        label = T) +
  NoLegend() +
  ggtitle("libraries 36, 49 mapped onto mixed genome")

DefaultAssay(hybrid_combined) <- "RNA"
FeaturePlot(hybrid_combined, features = "Gfap", max.cutoff = 10)
FeaturePlot(hybrid_combined, features = "Col1a1", max.cutoff = 400)
FeaturePlot(hybrid_combined, features = "Actb", max.cutoff = 50)
FeaturePlot(hybrid_combined, features = "MAP2", max.cutoff = 40)

FeaturePlot(hybrid_combined, features = c("GAD1", "SLC17A6"), max.cutoff = 15)
FeaturePlot(hybrid_combined, features = c("NEFM", "CUX2"), max.cutoff = 15)

sum(hybrid_combined@assays$RNA@counts@Dimnames[[1]] == "Gapdh")
sum(hybrid_combined@assays$RNA@counts@Dimnames[[1]] == "Col1a1")
sum(hybrid_combined@assays$RNA@counts[hybrid_combined@assays$RNA@counts@Dimnames[[1]] == "Gfap", ] > 0) #1809
write.table(x = 
              hybrid_combined@assays$RNA@counts@Dimnames[[2]][hybrid_combined$seurat_clusters == "2"],
            file = "rat_barcodes_cluster2.txt",
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = 
              hybrid_combined@assays$RNA@counts@Dimnames[[2]][hybrid_combined$seurat_clusters == "6"],
            file = "rat_barcodes_cluster6.txt",
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = 
              hybrid_combined@assays$RNA@counts@Dimnames[[2]][hybrid_combined$seurat_clusters %in% c("2", "6")],
            file = "rat_barcodes_cluster2_6.txt",
            quote = F, sep = "\t", row.names = F, col.names = F)

# load GRCh38 mapped only object
RNAseq_integrated_labeled$rat.ident <- "human"
head(RNAseq_integrated_labeled@assays$RNA@counts@Dimnames[[2]])
head(hybrid_combined@assays$RNA@counts@Dimnames[[2]])

rat_barcodes_cluster2 <- read_delim("rat_barcodes_cluster2.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    col_names = FALSE, trim_ws = TRUE, quote = "")
rat_barcodes_cluster2 <- unlist(rat_barcodes_cluster2)
rat_barcodes_cluster6 <- read_delim("rat_barcodes_cluster6.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    col_names = FALSE, trim_ws = TRUE)
rat_barcodes_cluster6 <- unlist(rat_barcodes_cluster6)

RNAseq_integrated_labeled$rat.ident[str_sub(RNAseq_integrated_labeled@assays$RNA@counts@Dimnames[[2]], end = -3L)
                                    %in% 
                                      str_sub(rat_barcodes_cluster2, end = -3L)] <- "rat_c2"
RNAseq_integrated_labeled$rat.ident[str_sub(RNAseq_integrated_labeled@assays$RNA@counts@Dimnames[[2]], end = -3L)
                                    %in% 
                                      str_sub(rat_barcodes_cluster6, end = -3L)] <- "rat_c6"
unique(RNAseq_integrated_labeled$rat.ident)
DimPlot(RNAseq_integrated_labeled, 
        group.by = "rat.ident", 
        cols = c("transparent", "transparent", "red"))
DimPlot(RNAseq_integrated_labeled, label = T) + 
  NoLegend() + 
  ggtitle("labeled by cell type")
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled$broad.cell.type == "NPC"] <- "rat_origin"
Idents(RNAseq_integrated_labeled) <- "broad.cell.type"
sum(RNAseq_integrated_labeled$broad.cell.type == "rat_origin")
# separate human from rat ####
RNAobj_sans_rat <- subset(RNAseq_integrated_labeled, broad.cell.type %in% c("GABA", "NEFM_neg_glut", "NEFM_pos_glut"))
DimPlot(RNAobj_sans_rat, label = T) + 
  NoLegend() + 
  ggtitle("labeled by cell type")
RNAobj_sans_rat <-  RunPCA(RNAobj_sans_rat, 
                       verbose = T, 
                       seed.use = 11)
RNAobj_sans_rat <- RunUMAP(RNAobj_sans_rat, 
                      reduction = "pca", 
                      dims = 1:30,
                      seed.use = 11)
RNAobj_sans_rat <- FindNeighbors(RNAobj_sans_rat, 
                            reduction = "pca", 
                            dims = 1:30)
RNAobj_sans_rat <- FindClusters(RNAobj_sans_rat, 
                           resolution = 0.5,
                           random.seed = 11)
DimPlot(RNAobj_sans_rat, label = T) + 
  NoLegend() + 
  ggtitle("removed rat astrocytes, then re-clustered")
DefaultAssay(RNAobj_sans_rat) <- "RNA"
FeaturePlot(RNAobj_sans_rat, features = c("GAD1", "SLC17A6"), max.cutoff = 15)
FeaturePlot(RNAobj_sans_rat, features = c("NEFM", "CUX2"), max.cutoff = 15)

