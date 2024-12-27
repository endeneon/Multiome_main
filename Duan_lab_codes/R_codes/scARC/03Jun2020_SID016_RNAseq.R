# 03 Jun 2020 Siwei
# Run multi-timepoint RNAseq use SID0, 1, and 6

# init
library(Seurat)
library(future)
library(readr)
library(sctransform)
library(stringr)
library(ggplot2)
library(cowplot)
library(patchwork)

# set parameters
options(future.globals.maxSize = 6442450944)

# setup multithread
plan("multisession", workers = 4) # should not use "multicore" here
plan("multisession", workers = 1) 

# use SeuratData for reference ########
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
InstallData("ifnb")

data("ifnb")
ifnb <- ifnb
ifnb.list <- SplitObject(ifnb, split.by = "stim")

# load scRNAseq data ###########
# read in SID0 data
SID0.RNAseq.data <- Read10X("../scRNA/SID0_scRNA/outs/filtered_feature_bc_matrix/")
SID0_demux_RNA_index <- read_table2("~/Galaxy_Drive/May15_Duan_017/demultiplex/SID0_scRNA_099_ASE/SID0_scRNA_099_ASE_cell_line_index.txt")
SID0.RNAseq.data <- SID0.RNAseq.data[ , colnames(SID0.RNAseq.data) %in% # save memory
                                          SID0_demux_RNA_index$BARCODE]

SID1.RNAseq.data <- Read10X("../scRNA/SID1_scRNA/outs/filtered_feature_bc_matrix/")
# SID1_demux_RNA_index <- read_table2("~/Galaxy_Drive/May15_Duan_017/demultiplex/SID1_scRNA_ASE_final.txt")
SID1_demux_RNA_index <- read_table2("~/Galaxy_Drive/May15_Duan_017/demultiplex/SID1_scRNA_ASE_099_final.txt")
SID1.RNAseq.data <- SID1.RNAseq.data[ , colnames(SID1.RNAseq.data) %in% # save memory
                                       SID1_demux_RNA_index$BARCODE]

SID6.RNAseq.data <- Read10X("../scRNA/SID6_scRNA/outs/filtered_feature_bc_matrix/")
SID6_demux_RNA_index <- read_table2("~/Galaxy_Drive/May15_Duan_017/demultiplex/SID6_scRNA_ASE_final.txt")
SID6.RNAseq.data <- SID6.RNAseq.data[ , colnames(SID6.RNAseq.data) %in% # save memory
                                       SID6_demux_RNA_index$BARCODE]

# use SCTransform wrapper, note that need to use SCT assay now
# remove confounding factors: MT
# make list of Seurat objects

RNASeq.object.list <- vector(mode = "list", length = 3)
RNASeq.object.list[[1]] <- CreateSeuratObject(counts = SID0.RNAseq.data, 
                                             project = "SID0")
RNASeq.object.list[[2]] <- CreateSeuratObject(counts = SID1.RNAseq.data, 
                                             project = "SID1")
RNASeq.object.list[[3]] <- CreateSeuratObject(counts = SID6.RNAseq.data, 
                                             project = "SID6")
names(RNASeq.object.list) <- c("SID0", "SID1", "SID6")
# remove unused objects to save memory
rm(SID0.RNAseq.data, SID1.RNAseq.data, SID6.RNAseq.data) 
   # SID0_demux_RNA_index, SID1_demux_RNA_index, SID6_demux_RNA_index)
# save.image(file = "RNASeq.multi.point.RData")

# assign individual identity to each cell
RNASeq.object.list[[1]] <- AddMetaData(RNASeq.object.list[[1]], 
                                       metadata = SID0_demux_RNA_index$BEST[
                                         match(colnames(RNASeq.object.list[[1]]), 
                                               SID0_demux_RNA_index$BARCODE)], 
                                       col.name = "indiv.identity")
RNASeq.object.list[[2]] <- AddMetaData(RNASeq.object.list[[2]], 
                                       metadata = SID1_demux_RNA_index$BEST[
                                         match(colnames(RNASeq.object.list[[2]]), 
                                               SID1_demux_RNA_index$BARCODE)], 
                                       col.name = "indiv.identity")
RNASeq.object.list[[3]] <- AddMetaData(RNASeq.object.list[[3]], 
                                       metadata = SID6_demux_RNA_index$BEST[
                                         match(colnames(RNASeq.object.list[[3]]), 
                                               SID6_demux_RNA_index$BARCODE)], 
                                       col.name = "indiv.identity")
rm(SID0_demux_RNA_index, SID1_demux_RNA_index, SID6_demux_RNA_index)

# lapply SCTransform to RNAseq objects
RNASeq.object.list <- lapply(X = RNASeq.object.list, FUN = function(x) {
  print(x)
  # store mitochondrial percentage in object meta data
  x <- PercentageFeatureSet(x, pattern = "^MT-|^Mt-", # need to include rat mt
                            col.name = "percent.mt") 
  x <- SCTransform(x, 
                   vars.to.regress = "percent.mt", 
                   verbose = T)
  
})
save.image(file = "RNASeq.multi.point.RData")

# perform integration, use SCT to normalise ##########
# The expected workflow for integratinge assays produced by SCTransform is 
# SelectIntegrationFeatures -> PrepSCTIntegration -> FindIntegrationAnchors.
# Note this is different from the traditional integrated analysis 
# documented on the Satija lab website. This protocol is extracted from
# the help document of Seurat 3 PrepSCTIntegration()

# select integration fetures
RNASeq.features <- SelectIntegrationFeatures(object.list = RNASeq.object.list, 
                                             nfeatures = 6000, 
                                             assay = rep("SCT", 3),
                                             verbose = T)
# integrate features
RNASeq.object.list <- PrepSCTIntegration(object.list = RNASeq.object.list, 
                                         anchor.features = RNASeq.features, 
                                         assay = "SCT")
# find anchors and integrate
RNASeq.anchors <- FindIntegrationAnchors(object.list = RNASeq.object.list, 
                                         dims = 1:20,
                                         # reference = c(1, 0, 0), # use SID0 as reference
                                         reduction = "cca", 
                                         anchor.features = RNASeq.features,
                                         normalization.method = "SCT", 
                                         verbose = T)
RNASeq.integrated <- IntegrateData(RNASeq.anchors)
# 
# RNASeq.combined <- IntegrateData(anchorset = RNASeq.object.list, 
#                                  dims = 1:20, 
#                                  normalization.method = "SCT", 
#                                  verbose = T)
save.image(file = "RNASeq.multi.point.RData")


# perform integrated analysis #######
DefaultAssay(RNASeq.integrated) <- "integrated"

# run standard workflow to generate visualisation + clustering
# save a backup copy here to see the difference of applying ScaleData()
# since the RNASeq.integrated was generated by SCTransform() it may not
# be necessary to apply ScaleData again?
# RNASeq.integrated.backup <- RNASeq.integrated
# RNASeq.integrated <- RNASeq.integrated.backup

RNASeq.integrated <- ScaleData(RNASeq.integrated, verbose = T)

RNASeq.integrated <- RunPCA(RNASeq.integrated, 
                            npcs = 30, verbose = T)
RNASeq.integrated <- RunUMAP(RNASeq.integrated, 
                             reduction = "pca", 
                             dims = 1:30, verbose = T)
RNASeq.integrated <- FindNeighbors(RNASeq.integrated, 
                                   reduction = "pca", 
                                   dims = 1:30, verbose = T)
RNASeq.integrated <- FindClusters(RNASeq.integrated, # need 14 communities
                                  resolution = 0.25, verbose = T)

# set rat astrocytes as separate entity
RNASeq.integrated$indiv.identity[RNASeq.integrated$seurat_clusters %in% "7"] <- "Rat astrocytes"
# change individual names
RNASeq.integrated$indiv.identity[RNASeq.integrated$indiv.identity %in% "SNG-CN_line_05"] <- "Line 05"
RNASeq.integrated$indiv.identity[RNASeq.integrated$indiv.identity %in% "SNG-CN_line_14"] <- "Line 14"
RNASeq.integrated$indiv.identity[RNASeq.integrated$indiv.identity %in% "SNG-iPS_line_15"] <- "Line 15"


# make plots
time_series_plot <- DimPlot(RNASeq.integrated, 
                            reduction = "umap", 
                            group.by = "orig.ident") 
time_series_plot[[1]]$layers[[1]]$aes_params$alpha <- 0.7
time_series_plot

DimPlot(RNASeq.integrated, reduction = "umap", label = T)

DimPlot(RNASeq.integrated, 
        reduction = "umap", 
        group.by = "indiv.identity")

DimPlot(RNASeq.integrated, 
        reduction = "umap", 
        group.by = c("orig.ident", "indiv.identity")) 

DimPlot(RNASeq.integrated, 
        reduction = "umap", 
        group.by = "indiv.identity", 
        split.by = "orig.ident", 
        cols = c("darkred", "orange", "darkblue", "black")) 


DefaultAssay(RNASeq.integrated) <- "SCT"
DefaultAssay(RNASeq.integrated) <- "integrated"
FeaturePlot(RNASeq.integrated, 
            features = c("RPL10"), 
            split.by = "orig.ident", 
            reduction = "umap")
FeaturePlot(RNASeq.integrated, 
            features = c("FOS"), 
            # cells = RNASeq.integrated$indiv.identity %in% "SNG-CN_line_05",
            split.by = c("orig.ident"), 
            reduction = "umap")
RNASeq.integrated.line.subset <- SubsetData(RNASeq.integrated, 
                                            cells = RNASeq.integrated$indiv.identity == "Line 14")
DefaultAssay(RNASeq.integrated.line.subset) <- "integrated"
FeaturePlot(RNASeq.integrated.line.subset, 
            features = c("FOS"), 
            # cells = RNASeq.integrated$indiv.identity %in% "SNG-CN_line_05",
            split.by = c("orig.ident"),
            reduction = "umap")



DefaultAssay(RNASeq.integrated) <- "RNA"
StackedVlnPlot(obj = RNASeq.integrated, 
               features = c( "GAD1", "GAD2", "PNOC", "HOXA5",
                            "SLC17A6", 
                            "DLG4",
                            "GLS",
                            "FOS",
                            "VIM", "NES", "SOX2", 
                            "S100b", "Slc1a3", "Gfap")) +
  coord_flip()

save.image(file = "RNASeq.multi.point.RData")


# all.gene.expression <- AverageExpression(object = RNASeq.integrated.line.subset, 
#                                          assays = "SCT",
#                                          features = c("GAPDH", "FOS", "BDNF"))
# all.gene.expression.1 <- AverageExpression(object = RNASeq.integrated.line.subset, 
#                                          assays = "SCT",
#                                          features = c("GAPDH", "FOS", "BDNF"))
# all.gene.expression.2 <- AverageExpression(object = RNASeq.integrated.line.subset, 
#                                            assays = "SCT",
#                                            features = c("GAPDH", "FOS", "BDNF"))
View(all.gene.expression.2[[1]])

# Subset each cell line to make DotPlot split by different time points
library(ggpubr)
# SID0
RNASeq.integrated.line.05 <- SubsetData(RNASeq.integrated, 
                                            cells = RNASeq.integrated$indiv.identity == "Line 05")
Idents(RNASeq.integrated.line.05)
DefaultAssay(RNASeq.integrated.line.05) <- "SCT"
line.05.dotplot <- DotPlot(RNASeq.integrated.line.05, features = c("FOS", "BDNF", "NPAS4", 
                                                                   "SP4", "TCF4"), 
                           cols = c("blue", "blue", "blue"), 
                           dot.scale = 6, 
                           split.by = "orig.ident") +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Line 05")

RNASeq.integrated.line.14 <- SubsetData(RNASeq.integrated, 
                                        cells = RNASeq.integrated$indiv.identity == "Line 14")
Idents(RNASeq.integrated.line.14)
DefaultAssay(RNASeq.integrated.line.14) <- "SCT"
line.14.dotplot <- DotPlot(RNASeq.integrated.line.14, features = c("FOS", "BDNF", "NPAS4", 
                                                                   "SP4", "TCF4"), 
                           cols = c("blue", "blue", "blue"), 
                           dot.scale = 6, 
                           split.by = "orig.ident") +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Line 14")

RNASeq.integrated.line.15 <- SubsetData(RNASeq.integrated, 
                                        cells = RNASeq.integrated$indiv.identity == "Line 15")
Idents(RNASeq.integrated.line.15)
DefaultAssay(RNASeq.integrated.line.15) <- "SCT"
line.15.dotplot <- DotPlot(RNASeq.integrated.line.15, features = c("FOS", "BDNF", "NPAS4", 
                                                                   "SP4", "TCF4"), 
                           cols = c("blue", "blue", "blue"), 
                           dot.scale = 6, 
                           split.by = "orig.ident") +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Line 15")

ggarrange(line.05.dotplot, 
          line.14.dotplot,
          line.15.dotplot, 
          ncol = 1)

rm(RNASeq.integrated.line.05)
rm(RNASeq.integrated.line.14)
rm(RNASeq.integrated.line.15)

#### PCA plot for 3x individuals x 3 time points x 2/3 cell types #####
# init constants

RNASeq.variable.features <- FindAllMarkers(RNASeq.integrated, assay = "SCT")
RNA.features.3celltypes <- unique(RNASeq.variable.features$gene)

names(new.cluster.ids) <- levels(RNASeq.integrated)
RNASeq.integrated.new.cluster.ids <- RenameIdents(RNASeq.integrated, new.cluster.ids)
RNASeq.integrated.new.cluster.ids <- AddMetaData(RNASeq.integrated.new.cluster.ids, 
                                                 metadata = Idents(RNASeq.integrated.new.cluster.ids), 
                                                 col.name = "cell.type.idents")
RNASeq.integrated.new.cluster.ids <- RNASeq.integrated.new.cluster.ids[ , RNASeq.integrated.new.cluster.ids$cell.type.idents %in% 
                                                                          c("Glut", "GABA", "NPC")] # drop all other cell types
RNASeq.variable.features <- FindAllMarkers(RNASeq.integrated.new.cluster.ids, assay = "SCT")
# Idents(RNASeq.integrated.subset.individual.timepoint)
new.cluster.ids <- c("Glut", "GABA", "Glut", "NPC",
                     "NPC", "Trans", "GABA", "Astro",
                     "Trans", "GABA", "Tail", "Tail",
                     "Tail", "Tail")
## 

RNASeq.integrated <- RNASeq.integrated.new.cluster.ids

for (indiv.id in c("Line 05", "Line 14", "Line 15")) { # cycle through individual ids
  print(indiv.id)
  RNASeq.integrated.subset.individual <- SubsetData(RNASeq.integrated, 
                                                    cells = RNASeq.integrated$indiv.identity == indiv.id)
  # attach metadata
  RNASeq.integrated.subset.individual <- AddMetaData(RNASeq.integrated.subset.individual, 
                                                     metadata = RNASeq.integrated$orig.ident[
                                                       RNASeq.integrated$indiv.identity %in% indiv.id], 
                                                     col.name = "orig.ident")
  RNASeq.integrated.subset.individual <- AddMetaData(RNASeq.integrated.subset.individual, 
                                                     metadata = RNASeq.integrated$seurat_clusters[
                                                       RNASeq.integrated$indiv.identity %in% indiv.id], 
                                                     col.name = "seurat_clusters")
  
  for (time.points in c("SID0", "SID1", "SID6")) {
    print(time.points)
    DefaultAssay(RNASeq.integrated.subset.individual) <- "SCT"
    RNASeq.integrated.subset.individual.timepoint <- 
      SubsetData(RNASeq.integrated.subset.individual, 
                 cells = RNASeq.integrated.subset.individual$orig.ident == "SID0")
    # attach metadata
    RNASeq.integrated.subset.individual.timepoint <- AddMetaData(RNASeq.integrated.subset.individual.timepoint, 
                                                                 metadata = RNASeq.integrated.subset.individual$orig.ident[
                                                                   RNASeq.integrated.subset.individual$orig.ident == "SID0"
                                                                   ], 
                                                                 col.name = "orig.ident")
    RNASeq.integrated.subset.individual.timepoint <- AddMetaData(RNASeq.integrated.subset.individual.timepoint, 
                                                                 metadata = RNASeq.integrated.subset.individual$seurat_clusters[
                                                                   RNASeq.integrated.subset.individual$orig.ident == "SID0"
                                                                   ], 
                                                                 col.name = "seurat_clusters")
    # attach new cluster identities
    
    # names(new.cluster.ids) <- levels(RNASeq.integrated.subset.individual.timepoint)
    RNASeq.integrated.subset.individual.timepoint <- RenameIdents(RNASeq.integrated.subset.individual.timepoint, 
                                                                  new.cluster.ids)
    DefaultAssay(RNASeq.integrated.subset.individual.timepoint) <- "SCT"
    exp.matrix <- AverageExpression(RNASeq.integrated.subset.individual.timepoint)[[2]]
    exp.matrix <- exp.matrix[rownames(exp.matrix) %in%
                               unique(RNASeq.variable.features$gene), ]
    # combine data frame
    
  }
}




# RNA.features.3celltypes <- rownames(RNASeq.variable.features)[RNASeq.variable.features$cluster %in% 
#                                                            c("Glut", "GABA", "NPC")]
# RNA.features.3celltypes <- RNASeq.variable.features$gene[RNASeq.variable.features$cell.type.idents %in% 
#                                                                 c("Glut", "GABA", "NPC")]

# RNASeq.integrated <- RNASeq.integrated.backup
RNASeq.integrated <- RNASeq.integrated.new.cluster.ids
RNASeq.integrated.plot <- RNASeq.integrated[, RNASeq.integrated$orig.ident == "SID0"]
DefaultAssay(RNASeq.integrated.plot) <- "RNA"
RNASeq.integrated.plot <- PercentageFeatureSet(RNASeq.integrated.plot, pattern = "^MT-|^Mt-", # need to include rat mt
                                    col.name = "percent.mt") 
RNASeq.integrated.plot <- SCTransform(RNASeq.integrated.plot, 
                           vars.to.regress = "percent.mt", 
                           verbose = T)
## make PCA and UMAP data and store
RNASeq.integrated.plot <- RunPCA(RNASeq.integrated.plot, verbose = T)
RNASeq.integrated.plot <- RunUMAP(RNASeq.integrated.plot, dims = 1:30, verbose = T)

RNASeq.integrated.plot <- FindNeighbors(RNASeq.integrated.plot, 
                             # k.param = 15, # set k here
                             dims = 1:30, 
                             verbose = T)
RNASeq.integrated.plot <- FindClusters(RNASeq.integrated.plot, 
                            resolution = 0.5, # adjust cluster number here, 0.5 ~ 14 clusters
                            verbose = T)

DimPlot(RNASeq.integrated.plot, reduction = "umap", label = T)

