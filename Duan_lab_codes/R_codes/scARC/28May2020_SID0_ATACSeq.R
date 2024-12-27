# 28 May 2020 Siwei
# Run single-timepoint ATACseq and RNASeq use SID0

# init
library(Seurat)
library(future)
library(readr)
library(sctransform)
library(stringr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(hdf5r)

# set parameters
options(future.globals.maxSize = 6442450944)
options(future.fork.enable = T) # force enable forks


# setup multithread
plan("multicore", workers = 4) # force multicore


plan("multisession", workers = 8) # should not use "multicore" here
plan("multisession", workers = 1) 
future::resetWorkers(x = 8)
# devtools::install_github("timoast/signac", ref = "develop")

# read in SID0 data
SID0.RNAseq.data <- Read10X("../scRNA/SID0_scRNA/outs/filtered_feature_bc_matrix/")
SID0.ATACseq.data <- Read10X_h5("../scATAC/SID0_scATAC/outs/filtered_peak_bc_matrix.h5")
SID0_demux_RNA_index <- read_table2("~/Galaxy_Drive/May15_Duan_017/demultiplex/SID0_scRNA/SID0_cell_line_index.txt")
SID0_demux_ATAC_index <- read_table2("~/Galaxy_Drive/May15_Duan_017/demultiplex/SID0_ATAC/SID0_ATAC_index.txt")

# do not clean up cells now since the singlecell.csv still contains all cells
# peak_count_matrix <- CreateGeneActivityMatrix(
#   peak.matrix = SID0.ATACseq.data, 
#   annotation.file = "~/2TB/Databases/hg38_rn6_hybrid/gencode.v30.Rn.6.gtf", 
#   seq.levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
#                  "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
#                  "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
#                  "chr19", "chr20", "chr21", "chr22", #human
#                  c(1:20)), #rat
#   upstream = 5000, 
#   keep.sparse = T,
#   include.body = F, verbose = T) 

# cleanup data
SID0_cleaned_RNAseq <- SID0.RNAseq.data[colnames(SID0.RNAseq.data) %in% 
                                          SID0_demux_RNA_index$BARCODE, ]
SID0_cleaned_ATACseq <- SID0.ATACseq.data[colnames(SID0.ATACseq.data) %in% 
                                            SID0_demux_ATAC_index$BARCODE, ]
# Try to load the cleaned count matrix since the raw matrix is too large
peak_count_matrix <- CreateGeneActivityMatrix(
  peak.matrix = SID0_cleaned_ATACseq,
  annotation.file = "~/2TB/Databases/hg38_rn6_hybrid/gencode.v30.Rn.6.gtf",
  seq.levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                 "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                 "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                 "chr19", "chr20", "chr21", "chr22" #human
                 , c(1:20)
                 ), #rat
  upstream = 2000,
  keep.sparse = T,
  include.body = F, verbose = T)

# Set up the ATAC object
SID0.ATAC.object <- CreateSeuratObject(counts = SID0_cleaned_ATACseq, 
                                       assay = "ATAC", 
                                       project = "10x_ATAC")
SID0.ATAC.object[["ACTIVITY"]] <- CreateAssayObject(counts = peak_count_matrix)
meta <- read.table("../scATAC/SID0_scATAC/outs/singlecell.csv", 
                   sep = ",", header = T, row.names = 1, stringsAsFactors = F)
meta <- meta[colnames(SID0.ATAC.object), ] # only take metadata of identified cells

# preprocess data
DefaultAssay(SID0.ATAC.object) <- "ACTIVITY"
SID0.ATAC.object <- FindVariableFeatures(SID0.ATAC.object)
SID0.ATAC.object <- NormalizeData(SID0.ATAC.object)
SID0.ATAC.object <- ScaleData(SID0.ATAC.object)

# use peaks that have at least 100 reads per cell
DefaultAssay(SID0.ATAC.object) <- "ATAC"
VariableFeatures(SID0.ATAC.object) <- names(which(Matrix::rowSums(SID0.ATAC.object) > 100))

SID0.ATAC.object <- RunLSI(SID0.ATAC.object, n = 50, scale.max = NULL)
SID0.ATAC.object <- RunUMAP(SID0.ATAC.object, reduction = "lsi", dims = 1:50)

# load corresponding RNASeq data
SID0.RNA.object <- readRDS(file = "SID0.RNASeq.RData")
SID0.RNA.object$tech <- "rna"

# make a combined map
p1 <- DimPlot(SID0.ATAC.object, reduction = "umap") +
  NoLegend() +
  ggtitle("SID0 scATAC-seq")
p2 <- DimPlot(SID0.RNA.object, 
              label = T, repel = T) +
  NoLegend() +
  ggtitle("SID0 scRNA-seq")
p1 + p2

DimPlot(SID0.ATAC.object, 
        reduction = "umap", 
        label = T, repel = T) +
  NoLegend() +
  ggtitle("SID0 scATAC-seq")

# identify corresponding anchors and transfer
SID0.transfer.anchors <- FindTransferAnchors(reference = SID0.RNA.object, 
                                             query = SID0.ATAC.object, 
                                             features = VariableFeatures(object = SID0.RNA.object), 
                                             reference.assay = "RNA", 
                                             query.assay = "ACTIVITY", 
                                             reduction = "cca")
# predict cell types and check efficiency
SID0.celltype.predictions <- TransferData(anchorset = SID0.transfer.anchors, 
                                          refdata = SID0.RNA.object$seurat_clusters, 
                                          weight.reduction = SID0.ATAC.object[["lsi"]])
SID0.ATAC.object <- AddMetaData(SID0.ATAC.object, 
                                metadata = SID0.celltype.predictions)
hist(SID0.ATAC.object$prediction.score.max) # check prediction score
table(SID0.ATAC.object$prediction.score.max > 0.5)

# subset ATAC-seq cells and only keep score > 0.5
SID0.ATAC.filtered <- subset(SID0.ATAC.object, 
                             subset = prediction.score.max > 0.5)
# match group colour
SID0.ATAC.filtered$predicted.id <- factor(SID0.ATAC.filtered$predicted.id, 
                                          levels = levels(SID0.RNA.object))

# make plot
p1 <- DimPlot(SID0.ATAC.filtered, group.by = "predicted.id", label = T, repel = T) +
  ggtitle("SID0 scATAC-seq cells") +
  scale_colour_hue(drop = F) +
  NoLegend()
p2 <- DimPlot(SID0.RNA.object, 
              label = T, repel = T) +
  NoLegend() +
  ggtitle("SID0 scRNA-seq cells")
p1 + p2


DimPlot(SID0.ATAC.filtered, 
        group.by = "predicted.id", 
        label = T, repel = T) +
  ggtitle("SID0 scATAC-seq cells") +
  scale_colour_hue(drop = F) +
  NoLegend()

#### plot cell line identities ######
SID0_umap <- FetchData(SID0.ATAC.filtered, vars = "UMAP_1")
# SID0_umap$UMAP_1 <- FetchData(SID0.ATAC.object, vars = "UMAP_1")
SID0_umap$UMAP_2 <- FetchData(SID0.ATAC.filtered, vars = "UMAP_2")
SID0_umap$UMAP_2 <- SID0_umap$UMAP_2$UMAP_2
SID0_umap$clusters <- FetchData(SID0.ATAC.filtered, vars = "predicted.id")
SID0_umap$clusters <- SID0_umap$clusters$predicted.id
SID0_umap$barcodes <- rownames(SID0.ATAC.filtered@reductions$umap)


SID0_cell_individuals <- base::merge(x = SID0_umap, 
                                     y = SID0_demux_ATAC_index, 
                                     by.x = "barcodes", 
                                     by.y = "BARCODE")


SID0_cell_individuals[SID0_cell_individuals$clusters %in% c(7), 8] <- "Rat Astrocytes"


ggplot(SID0_cell_individuals, aes(x = UMAP_1, y = UMAP_2, 
                                  colour = BEST)) +
  geom_point(size = 0.5) +
  scale_colour_manual(name = "Cell type", 
                      values = c("black", "green4", "darkblue")) +
  ggtitle("Unstimulated, scATAC-seq") +
  theme_classic()


##### separate glut, NPC, and GABA groups
SID0_NGN2_Glut_barcodes <- SID0_cell_individuals[SID0_cell_individuals$clusters %in% c(1,3,8,9,12), 1]
SID0_GABA_barcodes <- SID0_cell_individuals[SID0_cell_individuals$clusters %in% c(0,4,6), 1] 
SID0_NPC_barcodes <- SID0_cell_individuals[SID0_cell_individuals$clusters %in% c(2,5,13), 1] 

# separate lines
SID0_NGN2_Glut_barcodes_14 <- SID0_cell_individuals[(SID0_cell_individuals$clusters %in% c(1,3,8,9,12) &
                                                    (SID0_cell_individuals$BEST %in% "SNG-Line_14_NG_11923")), 1]
SID0_NGN2_Glut_barcodes_15 <- SID0_cell_individuals[(SID0_cell_individuals$clusters %in% c(1,3,8,9,12) &
                                                       (SID0_cell_individuals$BEST %in% "SNG-Line_15_NG_55047")), 1]
SID0_GABA_barcodes_14 <- SID0_cell_individuals[(SID0_cell_individuals$clusters %in% c(0,4,6) & 
                                                  (SID0_cell_individuals$BEST %in% "SNG-Line_14_NG_11923")), 1] 
SID0_GABA_barcodes_15 <- SID0_cell_individuals[(SID0_cell_individuals$clusters %in% c(0,4,6) & 
                                                  (SID0_cell_individuals$BEST %in% "SNG-Line_15_NG_55047")), 1] 


SID0_GABA_barcodes <- SID0_cell_individuals[SID0_cell_individuals$clusters %in% c(0,4,6), 1]


# write tables
write.table(SID0_NGN2_Glut_barcodes_14, 
            file = "SID0_NGN2_Glut_barcodes_14.txt", quote = F, 
            sep = "\t", row.names = F, col.names = F)
write.table(SID0_NGN2_Glut_barcodes_15, 
            file = "SID0_NGN2_Glut_barcodes_15.txt", quote = F, 
            sep = "\t", row.names = F, col.names = F)
write.table(SID0_GABA_barcodes_14, 
            file = "SID0_GABA_barcodes_14.txt", quote = F, 
            sep = "\t", row.names = F, col.names = F)
write.table(SID0_GABA_barcodes_15, 
            file = "SID0_GABA_barcodes_15.txt", quote = F, 
            sep = "\t", row.names = F, col.names = F)




write.table(SID0_NPC_barcodes, 
            file = "SID0_NPC_barcodes.txt", quote = F, 
            sep = "\t", row.names = F, col.names = F)


save.image(file = "SID0.scRNAseq.scATACseq.umap.RData")
