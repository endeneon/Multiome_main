# extract barcodes from anchored SID0/1/6
# for SNP and peak calling
# Use the large RNASeq.integrated object

options(future.fork.enable = T) # force enable forks

# setup multithread
plan("multicore", workers = 8) # force multicore

RNASeq_SID0 <- RNASeq.integrated[, RNASeq.integrated$orig.ident %in% "SID0"]
RNASeq_SID1 <- RNASeq.integrated[, RNASeq.integrated$orig.ident %in% "SID1"]
RNASeq_SID6 <- RNASeq.integrated[, RNASeq.integrated$orig.ident %in% "SID6"]

RNASeq_SID <- RNASeq_SID1
# read in scATAC-seq data, change for each time point
# SID.ATACseq.data <- Read10X_h5("../scATAC/SID0_scATAC/outs/filtered_peak_bc_matrix.h5")
SID.ATACseq.data <- Read10X_h5("../scATAC/SID1_scATAC/outs/filtered_peak_bc_matrix.h5")
# SID0_demux_RNA_index <- read_table2("~/Galaxy_Drive/May15_Duan_017/demultiplex/SID0_scRNA_099_ASE/SID0_scRNA_099_ASE_cell_line_index.txt")
SID_demux_ATAC_index <- read_table2("~/Galaxy_Drive/May15_Duan_017/demultiplex/SID1_ATAC/SID1_ATAC.best")
SID_demux_ATAC_index <- SID_demux_ATAC_index[, 1:6]
meta <- read.table("../scATAC/SID1_scATAC/outs/singlecell.csv", # 
                   sep = ",", header = T, row.names = 1, stringsAsFactors = F)
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

# DefaultAssay(RNASeq_SID0)
SID_cleaned_ATACseq <- SID.ATACseq.data[, colnames(SID.ATACseq.data) %in% 
                                            SID_demux_ATAC_index$BARCODE]
SID_cleaned_ATACseq <- SID_cleaned_ATACseq[, order(colnames(SID_cleaned_ATACseq))]
# make a barcode/cell identity data frame
SID_demux_ATAC_index_filtered <- SID_demux_ATAC_index[SID_demux_ATAC_index$BARCODE %in%
                                                         colnames(SID.ATACseq.data), ]
SID_demux_ATAC_index_filtered <- SID_demux_ATAC_index_filtered[order(SID_demux_ATAC_index_filtered$BARCODE), ]

# Try to load the cleaned count matrix since the raw matrix is too large
peak_count_matrix <- CreateGeneActivityMatrix(
  peak.matrix = SID_cleaned_ATACseq,
  annotation.file = "~/2TB/Databases/hg38_rn6_hybrid/gencode.v30.Rn.6.gtf",
  seq.levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                 "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                 "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                 "chr19", "chr20", "chr21", "chr22" #human
                 , c(1:20) #rat
  ), 
  upstream = 2000,
  keep.sparse = T,
  include.body = F, verbose = T)

# Set up the ATAC object
SID.ATAC.object <- CreateSeuratObject(counts = SID_cleaned_ATACseq, 
                                       assay = "ATAC", 
                                       project = "10x_ATAC")
SID.ATAC.object[["ACTIVITY"]] <- CreateAssayObject(counts = peak_count_matrix)

meta <- meta[order(colnames(SID_cleaned_ATACseq)), ] # only take metadata of identified cells



# preprocess data
DefaultAssay(SID.ATAC.object) <- "ACTIVITY"
SID.ATAC.object <- FindVariableFeatures(SID.ATAC.object, nfeatures = 10000)
SID.ATAC.object <- NormalizeData(SID.ATAC.object)
SID.ATAC.object <- ScaleData(SID.ATAC.object)

# use peaks that have at least 100 reads per cell
DefaultAssay(SID.ATAC.object) <- "ATAC"
VariableFeatures(SID.ATAC.object) <- names(which(Matrix::rowSums(SID.ATAC.object) > 100))

SID.ATAC.object <- RunLSI(SID.ATAC.object, n = 50, scale.max = NULL)
SID.ATAC.object <- RunUMAP(SID.ATAC.object, reduction = "lsi", dims = 1:50)

# # load corresponding RNASeq data
# SID.RNA.object <- readRDS(file = "SID0.RNASeq.RData")
# SID0.RNA.object$tech <- "rna"

# make a combined map
p1 <- DimPlot(SID.ATAC.object, reduction = "umap") +
  NoLegend() +
  ggtitle("SID0 scATAC-seq")
p2 <- DimPlot(RNASeq_SID, 
              label = T, repel = T) +
  NoLegend() +
  ggtitle("SID0 scRNA-seq")
p1 + p2

DimPlot(SID.ATAC.object, 
        reduction = "umap", 
        label = T, repel = T) +
  NoLegend() +
  ggtitle("SID0 scATAC-seq")

# identify corresponding anchors and transfer

RNASeq_SID$tech <- "rna"
DefaultAssay(RNASeq_SID) <- "SCT"
DefaultAssay(SID.ATAC.object) <- "ACTIVITY"
RNASeq_SID <- FindVariableFeatures(RNASeq_SID,
                                   nfeatures = 10000)
SID.transfer.anchors <- FindTransferAnchors(reference = RNASeq_SID, 
                                             query = SID.ATAC.object, 
                                             features = VariableFeatures(object = RNASeq_SID),
                                             reference.assay = "SCT",
                                             query.assay = "ACTIVITY", 
                                             reduction = "cca")
# predict cell types and check efficiency
SID.celltype.predictions <- TransferData(anchorset = SID.transfer.anchors, 
                                          refdata = RNASeq_SID$seurat_clusters, 
                                          weight.reduction = SID.ATAC.object[["lsi"]])
SIDD.ATAC.object <- AddMetaData(SID.ATAC.object, 
                                metadata = SID.celltype.predictions)
hist(SID.ATAC.object$prediction.score.max) # check prediction score
table(SID.ATAC.object$prediction.score.max > 0.5)

# subset ATAC-seq cells and only keep score > 0.5
SID.ATAC.filtered <- subset(SID.ATAC.object, 
                             subset = prediction.score.max > 0.5)
# match group colour
SID.ATAC.filtered$predicted.id <- factor(SID.ATAC.filtered$predicted.id, 
                                          levels = levels(RNASeq_SID))

# make plot
p1 <- DimPlot(SID.ATAC.filtered, group.by = "predicted.id", label = T, repel = T) +
  ggtitle("SID scATAC-seq cells") +
  scale_colour_hue(drop = F) +
  NoLegend()
p2 <- DimPlot(RNASeq_SID, 
              label = T, repel = T) +
  NoLegend() +
  ggtitle("SID scRNA-seq cells")
p1 + p2


DimPlot(SID.ATAC.filtered, 
        group.by = "predicted.id", 
        label = T, repel = T) +
  ggtitle("SID scATAC-seq cells") +
  scale_colour_hue(drop = F) +
  NoLegend()

#### plot cell line identities ######
SID_umap <- FetchData(SID.ATAC.filtered, vars = "UMAP_1")
# SID0_umap$UMAP_1 <- FetchData(SID0.ATAC.object, vars = "UMAP_1")
SID_umap$UMAP_2 <- FetchData(SID.ATAC.filtered, vars = "UMAP_2")
SID_umap$UMAP_2 <- SID_umap$UMAP_2$UMAP_2
SID_umap$clusters <- FetchData(SI0.ATAC.filtered, vars = "predicted.id")
SID_umap$clusters <- SID_umap$clusters$predicted.id
SI0_umap$barcodes <- rownames(SID.ATAC.filtered@reductions$umap)


SID_cell_individuals <- base::merge(x = SID_umap, 
                                     y = SID_demux_ATAC_index, 
                                     by.x = "barcodes", 
                                     by.y = "BARCODE")


SID_cell_individuals[SID_cell_individuals$clusters %in% c(7), 8] <- "Rat Astrocytes"


ggplot(SID_cell_individuals, aes(x = UMAP_1, y = UMAP_2, 
                                  colour = BEST)) +
  geom_point(size = 0.5) +
  scale_colour_manual(name = "Cell type", 
                      values = c("black", "green4", "darkblue")) +
  ggtitle("Unstimulated, scATAC-seq") +
  theme_classic()


##### separate glut, NPC, and GABA groups
SID_NGN2_Glut_barcodes <- SID_cell_individuals[SID_cell_individuals$clusters %in% c(1,6,9), 1]
SID_GABA_barcodes <- SID_cell_individuals[SID_cell_individuals$clusters %in% c(0,2), 1] 
SID_NPC_barcodes <- SID_cell_individuals[SID_cell_individuals$clusters %in% c(3,4), 1] 

# separate lines
SID0_NGN2_Glut_barcodes <- SID0_cell_individuals[(SID0_cell_individuals$clusters %in% c(1,3,8,9,12) &
                                                    (SID0_cel)), 1]



SID0_GABA_barcodes <- SID0_cell_individuals[SID0_cell_individuals$clusters %in% c(0,4,6), 1]


# write tables
write.table(SID0_NGN2_Glut_barcodes, 
            file = "SID0_NGN2_Glut_barcodes.txt", quote = F, 
            sep = "\t", row.names = F, col.names = F)
write.table(SID0_GABA_barcodes, 
            file = "SID0_GABA_barcodes.txt", quote = F, 
            sep = "\t", row.names = F, col.names = F)
write.table(SID0_NPC_barcodes, 
            file = "SID0_NPC_barcodes.txt", quote = F, 
            sep = "\t", row.names = F, col.names = F)
