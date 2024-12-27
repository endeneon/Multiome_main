# Chuxuan Li 12/10/2021
# comparing 1) using MACS2 on the entire ATACseq dataset, 2) using MACS2 on pre-identified
#cell types (merging the time points)

# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(readr)
library(ggplot2)
library(future)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)

set.seed(1886)
# set threads and parallelization
plan("multisession", workers = 2)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 11474836480)

# read in the combined object obtained by the -aggr method from 10x
aggr_raw <- Read10X_h5(filename = "/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/filtered_feature_bc_matrix.h5")

aggr <- CreateSeuratObject(
  counts = aggr_raw$`Gene Expression`,
  assay = "RNA"
)

peaks <- aggr_raw$Peaks

# choose the peaks that only start with chr
peaks_subsetted <- peaks[str_sub(peaks@Dimnames[[1]], end = 3L) %in% "chr", ]

# create chromatin assay without filtering
chrom_assay <- CreateChromatinAssay(
  counts = peaks_subsetted,
  sep = c(":", "-"),
  #genome = 'hg38',
  fragments = '/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/atac_fragments.tsv.gz',
  annotation = annot_aggr_signac_ucsc
  # min.cells = 10,
  # min.features = 200
)

aggr[["ATAC"]] <- chrom_assay

aggr$time.group.ident <- str_sub(aggr@assays$ATAC@counts@Dimnames[[2]], 
                                 start = -1L)

aggr$time.ident <- NA
aggr$time.ident[aggr$time.group.ident %in% c("1", "4")] <- "0hr"
aggr$time.ident[aggr$time.group.ident %in% c("2", "5")] <- "1hr"  
aggr$time.ident[aggr$time.group.ident %in% c("3", "6")] <- "6hr"  
# check if anything is unlabeled
sum(is.na(aggr$time.ident))


# QC #####
DefaultAssay(aggr) <- "ATAC"
aggr <- NucleosomeSignal(aggr)
aggr <- TSSEnrichment(aggr, 
                      fast = F # if True, cannot plot the accessibility profile at the TSS
)

# using fewer filters filtered object to compare with filtered results
aggr_filtered <- subset(x = aggr,
                        subset = nFeature_ATAC < 20000 &
                          #pct_reads_in_peaks > 15 &
                          #blacklist_ratio < 0.05 &
                          nucleosome_signal < 4 &
                          TSS.enrichment > 2)


# call peaks using MACS2
peaks <- CallPeaks(aggr_filtered, macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(aggr_filtered),
  features = peaks,
  cells = colnames(aggr_filtered)
)
fragpath <- "/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/atac_fragments.tsv.gz"
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/annot_aggr_signac_EnsDb_UCSC.RData")

# create a new assay using the MACS2 peak set and add it to the Seurat object
aggr_filtered[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annot_aggr_signac_ucsc
)

aggr_filtered$time.group.ident <- str_sub(aggr_filtered@assays$ATAC@counts@Dimnames[[2]], 
                                 start = -1L)


# dimensional reduction ####
DefaultAssay(aggr_filtered) <- "peaks"

aggr_filtered <- RunTFIDF(aggr_filtered)
aggr_filtered <- FindTopFeatures(aggr_filtered, min.cutoff = 'q0')
aggr_filtered <- RunSVD(aggr_filtered)

save("aggr_filtered", file = "ATAC_RNA_combined_10x_aggregated_data_filtered.RData")

DefaultAssay(aggr_filtered) <- "RNA"
aggr_filtered <- PercentageFeatureSet(aggr_filtered, pattern = "^MT-", col.name = "percent.mt")

aggr_filtered <- SCTransform(aggr_filtered,
                             vars.to.regress = "percent.mt",
                             method = "glmGamPoi",
                             variable.features.n = 9000,
                             verbose = T,
                             seed.use = 42)
aggr_filtered <- RunPCA(aggr_filtered)

# build a joint neighbor graph using both assays
aggr_filtered <- FindMultiModalNeighbors(
  object = aggr_filtered,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
aggr_filtered <- RunUMAP(
  object = aggr_filtered,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)


aggr_filtered <- FindNeighbors(object = aggr_filtered, 
                               reduction = 'lsi', 
                               dims = 2:30) # remove the first component
aggr_filtered <- FindClusters(object = aggr_filtered, 
                              verbose = FALSE, 
                              graph.name = "peaks_snn",
                              algorithm = 3, 
                              resolution = 0.5,
                              random.seed = 99)

# integrate with RNAseq dataset ####
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/integrated_with_all_labels.RData")
# get the barcodes for each cluster from the RNAseq object
DefaultAssay(aggr_filtered) <- "peaks"
aggr_filtered$RNA.cluster.ident <- NA
aggr_filtered$trimmed.barcodes <- str_sub(aggr_filtered@assays$ATAC@data@Dimnames[[2]],
                                          end = -3L)

for (i in 1:length(levels(integrated_renamed_1$seurat_clusters))){
  barcodes <- str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]][integrated_renamed_1$seurat_clusters %in% 
                                                                             as.character(i)], end = -7L)
  
  aggr_filtered$RNA.cluster.ident[aggr_filtered$trimmed.barcodes 
                                  %in% barcodes] <- as.character(i) 
}

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = aggr_filtered, 
        label = TRUE, 
        repel = T,
        reduction = "umap",
        group.by = "RNA.cluster.ident") + 
  ggtitle("clustering of ATACseq data with RNAseq projection")

