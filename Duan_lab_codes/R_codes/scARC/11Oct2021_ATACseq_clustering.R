## Chuxuan Li 10/11/2021
# ATACseq merging dataset method 1: use 10x aggregate function
# clustering of ATACseq data


# init
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
library(stringr)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)


set.seed(1886)
# set threads and parallelization
plan("multisession", workers = 2)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 11474836480)

# # Read files, separate the ATACseq data from the .h5 matrix
# g_2_0_raw <- 
#   Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_2_0/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information 
# g_2_0_raw_pks <- 
#   CreateSeuratObject(counts = g_2_0_raw$Peaks,
#                      project = "g_2_0_raw_pks")
# 
# g_2_1_raw <- 
#   Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_2_1/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information
# g_2_1_raw_pks <- 
#   CreateSeuratObject(counts = g_2_1_raw$Peaks,
#                      project = "g_2_1_raw_pks")
# 
# g_2_6_raw <- 
#   Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_2_6/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information
# g_2_6_raw_pks <- 
#   CreateSeuratObject(counts = g_2_6_raw$Peaks,
#                      project = "g_2_6_raw_pks")
# 
# 
# g_8_0_raw <- 
#   Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_8_0/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information
# g_8_0_raw_pks <- 
#   CreateSeuratObject(counts = g_8_0_raw$Peaks,
#                      project = "g_8_0_raw_pks")
# 
# g_8_1_raw <- 
#   Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_8_1/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information
# g_8_1_raw_pks <- 
#   CreateSeuratObject(counts = g_8_1_raw$Peaks,
#                      project = "g_8_1_raw_pks")
# 
# g_8_6_raw <- 
#   Read10X_h5(filename = "../hg38_Rnor6_mixed/libraries_8_6/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information
# g_8_6_raw_pks <- 
#   CreateSeuratObject(counts = g_8_6_raw$Peaks,
#                      project = "g_8_6_raw_pks")

# read in the combined object obtained by the -aggr method from 10x
aggr_raw <- Read10X_h5(filename = "/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/filtered_feature_bc_matrix.h5")
# load ATACseq information
aggr_raw_pks <- 
  CreateSeuratObject(counts = aggr_raw$Peaks,
                     project = "aggregated_raw_pks")

chrom_assay <- CreateChromatinAssay(
  counts = aggr_raw$Peaks,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/atac_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

aggr_signac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

aggr_signac$group_time.ident <- str_sub(aggr_signac@assays$peaks@counts@Dimnames[[2]], 
                                        start = -1L)

aggr_signac$time.ident <- NA
aggr_signac$time.ident[aggr_signac$group_time.ident %in% c("1", "4")] <- "0hr"
aggr_signac$time.ident[aggr_signac$group_time.ident %in% c("2", "5")] <- "1hr"  
aggr_signac$time.ident[aggr_signac$group_time.ident %in% c("3", "6")] <- "6hr"  
# check if anything is unlabeled
sum(is.na(aggr_signac$time.ident))

# annot_aggr_signac <-
#   GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86,
#                       standard.chromosomes = T)
# make a txdb from Gencode database
# txdb.v32 <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# 
# exon_txdb_v32 <- transcripts(txdb.v32)
Annotation(aggr_signac) <- annot_aggr_signac_ucsc
#seqlevelsStyle(annot_aggr_signac) <- 'UCSC' # cannot  use FTP service


# QC #####
aggr_signac <- NucleosomeSignal(aggr_signac)
aggr_signac <- TSSEnrichment(aggr_signac, 
                             fast = F # if True, cannot plot the accessibility profile at the TSS
)

# look at TSS enrichment plot
aggr_signac$high.tss <- ifelse(aggr_signac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(aggr_signac, group.by = 'high.tss') + NoLegend()
TSSPlot(aggr_signac) + NoLegend()

aggr_signac$nucleosome_group <- ifelse(aggr_signac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = aggr_signac, group.by = 'nucleosome_group')

# peak_region_fragments = Col sum of the peak matrix
# passed_filters = total readsi (not just reads in peaks) per cell
# nFeatures = number of unique reads in one cell
# nCount = number of all reads in one cell
# calculate the percent of reads that are in peaks
peaks_enrichment <- TSSEnrichment(aggr_signac,
                                  tss.positions = aggr_signac@assays$peaks@ranges,
                                  fast = F)
aggr_signac$peak_region_fragments <- peaks_enrichment$TSS.enrichment

aggr_signac$pct_reads_in_peaks <- aggr_signac$peak_region_fragments / 
  aggr_signac$nCount_peaks * 100


# get blacklisted regions from existing database
blacklisted_regions <- 
  read_delim("~/Data/Databases/Genomes/hg38/hg38_blacklisted_regions.bed", 
             delim = "\t", escape_double = FALSE, 
             col_names = FALSE, trim_ws = TRUE)

# turn bed file into Granges
blacklisted_granges <- GRanges(seqnames = blacklisted_regions$X1,
                               ranges = IRanges(start = blacklisted_regions$X2,
                                                end = blacklisted_regions$X3))

# get number of fragments in the blacklisted regions using TSSenrichment tss.positions argument
blacklisted_enrichment <- TSSEnrichment(aggr_signac, 
                                        tss.positions = blacklisted_granges,
                                        fast = F) # if True, cannot plot the accessibility profile at the TSS)
aggr_signac$blacklist_region_fragments <- blacklisted_enrichment$TSS.enrichment

# calculate the blacklisted fragments in all reads that are in peaks
aggr_signac$blacklist_ratio <- 
  aggr_signac$blacklist_region_fragments / aggr_signac$peak_region_fragments


# 
# aggr_signac$pct_reads_in_peaks <- 
#   colSums(aggr_signac@assays$peaks@counts) / 
#   aggr_signac$nCount_peaks * 100

#save("aggr_signac", file = "aggregated_signac_object.RData")

VlnPlot(
  object = aggr_signac,
  features = c('nFeature_peaks','peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.05,
  ncol = 5
)


# using filtered object to look at RNAseq projection #####
# remove outliers
load("/home/cli/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/aggregated_signac_object.RData")

aggr_signac_filtered <- subset(
  x = aggr_signac,
  subset = nFeature_peaks > 1000 & 
    nFeature_peaks < 20000 &
    #pct_reads_in_peaks > 15 &
    #blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
# calculate some distribution
median(aggr_signac$nFeature_peaks) # 2504
hist(aggr_signac$nFeature_peaks, 
     breaks = 1000,
     xlim = c(0, 10000))
hist(aggr_signac$nucleosome_signal,
     breaks = 100)
hist(aggr_signac$TSS.enrichment,
     breaks = 1000,
     xlim = c(0, 10))
hist(aggr_signac$blacklist_ratio, # peak at approx 0.2
     breaks = 1000,
     xlim = c(0, 2))

# linear reduction
aggr_signac_filtered <- RunTFIDF(aggr_signac_filtered)
aggr_signac_filtered <- FindTopFeatures(aggr_signac_filtered, min.cutoff = 'q0')
aggr_signac_filtered <- RunSVD(aggr_signac_filtered)

# check which components are associated with technical (nonbiological) variance (sequencing depth)
DepthCor(aggr_signac_filtered)

# clustering and downstream analysis
aggr_signac_filtered <- RunUMAP(object = aggr_signac_filtered, 
                                reduction = 'lsi', 
                                dims = 2:30, 
                                seed.use = 1889
                                ) # remove the first component
aggr_signac_filtered <- FindNeighbors(object = aggr_signac_filtered, 
                                      reduction = 'lsi', 
                                      dims = 2:30) # remove the first component
aggr_signac_filtered <- FindClusters(object = aggr_signac_filtered, 
                                     verbose = FALSE, 
                                     algorithm = 3, 
                                     resolution = 0.5,
                                     random.seed = 1889)

DimPlot(object = aggr_signac_filtered, label = TRUE) + 
  NoLegend() +
  ggtitle("clustering of ATACseq data using aggregated dataset")


save("aggr_signac_filtered", file = "aggr_signac_filtered.RData")


# load the RNAseq object
integrated <- readRDS("/home/cli/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/integrated_by_Seurat.rds")

# get the barcodes for each cluster from the RNAseq object
#integrated <- integrated_by_Seurat

aggr_signac_filtered$RNA.cluster.ident <- NA
aggr_signac_filtered$trimmed.barcodes <- str_sub(aggr_signac_filtered@assays$peaks@data@Dimnames[[2]],
                                                 end = -3L)

for (i in 1:length(levels(integrated$seurat_clusters))){
  barcodes <- str_sub(integrated@assays$RNA@counts@Dimnames[[2]][integrated$seurat_clusters %in% 
                                                           as.character(i)], end = -7L)
  
  aggr_signac_filtered$RNA.cluster.ident[aggr_signac_filtered$trimmed.barcodes 
                                         %in% barcodes] <- as.character(i) 
}

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = aggr_signac_filtered, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.cluster.ident") + 
  ggtitle("clustering of ATACseq data with RNAseq projection")


# project ATACseq onto RNAseq
integrated$ATAC.cluster.ident <- NA
integrated$trimmed.barcodes <- str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -7L)

for (i in 1:length(levels(aggr_signac_filtered$seurat_clusters))){
  barcodes <- str_sub(aggr_signac_filtered@assays$peaks@counts@Dimnames[[2]][aggr_signac_filtered$seurat_clusters %in% 
                                                                   as.character(i)], end = -3L)
  
  integrated$ATAC.cluster.ident[integrated$trimmed.barcodes 
                                         %in% barcodes] <- as.character(i) 
}

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = integrated, 
        label = TRUE, 
        repel = T,
        group.by = "ATAC.cluster.ident") + 
  ggtitle("clustering of RNAseq data with ATACseq projection")


# add an attribute onto the RNAseq data from ATACseq nCount and nFeature
integrated$ATAC.nCount <- NA
integrated$ATAC.nFeature <- NA
integrated$ATAC.nCount[integrated$trimmed.barcodes %in% aggr_signac_filtered$trimmed.barcodes] <- 
  aggr_signac_filtered$nCount_peaks[integrated$trimmed.barcodes %in% aggr_signac_filtered$trimmed.barcodes]
integrated$ATAC.nFeature[integrated$trimmed.barcodes %in% aggr_signac_filtered$trimmed.barcodes] <- 
  aggr_signac_filtered$nFeature_peaks[integrated$trimmed.barcodes %in% aggr_signac_filtered$trimmed.barcodes]


# determine the cutoff
hist(x = integrated$ATAC.nCount,
     breaks = 1000)
ncount_feat <- FeaturePlot(integrated, 
                           features = "ATAC.nCount",
                           max.cutoff = 15000,
                           cols = c("grey", "red3")
                          )
ncount_feat[[1]]$layers[[1]]$aes_params$alpha = .3
ncount_feat

hist(x = integrated$ATAC.nFeature,
     breaks = 1000)
nfeature_feat <- FeaturePlot(integrated, 
            features = "ATAC.nFeature",
            max.cutoff = 10000,
            cols = c("grey", "red3"))
nfeature_feat[[1]]$layers[[1]]$aes_params$alpha = .3
nfeature_feat


# check how many cells in RNAseq are not labeled by ATACseq
sum(is.na(integrated$ATAC.nFeature))/length(integrated$ATAC.nFeature)
# check how many cells in ATACseq are not labeled by RNAseq
sum(is.na(aggr_signac_filtered$RNA.cluster.ident))/length(aggr_signac_filtered$RNA.cluster.ident)

# without QC for ATACseq
sum(str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -7L) %in% 
         str_sub(aggr_signac@assays$peaks@data@Dimnames[[2]],end = -3L)) / 
  length(integrated@assays$RNA@data@Dimnames[[2]])

sum(str_sub(aggr_signac@assays$peaks@data@Dimnames[[2]], end = -3L) %in% 
          str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -7L))/ 
  length(aggr_signac@assays$peaks@data@Dimnames[[2]])



# using fewer filters filtered object to compare with filtered results #####
aggr_signac_unfiltered <- subset(x = aggr_signac,
                                 subset = nFeature_peaks < 20000 &
                                   #pct_reads_in_peaks > 15 &
                                   #blacklist_ratio < 0.05 &
                                   nucleosome_signal < 4 &
                                   TSS.enrichment > 2
)

sum(aggr_signac$nFeature_peaks < 20000 &
  aggr_signac$nucleosome_signal < 4 &
  aggr_signac$TSS.enrichment > 2)

sum(aggr_signac$nFeature_peaks < 20000 &
         aggr_signac$nFeature_peaks > 1000 & 
         aggr_signac$nucleosome_signal < 4 &
         aggr_signac$TSS.enrichment > 2)

# linear reduction
aggr_signac_unfiltered <- RunTFIDF(aggr_signac_unfiltered)
aggr_signac_unfiltered <- FindTopFeatures(aggr_signac_unfiltered, min.cutoff = 'q0')
aggr_signac_unfiltered <- RunSVD(aggr_signac_unfiltered)

# check which components are associated with technical (nonbiological) variance (sequencing depth)
DepthCor(aggr_signac_unfiltered)

# clustering and downstream analysis
aggr_signac_unfiltered <- RunUMAP(object = aggr_signac_unfiltered, 
                                reduction = 'lsi', 
                                dims = 2:30, 
                                seed.use = 99
) # remove the first component
aggr_signac_unfiltered <- FindNeighbors(object = aggr_signac_unfiltered, 
                                      reduction = 'lsi', 
                                      dims = 2:30) # remove the first component
aggr_signac_unfiltered <- FindClusters(object = aggr_signac_unfiltered, 
                                     verbose = FALSE, 
                                     algorithm = 3, 
                                     resolution = 0.5,
                                     random.seed = 99)

DimPlot(object = aggr_signac_unfiltered, label = TRUE) + 
  NoLegend() +
  ggtitle("clustering of ATACseq data using aggregated dataset\n (TSS enrichment > 2)")


save("aggr_signac_unfiltered", file = "aggr_signac_unfiltered.RData")


# get the barcodes for each cluster from the RNAseq object
integrated <- readRDS("/home/cli/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/integrated_by_Seurat.rds")

aggr_signac_unfiltered$RNA.cluster.ident <- NA
aggr_signac_unfiltered$trimmed.barcodes <- str_sub(aggr_signac_unfiltered@assays$peaks@data@Dimnames[[2]],
                                                 end = -3L)

for (i in 1:length(levels(integrated_renamed_1$seurat_clusters))){
  barcodes <- str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]][integrated_renamed_1$seurat_clusters %in% 
                                                                   as.character(i)], end = -7L)
  
  aggr_signac_unfiltered$RNA.cluster.ident[aggr_signac_unfiltered$trimmed.barcodes 
                                         %in% barcodes] <- as.character(i) 
}

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = aggr_signac_unfiltered, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.cluster.ident") + 
  ggtitle("clustering of ATACseq data with RNAseq projection\n (TSS enrichment > 2)")


# project ATACseq onto RNAseq
integrated_renamed_1$ATAC.cluster.ident <- NA
integrated_renamed_1$trimmed.barcodes <- str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]], end = -7L)

for (i in 1:length(levels(aggr_signac_unfiltered$seurat_clusters))){
  barcodes <- str_sub(aggr_signac_unfiltered@assays$peaks@counts@Dimnames[[2]][aggr_signac_unfiltered$seurat_clusters %in% 
                                                                               as.character(i)], end = -3L)
  
  integrated_renamed_1$ATAC.cluster.ident[integrated_renamed_1$trimmed.barcodes 
                                %in% barcodes] <- as.character(i) 
}

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = integrated_renamed_1, 
        label = TRUE, 
        repel = T,
        group.by = "ATAC.cluster.ident") + 
  ggtitle("clustering of RNAseq data with ATACseq projection\n (TSS enrichment > 2)")


# add an attribute onto the RNAseq data from ATACseq nCount and nFeature
integrated_renamed_1$ATAC.nCount <- NA
integrated_renamed_1$ATAC.nFeature <- NA
integrated_renamed_1$ATAC.nCount[integrated_renamed_1$trimmed.barcodes %in% aggr_signac_unfiltered$trimmed.barcodes] <- 
  aggr_signac_unfiltered$nCount_peaks[integrated_renamed_1$trimmed.barcodes %in% aggr_signac_unfiltered$trimmed.barcodes]
integrated_renamed_1$ATAC.nFeature[integrated_renamed_1$trimmed.barcodes %in% aggr_signac_unfiltered$trimmed.barcodes] <- 
  aggr_signac_unfiltered$nFeature_peaks[integrated_renamed_1$trimmed.barcodes %in% aggr_signac_unfiltered$trimmed.barcodes]


# determine the cutoff
hist(x = integrated_renamed_1$ATAC.nCount,
     breaks = 1000,
     xlim = c(0, 50000))
ncount_feat <- FeaturePlot(integrated_renamed_1, 
                           features = "ATAC.nCount",
                           max.cutoff = 25000,
                           cols = c("grey", "red3")
)
ncount_feat[[1]]$layers[[1]]$aes_params$alpha = .3
ncount_feat

hist(x = integrated_renamed_1$ATAC.nFeature,
     breaks = 1000,
     xlim = c(0, 20000))
nfeature_feat <- FeaturePlot(integrated_renamed_1, 
                             features = "ATAC.nFeature",
                             max.cutoff = 16000,
                             cols = c("grey", "red3"))
nfeature_feat[[1]]$layers[[1]]$aes_params$alpha = .3
nfeature_feat


# check how many cells in RNAseq are not labeled by ATACseq
sum(is.na(integrated_renamed_1$ATAC.nFeature))/length(integrated_renamed_1$ATAC.nFeature)
# check how many cells in ATACseq are not labeled by RNAseq
sum(is.na(aggr_signac_unfiltered$RNA.cluster.ident))/length(aggr_signac_unfiltered$RNA.cluster.ident)

# without QC for ATACseq
sum(str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]], end = -7L) %in% 
      str_sub(aggr_signac@assays$peaks@data@Dimnames[[2]],end = -3L)) / 
  length(integrated_renamed_1@assays$RNA@data@Dimnames[[2]])

sum(str_sub(aggr_signac@assays$peaks@data@Dimnames[[2]], end = -3L) %in% 
      str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]], end = -7L))/ 
  length(aggr_signac@assays$peaks@data@Dimnames[[2]])


# check duplicate barcodes 
length(aggr_signac_unfiltered$trimmed.barcodes) #total number of barcodes (ATAC): 80970
length(unique(aggr_signac_unfiltered$trimmed.barcodes)) # unique num of barcodes (ATAC): 77183
length(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]]) #total num of barcodes (RNA): 83093
length(unique(integrated_renamed_1$trimmed.barcodes)) # unique num of barcodes (RNA): 79159
77183/80970 # unique ATAC barcodes to all ATAC barcodes ratio: 0.9532296
79159/83093 # unique RNA barcodes to all RNA barcodes ratio: 0.9526555

sum(is.na(aggr_signac_unfiltered$RNA.cluster.ident))

# load labeled RNAseq objects
load("/home/cli/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/integrated_labeled_updated.RData")

# create another attribute to project RNAseq cell type ident onto the ATACseq data
aggr_signac_unfiltered$RNA.cell.type.ident <- NA
# find all barcodes corresponding to each cell type

for (i in levels(integrated_renamed_1@active.ident)){
  barcodes <- str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]][integrated_renamed_1@active.ident %in% i], 
                      end = -7L)
  
  aggr_signac_unfiltered$RNA.cell.type.ident[aggr_signac_unfiltered$trimmed.barcodes %in% barcodes] <- i
}

save("aggr_signac_unfiltered", file = "aggr_ATACseq_with_all_labels.RData")

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = aggr_signac_unfiltered, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.cell.type.ident") + 
  ggtitle("clustering of ATACseq data with cell type projection")


# plot peak tracks as an example
CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr11-22334000-22340000", # SLC17A6
  extend.upstream = 1,
  extend.downstream = 1
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr19-49427411-49431451", # SLC17A7
  extend.upstream = 10000,
  extend.downstream = 1
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr10-26204609-26224140" # GAD2
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr2-190871076-190883628" # GLS
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr7-101805750-101822128" # CUX1
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr12-111026276-111037993" # CUX2
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr14-28759074-28770857", # FOXG1
  extend.upstream = 10000
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr3-181703353-181714656" # SOX2
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr1-156675466-156688433" # NES
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr14-75269196-75280682" # FOS
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr5-138456812-138468286" # EGR1
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr11-27719113-27734038" # BDNF
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr6-6005375-6014978" # NRN1
)

CoveragePlot(
  object = aggr_signac_unfiltered,
  region = "chr7-101163942-101173124" # VGF
)

# separate time points, look at response markers for each time point
aggr_signac_0hr <- subset(aggr_signac_unfiltered, subset = time.ident == "0hr")
aggr_signac_1hr <- subset(aggr_signac_unfiltered, subset = time.ident == "1hr")
aggr_signac_6hr <- subset(aggr_signac_unfiltered, subset = time.ident == "6hr")

CoveragePlot(
  object = aggr_signac_0hr,
  region = "chr14-75269196-75280682" # FOS
) + ggtitle("0hr")

CoveragePlot(
  object = aggr_signac_1hr,
  region = "chr14-75269196-75280682" # FOS
) + ggtitle("1hr")

CoveragePlot(
  object = aggr_signac_6hr,
  region = "chr14-75269196-75280682" # FOS
) + ggtitle("6hr")



CoveragePlot(
  object = aggr_signac_0hr,
  region = "chr5-138456812-138468286" # EGR1
) + ggtitle("0hr")

CoveragePlot(
  object = aggr_signac_1hr,
  region = "chr5-138456812-138468286" # EGR1
) + ggtitle("1hr")

CoveragePlot(
  object = aggr_signac_6hr,
  region = "chr5-138456812-138468286" # EGR1
) + ggtitle("6hr")


CoveragePlot(
  object = aggr_signac_0hr,
  region = "chr11-27719113-27734038" # BDNF
) + ggtitle("0hr")

CoveragePlot(
  object = aggr_signac_1hr,
  region = "chr11-27719113-27734038" # BDNF
) + ggtitle("1hr")

CoveragePlot(
  object = aggr_signac_6hr,
  region = "chr11-27719113-27734038" # BDNF
) + ggtitle("6hr")


CoveragePlot(
  object = aggr_signac_0hr,
  region = "chr6-6005375-6014978" # NRN1
) + ggtitle("0hr")

CoveragePlot(
  object = aggr_signac_1hr,
  region = "chr6-6005375-6014978" # NRN1
) + ggtitle("1hr")

CoveragePlot(
  object = aggr_signac_6hr,
  region = "chr6-6005375-6014978" # NRN1
) + ggtitle("6hr")


CoveragePlot(
  object = aggr_signac_0hr,
  region = "chr7-101163942-101173124" # VGF
) + ggtitle("0hr")

CoveragePlot(
  object = aggr_signac_1hr,
  region = "chr7-101163942-101173124" # VGF
) + ggtitle("1hr")

CoveragePlot(
  object = aggr_signac_6hr,
  region = "chr7-101163942-101173124" # VGF
) + ggtitle("6hr")

# link peaks to genes 
# DefaultAssay(aggr_signac_unfiltered) <- "peaks"
# aggr_hg38_matched <- aggr_signac_unfiltered
# # aggr_hg38_matched@assays$peaks <- aggr_hg38_matched@assays$peaks[!(aggr_signac_unfiltered@assays$peaks@ranges@seqnames@values %in% 
# #                                                                    c("KI270728.1", "KI270727.1", 
# #                                                                      "GL000205.2", "GL000195.1",
# #                                                                      "GL000219.1", "KI270734.1", 
# #                                                                      "KI270731.1", "KI270711.1", "KI270713.1"))]
# # aggr_hg38_matched@assays$peaks <- aggr_hg38_matched@assays$peaks[aggr_signac_unfiltered@assays$peaks@counts@Dimnames[[1]] %in% 
# #                                                                    (str_sub(aggr_signac_unfiltered@assays$peaks@counts@Dimnames[[1]],
# #                                                                             end = 4L) %in% "chr")]
# aggr_hg38_matched@assays$peaks <- aggr_hg38_matched@assays$peaks[aggr_signac_unfiltered@assays$peaks@counts@Dimnames[[1]] %in% 
#                                                                    (str_sub(aggr_signac_unfiltered@assays$peaks@counts@Dimnames[[1]],
#                                                                             end = 4L) %in% "chr")]
# 
# aggr_hg38_matched@assays$peaks@counts[ ,(aggr_signac_unfiltered@assays$peaks@ranges@seqnames@values %in% 
#                                   c("KI270728.1", "KI270727.1",
#                                     "GL000205.2", "GL000195.1",
#                                     "GL000219.1", "KI270734.1",
#                                     "KI270731.1", "KI270711.1", "KI270713.1"))] <- NULL
# aggr_hg38_matched <- RegionStats(aggr_hg38_matched,
#                                  genome = BSgenome.Hsapiens.UCSC.hg38)

