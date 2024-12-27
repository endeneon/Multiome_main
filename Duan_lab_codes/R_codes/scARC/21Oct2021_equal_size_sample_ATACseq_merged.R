## Chuxuan Li 10/21/2021
# ATACseq merging dataset method 3: use combined equal sized samples
# no normalization, starts from dimensional reduction

# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
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
library(GenomicRanges)
library(AnnotationDbi)
library(patchwork)

library(stringr)

# library(reticulate)

# conda_list()
# use_condaenv("rstudio")


set.seed(1127)
# set threads and parallelization
plan("multisession", workers = 6)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

# debug
# 
# aggr_raw <- Read10X_h5(filename = "/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/filtered_feature_bc_matrix.h5")
# # load ATACseq information
# aggr_raw_pks <- 
#   CreateSeuratObject(counts = aggr_raw$Peaks,
#                      project = "aggregated_raw_pks")
# 
# chrom_assay <- CreateChromatinAssay(
#   counts = aggr_raw$Peaks,
#   sep = c(":", "-"),
#   #genome = 'hg38',
#   fragments = '/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/atac_fragments.tsv.gz',
#   min.cells = 10,
#   min.features = 200
# )
# 
# norm_signac <- CreateSeuratObject(
#   counts = chrom_assay,
#   assay = "peaks")
# 
# m3 <- str_sub(unique(gsub("(^.*):.*", "\\1", aggr_raw_pks@assays$RNA@counts@Dimnames[[1]])), 1L, 2L)
# GenBankAccn_prefixes<- c("GL",  # alt-scaffold or unlocalized-scaffold
#                           "JH",  # alt-scaffold
#                           "KB",  # alt-scaffold
#                           "KI",  # alt-scaffold or unlocalized-scaffold
#                           "KN",  # fix-patch or novel-patch
#                           "KQ",  # fix-patch or novel-patch
#                           "KV",  # fix-patch or novel-patch
#                           "KZ")  # fix-patch or novel-patch
# m32 <- match(m3, GenBankAccn_prefixes)
# m32


normalized_raw <- 
  Read10X_h5(filename = 
               "/home/cli/NVME/scARC_Duan_018/hg38_mapped_use_Glut_GABA_bed/sample_2_8_aggr_use_Glut_GABA_bed_13Oct2021/outs/filtered_feature_bc_matrix.h5")

chrom_assay <- CreateChromatinAssay(
  counts = normalized_raw$Peaks,
  sep = c(":", "-"),
  #genome = 'hg38',
  fragments = '/home/cli/NVME/scARC_Duan_018/hg38_mapped_use_Glut_GABA_bed/sample_2_8_aggr_use_Glut_GABA_bed_13Oct2021/outs/atac_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

norm_signac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)

norm_signac$group_time.ident <- str_sub(norm_signac@assays$peaks@counts@Dimnames[[2]], 
                                        start = -1L)

norm_signac$time.ident <- NA
norm_signac$time.ident[norm_signac$group_time.ident %in% c("1", "4")] <- "0hr"
norm_signac$time.ident[norm_signac$group_time.ident %in% c("2", "5")] <- "1hr"  
norm_signac$time.ident[norm_signac$group_time.ident %in% c("3", "6")] <- "6hr"  
# check if anything is unlabeled
sum(is.na(norm_signac$time.ident))

annot_norm_signac <-
  GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86,
                      standard.chromosomes = T)
# make a txdb from Gencode database
# txdb.v32 <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# 
# exon_txdb_v32 <- transcripts(txdb.v32)
Annotation(norm_signac) <- annot_norm_signac_ucsc
#seqlevelsStyle(annot_norm_signac) <- 'UCSC' # cannot  use FTP service


# QC #####
norm_signac <- NucleosomeSignal(norm_signac)
norm_signac <- TSSEnrichment(norm_signac, 
                             fast = F # if True, cannot plot the accessibility profile at the TSS
)

# look at TSS enrichment plot
norm_signac$high.tss <- ifelse(norm_signac$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(norm_signac, group.by = 'high.tss') + NoLegend()
TSSPlot(norm_signac) + NoLegend()

norm_signac$nucleosome_group <- ifelse(norm_signac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = norm_signac, 
                  group.by = 'nucleosome_group')

VlnPlot(
  object = norm_signac,
  features = c('nFeature_peaks', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.05,
  ncol = 3
)

# calculate some distribution
median(norm_signac$nFeature_peaks) # 3009
hist(norm_signac$nFeature_peaks, 
     breaks = 1000,
     xlim = c(0, 20000))
hist(norm_signac$nucleosome_signal, # centered at ~0.4
     breaks = 100)
hist(norm_signac$TSS.enrichment, # centered at ~4
     breaks = 1000,
     xlim = c(0, 10))


norm_signac_filtered <- subset(x = norm_signac,
                               subset = nFeature_peaks < 20000 & 
                                 nucleosome_signal < 4 &
                                 TSS.enrichment > 2)

# linear reduction
norm_signac_filtered <- RunTFIDF(norm_signac_filtered)
norm_signac_filtered <- FindTopFeatures(norm_signac_filtered, min.cutoff = 'q0')
norm_signac_filtered <- RunSVD(norm_signac_filtered)

# check which components are associated with technical (nonbiological) variance (sequencing depth)
DepthCor(norm_signac_filtered)

# clustering and downstream analysis
norm_signac_filtered <- RunUMAP(object = norm_signac_filtered, 
                                  reduction = 'lsi', 
                                  dims = 2:30, 
                                  seed.use = 99
) # remove the first component
norm_signac_filtered <- FindNeighbors(object = norm_signac_filtered, 
                                        reduction = 'lsi', 
                                        dims = 2:30) # remove the first component
norm_signac_filtered <- FindClusters(object = norm_signac_filtered, 
                                       verbose = FALSE, 
                                       algorithm = 3, 
                                       resolution = 0.5,
                                       random.seed = 99)

DimPlot(object = norm_signac_filtered, label = TRUE) + 
  NoLegend() +
  ggtitle("clustering of ATACseq data using aggregated dataset\n (TSS enrichment > 2)")


save("norm_signac_filtered", file = "norm_signac_filtered.RData")


# get the barcodes for each cluster from the RNAseq object
integrated <- readRDS("/home/cli/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/integrated_by_Seurat.rds")

norm_signac_filtered$RNA.cluster.ident <- NA
norm_signac_filtered$trimmed.barcodes <- str_sub(norm_signac_filtered@assays$peaks@data@Dimnames[[2]],
                                                   end = -3L)

for (i in 1:length(levels(integrated$seurat_clusters))){
  barcodes <- str_sub(integrated@assays$RNA@counts@Dimnames[[2]][integrated$seurat_clusters %in% 
                                                                   as.character(i)], end = -7L)
  
  norm_signac_filtered$RNA.cluster.ident[norm_signac_filtered$trimmed.barcodes 
                                           %in% barcodes] <- as.character(i) 
}

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = norm_signac_filtered, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.cluster.ident") + 
  ggtitle("clustering of ATACseq data with RNAseq projection\n (TSS enrichment > 2)")


# project ATACseq onto RNAseq
integrated$ATAC.cluster.ident <- NA
integrated$trimmed.barcodes <- str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -7L)

for (i in 1:length(levels(norm_signac_filtered$seurat_clusters))){
  barcodes <- str_sub(norm_signac_filtered@assays$peaks@counts@Dimnames[[2]][norm_signac_filtered$seurat_clusters %in% 
                                                                                 as.character(i)], end = -3L)
  
  integrated$ATAC.cluster.ident[integrated$trimmed.barcodes 
                                %in% barcodes] <- as.character(i) 
}

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = integrated, 
        label = TRUE, 
        repel = T,
        group.by = "ATAC.cluster.ident") + 
  ggtitle("clustering of RNAseq data with ATACseq projection\n (TSS enrichment > 2)")


# add an attribute onto the RNAseq data from ATACseq nCount and nFeature
integrated$ATAC.nCount <- NA
integrated$ATAC.nFeature <- NA
integrated$ATAC.nCount[integrated$trimmed.barcodes %in% norm_signac_filtered$trimmed.barcodes] <- 
  norm_signac_filtered$nCount_peaks[integrated$trimmed.barcodes %in% norm_signac_filtered$trimmed.barcodes]
integrated$ATAC.nFeature[integrated$trimmed.barcodes %in% norm_signac_filtered$trimmed.barcodes] <- 
  norm_signac_filtered$nFeature_peaks[integrated$trimmed.barcodes %in% norm_signac_filtered$trimmed.barcodes]


# determine the cutoff
hist(x = integrated$ATAC.nCount,
     breaks = 1000,
     xlim = c(0, 50000))
ncount_feat <- FeaturePlot(integrated, 
                           features = "ATAC.nCount",
                           max.cutoff = 25000,
                           cols = c("grey", "red3")
)
ncount_feat[[1]]$layers[[1]]$aes_params$alpha = .3
ncount_feat

hist(x = integrated$ATAC.nFeature,
     breaks = 1000,
     xlim = c(0, 20000))
nfeature_feat <- FeaturePlot(integrated, 
                             features = "ATAC.nFeature",
                             max.cutoff = 16000,
                             cols = c("grey", "red3"))
nfeature_feat[[1]]$layers[[1]]$aes_params$alpha = .3
nfeature_feat


# check how many cells in RNAseq are not labeled by ATACseq
sum(is.na(integrated$ATAC.nFeature))/length(integrated$ATAC.nFeature)
# check how many cells in ATACseq are not labeled by RNAseq
sum(is.na(norm_signac_filtered$RNA.cluster.ident))/length(norm_signac_filtered$RNA.cluster.ident)

# without QC for ATACseq
sum(str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -7L) %in% 
      str_sub(norm_signac@assays$peaks@data@Dimnames[[2]],end = -3L)) / 
  length(integrated@assays$RNA@data@Dimnames[[2]])

sum(str_sub(norm_signac@assays$peaks@data@Dimnames[[2]], end = -3L) %in% 
      str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -7L))/ 
  length(norm_signac@assays$peaks@data@Dimnames[[2]])


# check duplicate barcodes 
length(norm_signac_filtered$trimmed.barcodes) #total number of barcodes (ATAC): 80970
length(unique(norm_signac_filtered$trimmed.barcodes)) # unique num of barcodes (ATAC): 77183
length(integrated@assays$RNA@counts@Dimnames[[2]]) #total num of barcodes (RNA): 83093
length(unique(integrated$trimmed.barcodes)) # unique num of barcodes (RNA): 79159
77183/80970 # unique ATAC barcodes to all ATAC barcodes ratio: 0.9532296
79159/83093 # unique RNA barcodes to all RNA barcodes ratio: 0.9526555


# load labeled RNAseq objects
load("/home/cli/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/integrated_labeled.RData")

# create another attribute to project RNAseq cell type ident onto the ATACseq data
norm_signac_filtered$RNA.cell.type.ident <- NA
# find all barcodes corresponding to each cell type

for (i in levels(integrated_renamed@active.ident)){
  barcodes <- str_sub(integrated_renamed@assays$RNA@counts@Dimnames[[2]][integrated_renamed@active.ident %in% i], 
                      end = -7L)
  
  norm_signac_filtered$RNA.cell.type.ident[norm_signac_filtered$trimmed.barcodes %in% barcodes] <- i
}


# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = norm_signac_filtered, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.cell.type.ident") + 
  ggtitle("clustering of ATACseq data with cell type projection")




# plot peak tracks as an example
CoveragePlot(
  object = norm_signac_filtered,
  region = "chr11-22334000-22340000", # SLC17A6
  extend.upstream = 1,
  extend.downstream = 1
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr19-49427411-49431451", # SLC17A7
  extend.upstream = 10000,
  extend.downstream = 1
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr10-26204609-26224140" # GAD2
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr2-190871076-190883628" # GLS
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr7-101805750-101822128" # CUX1
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr12-111026276-111037993" # CUX2
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr14-28759074-28770857", # FOXG1
  extend.upstream = 10000
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr3-181703353-181714656" # SOX2
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr1-156675466-156688433" # NES
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr14-75269196-75280682" # FOS
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr5-138456812-138468286" # EGR1
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr11-27719113-27734038" # BDNF
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr6-6005375-6014978" # NRN1
)

CoveragePlot(
  object = norm_signac_filtered,
  region = "chr7-101163942-101173124" # VGF
)

# separate time points, look at response markers for each time point
norm_signac_0hr <- subset(norm_signac_filtered, subset = time.ident == "0hr")
norm_signac_1hr <- subset(norm_signac_filtered, subset = time.ident == "1hr")
norm_signac_6hr <- subset(norm_signac_filtered, subset = time.ident == "6hr")

CoveragePlot(
  object = norm_signac_0hr,
  region = "chr14-75269196-75280682" # FOS
) + ggtitle("0hr")

CoveragePlot(
  object = norm_signac_1hr,
  region = "chr14-75269196-75280682" # FOS
) + ggtitle("1hr")

CoveragePlot(
  object = norm_signac_6hr,
  region = "chr14-75269196-75280682" # FOS
) + ggtitle("6hr")



CoveragePlot(
  object = norm_signac_0hr,
  region = "chr5-138456812-138468286" # EGR1
) + ggtitle("0hr")

CoveragePlot(
  object = norm_signac_1hr,
  region = "chr5-138456812-138468286" # EGR1
) + ggtitle("1hr")

CoveragePlot(
  object = norm_signac_6hr,
  region = "chr5-138456812-138468286" # EGR1
) + ggtitle("6hr")


CoveragePlot(
  object = norm_signac_0hr,
  region = "chr11-27719113-27734038" # BDNF
) + ggtitle("0hr")

CoveragePlot(
  object = norm_signac_1hr,
  region = "chr11-27719113-27734038" # BDNF
) + ggtitle("1hr")

CoveragePlot(
  object = norm_signac_6hr,
  region = "chr11-27719113-27734038" # BDNF
) + ggtitle("6hr")


CoveragePlot(
  object = norm_signac_0hr,
  region = "chr6-6005375-6014978" # NRN1
) + ggtitle("0hr")

CoveragePlot(
  object = norm_signac_1hr,
  region = "chr6-6005375-6014978" # NRN1
) + ggtitle("1hr")

CoveragePlot(
  object = norm_signac_6hr,
  region = "chr6-6005375-6014978" # NRN1
) + ggtitle("6hr")


CoveragePlot(
  object = norm_signac_0hr,
  region = "chr7-101163942-101173124" # VGF
) + ggtitle("0hr")

CoveragePlot(
  object = norm_signac_1hr,
  region = "chr7-101163942-101173124" # VGF
) + ggtitle("1hr")

CoveragePlot(
  object = norm_signac_6hr,
  region = "chr7-101163942-101173124" # VGF
) + ggtitle("6hr")

