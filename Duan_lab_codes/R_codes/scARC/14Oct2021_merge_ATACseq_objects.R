## Chuxuan Li 10/14/2021
# ATACseq merging dataset method 2: use individual ones and use Signac 
# clustering of ATACseq data


# init
library(Seurat)
#library(sctransform)
#library(glmGamPoi)
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
library(uwot)

# library(reticulate)

# conda_list()
# use_condaenv("rstudio")


set.seed(1886)
# set threads and parallelization
plan("multisession", workers = 6)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

# read peak file from the aggregated dataset
aggregated_bed <- read.table(
  file = "/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/atac_peaks_cleaned.bed",
  col.names = c("chr", "start", "end")
)
# convert into granges object
gr <- makeGRangesFromDataFrame(aggregated_bed)
# filter bad peaks
peakwidths <- width(gr)
hist(peakwidths,
     breaks = 1000)
gr <- gr[peakwidths  < 1100 & peakwidths > 20]
gr

# read in csv files containing metadata
meta_lst <- vector(mode = "list", length = 6L)
path_list <- c("/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_0/outs/per_barcode_metrics.csv",
               "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_1/outs/per_barcode_metrics.csv",
               "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_6/outs/per_barcode_metrics.csv",
               "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_0/outs/per_barcode_metrics.csv",
               "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_1/outs/per_barcode_metrics.csv",
               "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_6/outs/per_barcode_metrics.csv")

for (i in 1:length(path_list)) {
  # read in csv file containing barcodes
  md <- read.table(
    file = path_list[i],
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )
  # QC
  #print(md$atac_fragments)
  md <- md[md$atac_fragments > 500, ]
  # put in list
  meta_lst[[i]] <- md
}

# read in the fragment files in a loop
frag_lst <- vector(mode = "list", length = 6L)
path_list <- c("/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_0/outs/atac_fragments.tsv.gz",
               "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_1/outs/atac_fragments.tsv.gz",
               "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_6/outs/atac_fragments.tsv.gz",
               "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_0/outs/atac_fragments.tsv.gz",
               "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_1/outs/atac_fragments.tsv.gz",
               "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_6/outs/atac_fragments.tsv.gz")

for (i in 1:length(path_list)) {
  # read in csv file containing barcodes
  fg <- CreateFragmentObject(
    path = path_list[i],
    cells = rownames(meta_lst[[i]]))
  frag_lst[[i]] <- fg
}

# quantify peaks by creating a matrix of peaks * cell for each sample, then create Seurat objects
obj_lst <- vector(mode = "list", length = 6L)
ori_labels <- vector(mode = "list", length = 6L)
for (i in 1:length(frag_lst)){
  f <- frag_lst[[i]]
  m <- meta_lst[[i]]
  p <- path_list[i]
  ori <- gsub("^.*([0-9]_[0-9])/.*", "\\1", p)
  #print(ori)
  ori_labels[i] <- ori
  counts <- FeatureMatrix(
    fragments = f,
    features = gr,
    cells = rownames(f)
  )
  assay <- CreateChromatinAssay(counts, fragments = f)
  obj <- CreateSeuratObject(assay, assay = "ATAC", meta.data = m)
  obj$origin <- ori
  obj_lst[[i]] <- obj
}

l <- as.list(obj_lst[2:6])

# merge the datasets
merged_signac <- merge(
  x = obj_lst[[1]],
  y = l,
  add.cell.ids = c("2_0", "2_1", "2_6", "8_0", "8_1", "8_6")
)

# check the chromatin assay 
merged_signac[["ATAC"]] 
# ChromatinAssay data with 111255 features for 2485126 cells
# Variable features: 0 
# Genome: 
#   Annotation present: FALSE 
# Motifs present: FALSE 
# Fragment files: 6 

# do normalization
merged_signac_processed <- RunTFIDF(merged_signac)
merged_signac_processed <- FindTopFeatures(merged_signac_processed, min.cutoff = 20)
merged_signac_processed <- RunSVD(merged_signac_processed)

# check which components are associated with technical (nonbiological) variance (sequencing depth)
DepthCor(merged_signac)

#save("merged_signac", file = "merged_signac.RData")


# clustering and downstream analysis
merged_signac_processed <- RunUMAP(object = merged_signac_processed, 
                                   reduction = 'lsi',
                                   dims = 2:30, 
                                   # umap.method = "umap-learn",
                                   #umap.method = "uwot-learn",
                                   #metric = "correlation",
                                   seed.use = 1234) # remove the first component
merged_signac_processed <- FindNeighbors(object = merged_signac_processed, 
                                         reduction = 'lsi', 
                                         dims = 2:30) # remove the first component
merged_signac_processed <- FindClusters(object = merged_signac_processed, 
                                        verbose = FALSE, 
                                        algorithm = 3, 
                                        resolution = 0.5,
                                        random.seed = 1889)

DimPlot(object = merged_signac_filtered, label = TRUE) + 
  NoLegend() +
  ggtitle("clustering of ATACseq data using Signac merge")

# check how well the objects are blended
DimPlot(merged_signac, group.by = 'origin', pt.size = 0.1)

save("merged_signac_processed", file = "merged_signac_processed.RData")


# load labeled RNAseq objects
load("/home/cli/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/integrated_labeled_updated.RData")

# create another attribute to project RNAseq cell type ident onto the ATACseq data
merged_signac$RNA.cell.type.ident <- NA
# find all barcodes corresponding to each cell type

for (i in levels(integrated_renamed@active.ident)){
  barcodes <- str_sub(merged_signac@assays$RNA@counts@Dimnames[[2]][integrated_renamed@active.ident %in% i], 
                      end = -7L)
  
  merged_signac$RNA.cell.type.ident[merged_signac$trimmed.barcodes %in% barcodes] <- i
}


# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = merged_signac, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.cell.type.ident") + 
  ggtitle("clustering of ATACseq data with cell type projection")




# plot peak tracks as an example
CoveragePlot(
  object = merged_signac,
  region = "chr11-22334000-22340000", # SLC17A6
  extend.upstream = 1,
  extend.downstream = 1
)

CoveragePlot(
  object = merged_signac,
  region = "chr19-49427411-49431451", # SLC17A7
  extend.upstream = 10000,
  extend.downstream = 1
)

CoveragePlot(
  object = merged_signac,
  region = "chr10-26204609-26224140" # GAD2
)

CoveragePlot(
  object = merged_signac,
  region = "chr2-190871076-190883628" # GLS
)

CoveragePlot(
  object = merged_signac,
  region = "chr7-101805750-101822128" # CUX1
)

CoveragePlot(
  object = merged_signac,
  region = "chr12-111026276-111037993" # CUX2
)

CoveragePlot(
  object = merged_signac,
  region = "chr14-28759074-28770857", # FOXG1
  extend.upstream = 10000
)

CoveragePlot(
  object = merged_signac,
  region = "chr3-181703353-181714656" # SOX2
)

CoveragePlot(
  object = merged_signac,
  region = "chr1-156675466-156688433" # NES
)

CoveragePlot(
  object = merged_signac,
  region = "chr14-75269196-75280682" # FOS
)

CoveragePlot(
  object = merged_signac,
  region = "chr5-138456812-138468286" # EGR1
)

CoveragePlot(
  object = merged_signac,
  region = "chr11-27719113-27734038" # BDNF
)

CoveragePlot(
  object = merged_signac,
  region = "chr6-6005375-6014978" # NRN1
)

CoveragePlot(
  object = merged_signac,
  region = "chr7-101163942-101173124" # VGF
)

# separate time points, look at response markers for each time point
merged_signac_0hr <- subset(merged_signac, subset = time.ident == "0hr")
merged_signac_1hr <- subset(merged_signac, subset = time.ident == "1hr")
merged_signac_6hr <- subset(merged_signac, subset = time.ident == "6hr")

CoveragePlot(
  object = merged_signac_0hr,
  region = "chr14-75269196-75280682" # FOS
) + ggtitle("0hr")

CoveragePlot(
  object = merged_signac_1hr,
  region = "chr14-75269196-75280682" # FOS
) + ggtitle("1hr")

CoveragePlot(
  object = merged_signac_6hr,
  region = "chr14-75269196-75280682" # FOS
) + ggtitle("6hr")



CoveragePlot(
  object = merged_signac_0hr,
  region = "chr5-138456812-138468286" # EGR1
) + ggtitle("0hr")

CoveragePlot(
  object = merged_signac_1hr,
  region = "chr5-138456812-138468286" # EGR1
) + ggtitle("1hr")

CoveragePlot(
  object = merged_signac_6hr,
  region = "chr5-138456812-138468286" # EGR1
) + ggtitle("6hr")


CoveragePlot(
  object = merged_signac_0hr,
  region = "chr11-27719113-27734038" # BDNF
) + ggtitle("0hr")

CoveragePlot(
  object = merged_signac_1hr,
  region = "chr11-27719113-27734038" # BDNF
) + ggtitle("1hr")

CoveragePlot(
  object = merged_signac_6hr,
  region = "chr11-27719113-27734038" # BDNF
) + ggtitle("6hr")


CoveragePlot(
  object = merged_signac_0hr,
  region = "chr6-6005375-6014978" # NRN1
) + ggtitle("0hr")

CoveragePlot(
  object = merged_signac_1hr,
  region = "chr6-6005375-6014978" # NRN1
) + ggtitle("1hr")

CoveragePlot(
  object = merged_signac_6hr,
  region = "chr6-6005375-6014978" # NRN1
) + ggtitle("6hr")


CoveragePlot(
  object = merged_signac_0hr,
  region = "chr7-101163942-101173124" # VGF
) + ggtitle("0hr")

CoveragePlot(
  object = merged_signac_1hr,
  region = "chr7-101163942-101173124" # VGF
) + ggtitle("1hr")

CoveragePlot(
  object = merged_signac_6hr,
  region = "chr7-101163942-101173124" # VGF
) + ggtitle("6hr")

