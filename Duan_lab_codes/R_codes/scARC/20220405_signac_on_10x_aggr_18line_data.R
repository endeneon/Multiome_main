# Chuxuan Li 04/05/2022
# read and process 10x-aggr aggregated 18-line ATACseq data 

# init####
library(Seurat)
library(Signac)
library(sctransform)
library(glmGamPoi)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

library(patchwork)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)
library(stringr)
library(future)

plan("multisession", workers = 2)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# load data ####
# load ATAC filtered bc matrix
counts <- Read10X_h5("//data/FASTQ/Duan_Project_024/hybrid_output/hybrid_aggr_5groups_no2ndary/outs/filtered_feature_bc_matrix.h5")
frag <- "/data/FASTQ/Duan_Project_024/hybrid_output/hybrid_aggr_5groups_no2ndary/outs/atac_fragments.tsv.gz"
# load RNA object
load("../Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")
# load annotation
load("../R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")

ATAC <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = frag,
  annotation = ens_use
)
ATAC
# sum(rowSums(ATAC@data))
# 2206788642
# sum(rowSums(ATAC@counts))
# 2206788642

# remove rat genes ####
ATAC <- CreateSeuratObject(counts = ATAC, assay = "ATAC")
use.intv <- ATAC@assays$ATAC@counts@Dimnames[[1]][str_detect(ATAC@assays$ATAC@counts@Dimnames[[1]], "^chr")]
counts <- ATAC@assays$ATAC@counts[(which(rownames(ATAC@assays$ATAC@counts) %in% use.intv)), ]
ATAC <- subset(ATAC, features = rownames(counts))

# assign group idents to ATACseq obj first
orig.idents <- as.character(sort(as.numeric(unique(str_extract(ATAC@assays$ATAC@counts@Dimnames[[2]], "[0-9]+$")))))
orig.idents
group.idents <- list.files("/data/FASTQ/Duan_Project_024/hybrid_output/", pattern = "filtered_feature_bc_matrix.h5", recursive = T,
           full.names = T)
group.idents <- group.idents[1:15]
group.idents <- str_extract(group.idents, "[0-9]+-[0|1|6]")
group.idents
ATAC$group.ident <- "unknown"

for (i in 1:length(orig.idents)){
  o <- orig.idents[i]
  print(o)
  print(group.idents[i])
  ATAC$orig.ident[str_extract(ATAC@assays$ATAC@counts@Dimnames[[2]], "[0-9]+$") == o] <- o
  ATAC$group.ident[str_extract(ATAC@assays$ATAC@counts@Dimnames[[2]], "[0-9]+$") == o] <- group.idents[i]
}

unique(ATAC$group.ident)

# map to RNAseq cells ####
ATAC$cells.in.RNA <- "false"
unique(integrated_labeled$orig.ident)
for (i in 1:length(group.idents)){
  print(group.idents[i])
  subRNA <- subset(integrated_labeled, orig.ident == group.idents[i])

  ATAC$cells.in.RNA[str_extract(ATAC@assays$ATAC@counts@Dimnames[[2]], "^[A-Z]+-") %in%
                      str_extract(integrated_labeled@assays$RNA@counts@Dimnames[[2]], "^[A-Z]+-") &
                      ATAC$group.ident == group.idents[i]] <- "true"
}

unique(ATAC$cells.in.RNA)
sum(ATAC$cells.in.RNA == "true") #146412
sum(ATAC$cells.in.RNA == "false") #48519
ATAC <- subset(ATAC, cells.in.RNA == "true")

# # select only cells that are mapped to human barcodes ####
# # read barcode .best files
setwd("/data/FASTQ/Duan_Project_024/hybrid_output")
pathlist <- list.files(path = ".", full.names = T, pattern = ".best.tsv", recursive = T)
pathlist <- pathlist[str_detect(pathlist, "demux_0.01", T)]
barcode_list <- vector(mode = "list", length = length(pathlist))
for (i in 1:length(pathlist)){
  barcode_list[[i]] <- read.delim(pathlist[i], header = F, row.names = NULL)
  colnames(barcode_list[[i]]) <- c("barcode", "line")
}
ATAC$cell.line.ident <- "unmatched"
group.idents
pathlist
# sep lines in each group's barcode list, then assign line by barcode
for (i in 1:length(barcode_list)){
  # total 15 barcode lists for 15 libraries, each lib has 3-4 lines
  lines <- unique(barcode_list[[i]]$line)
  print(lines)
  # initialize
  for (j in 1:length(lines)){
    # each lib has a list of lines and corresponding barcodes
    line_spec_barcodes <- barcode_list[[i]]$barcode[barcode_list[[i]]$line == lines[j]]
    ATAC$cell.line.ident[str_extract(ATAC@assays$ATAC@counts@Dimnames[[2]], "^[A-Z]+") %in%
                         str_extract(line_spec_barcodes, "^[A-Z]+") &
                         ATAC$group.ident == group.idents[i]] <- lines[j]
  }
  print(unique(ATAC$cell.line.ident))
}
print(unique(ATAC$cell.line.ident))
ATAC <- subset(ATAC, cell.line.ident != "unmatched")

setwd("~/NVME/scARC_Duan_018/Duan_project_024_ATAC")
save(ATAC, file = "18line_10xaggred_removed_rat_genes_mapped_to_demuxed_barcodes.RData")

# QC ####
# compute nucleosome signal score per cell
ATAC <- NucleosomeSignal(object = ATAC)

# compute TSS enrichment score per cell
ATAC <- TSSEnrichment(object = ATAC#, fast = T)
)
ATAC <- subset(
  x = ATAC,
  subset = nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# dimensional reduction ####
ATAC <- RunTFIDF(ATAC)
ATAC <- FindTopFeatures(ATAC, min.cutoff = 'q0')
ATAC <- RunSVD(ATAC)
ATAC <- RunUMAP(object = ATAC, reduction = 'lsi', dims = 2:30, seed.use = 42)
ATAC <- FindNeighbors(object = ATAC, reduction = 'lsi', dims = 2:30)
ATAC <- FindClusters(object = ATAC, verbose = FALSE, 
                     algorithm = 3, random.seed = 42)
# assign time point ident
ATAC$time.ident <- paste0(str_extract(ATAC$group.ident, "[0|1|6]$"), "hr")
unique(ATAC$time.ident)

DimPlot(object = ATAC, label = TRUE) + 
  NoLegend() +
  ggtitle("clustering of 18-line ATACseq data")

DimPlot(object = ATAC, label = F, group.by = "group.ident") + 
  ggtitle("by library")

DimPlot(object = ATAC, label = TRUE, group.by = "time.ident") + 
  NoLegend() +
  ggtitle("by time")

save(ATAC, file = "18line_10xaggred_normalized_clustered_obj.RData")

# gene activity ####
gene.activities <- GeneActivity(ATAC)
ATAC[["gact"]] <- CreateAssayObject(counts = gene.activities)
ATAC <- NormalizeData(
  object = ATAC,
  assay = 'gact',
  normalization.method = 'LogNormalize',
  scale.factor = median(ATAC$nCount_ATAC)
)

FeaturePlot(
  object = ATAC,
  features = c("GAD1", "GAD2", "SLC17A6", "SLC17A7"),
  # pt.size = 0.1,
  max.cutoff = "q95",
  ncol = 2)

FeaturePlot(
  object = ATAC,
  features = c("NEFM", "CUX2"),
  pt.size = 0.1,
  max.cutoff = "q95",
  ncol = 2)

subtype_markers <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
                     "EBF1", "SERTAD4", # striatal
                     #"FOXG1", # forebrain
                     "TBR1", 
                     #"FOXP2", 
                     "TLE4", # pallial glutamatergic
                     "FEZF2", "ADCYAP1", "TCERG1L", 
                     #"SEMA3E", # subcerebral
                     "POU3F2", "CUX1", "BCL11B",  # cortical
                     "LHX2", # general cortex
                     #"EOMES", # hindbrain
                     "NPY", "SST", "DLX2", 
                     #"DLX5", # inhibitory
                     "SATB2", "CBLN2", "CUX2", "NEFM", # excitatory
                     #"VIM", 
                     "SOX2", "NES" #NPC
                     )

StackedVlnPlot(ATAC, subtype_markers)

# apply RNAseq cell types ####
types <- as.vector(unique(integrated_labeled$cell.type))
ATAC$RNA.cell.type <- "unknown"
for (i in 1:length(types)){
  print(types[i])
  subRNA <- subset(integrated_labeled, cell.type == types[i])
  ATAC$RNA.cell.type[str_extract(ATAC@assays$ATAC@counts@Dimnames[[2]], "^[A-Z]+-") %in% 
                       str_extract(subRNA@assays$RNA@counts@Dimnames[[2]], "^[A-Z]+-")] <- types[i]
}
unique(ATAC$RNA.cell.type)
sum(ATAC$RNA.cell.type == "unknown") #8899
DimPlot(object = ATAC, 
        label = TRUE, 
        group.by = "RNA.cell.type") + 
  # NoLegend() +
  ggtitle("ATACseq data projected by RNAseq cell types")


ATAC$cell.type <- ATAC$RNA.cell.type
ATAC$cell.type[ATAC$RNA.cell.type %in% c("SST_pos_GABA", "SEMA3E_pos_GABA", "GABA")] <- "GABA"
ATAC$cell.type[ATAC$RNA.cell.type == "unknown glut"] <- "unknown"

ATAC$RNA.cell.type[ATAC$RNA.cell.type == "unknown glut"] <- "unknown"
ATAC$RNA.cell.type <-
  factor(ATAC$RNA.cell.type,
         levels = c("NEFM_pos_glut",
                    "GABA",
                    "NEFM_neg_glut", 
                    "SEMA3E_pos_GABA",
                    "unknown",
                    "SST_pos_GABA"))
DimPlot(object = ATAC, 
        label = F, 
        cols = brewer.pal(n = 6,
                          name = "Dark2"),
        group.by = "RNA.cell.type") + 
  # NoLegend() +
  ggtitle("ATACseq data projected by RNAseq cell types")

sum(ATAC$cell.type == "unknown") #17867
sum(ATAC$cell.type == "GABA") #52672
sum(ATAC$cell.type == "NEFM_neg_glut") #40477
sum(ATAC$cell.type == "NEFM_pos_glut") #33554


# assign time x cell type ident
types <- unique(ATAC$cell.type)
times <- unique(ATAC$time.ident)
ATAC$timextype.ident <- "NA"

for (i in 1:length(types)){
  print(types[i])
  for (j in 1:length(times)){
    print(times[j])
    ATAC$timextype.ident[ATAC$time.ident == times[j] &
                           ATAC$cell.type == types[i]] <- 
      paste(types[i], times[j], sep = "_")
  }
}
unique(ATAC$timextype.ident)

save(ATAC, file = "18line_10xaggred_labeled_by_RNA_obj.RData")

# plot other dimplots
DimPlot(object = ATAC, label = F, 
        group.by = "time.ident") +
  ggtitle("By time")
DimPlot(object = ATAC, label = F, 
        group.by = "group.ident") +
  ggtitle("By library")

# vlnplot for percent mt
sum(str_detect(rownames(multiomic_obj), "^chrM")) #0
sum(str_detect(rownames(multiomic_obj), "^MT-")) #0

multiomic_obj[["percent.mt"]] <- PercentageFeatureSet(multiomic_obj, assay = "ATAC",
                                                   pattern = "^MT-")
DefaultAssay(multiomic_obj)
Idents(multiomic_obj) <- "rat.ident"
VlnPlot(multiomic_obj, features = c("nFeature_ATAC", "nCount_ATAC", "percent.mt"), 
        ncol = 3, fill.by = "feature", pt.size = 0.1)

sum(integrated$percent.mt > 10) #20731
sum(integrated$percent.mt > 20) #750