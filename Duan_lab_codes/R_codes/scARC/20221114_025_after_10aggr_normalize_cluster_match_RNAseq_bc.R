# Chuxuan Li 11/14/2022
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
set.seed(2001)
options(future.globals.maxSize = 229496729600)

# load data ####
# load ATAC filtered bc matrix
counts <- Read10X_h5("~/NVME/scARC_Duan_025_GRCh38/Duan_Project_aggregated/outs/filtered_feature_bc_matrix.h5")
frag <- "~/NVME/scARC_Duan_025_GRCh38/Duan_Project_aggregated/outs/atac_fragments.tsv.gz"
# load RNA object
load("~/NVME/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/integrated_labeled.RData")
# load annotation
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")

ATAC <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = frag,
  annotation = ens_use
)
ATAC # 289272 features for 189397 cells
# sum(rowSums(ATAC@data))
# 2066737480
# sum(rowSums(ATAC@counts))
# 2066737480

# remove rat genes ####
ATAC <- CreateSeuratObject(counts = ATAC, assay = "ATAC") #23.7 GB
use.intv <- ATAC@assays$ATAC@counts@Dimnames[[1]][str_detect(ATAC@assays$ATAC@counts@Dimnames[[1]], "^chr")]
counts <- ATAC@assays$ATAC@counts[(which(rownames(ATAC@assays$ATAC@counts) %in% use.intv)), ]
ATAC <- subset(ATAC, features = rownames(counts))

# assign group idents to ATACseq obj
orig.idents <- as.character(sort(as.numeric(unique(str_extract(ATAC@assays$ATAC@counts@Dimnames[[2]], "[0-9]+$")))))
orig.idents
ATAC$orig.ident <- str_extract(ATAC@assays$ATAC@counts@Dimnames[[2]], "[0-9]+$")
group.idents <- list.files("/nvmefs/scARC_Duan_025_GRCh38/", pattern = "filtered_feature_bc_matrix.h5", recursive = T,
                           full.names = T)
group.idents <- group.idents[2:16]
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
                      str_extract(subRNA@assays$RNA@counts@Dimnames[[2]], "^[A-Z]+-") &
                      ATAC$group.ident == group.idents[i]] <- "true"
}
unique(ATAC$cells.in.RNA)
sum(ATAC$cells.in.RNA == "true") #147551
sum(ATAC$cells.in.RNA == "false") #50186
ATAC_clean <- subset(ATAC, cells.in.RNA == "true")

# # select only cells that are mapped to human barcodes ####
ATAC_clean$cell.line.ident <- "unmatched"
groups <- unique(integrated_labeled$orig.ident)

# sep lines in each group's barcode list, then assign line by barcode
for (g in groups) {
  cat("group: ", g, "\n")
  subRNA <- subset(integrated_labeled, orig.ident == g)
  lines <- unique(subRNA$cell.line.ident)
  for (l in lines) {
    subsubRNA <- subset(subRNA, cell.line.ident == l)
    bc <- subsubRNA@assays$RNA@counts@Dimnames[[2]]
    bc <- str_remove(bc, "-[0-9]_[0-9]+")
    print(sum(ATAC_clean$group.ident == g &
          str_remove(ATAC_clean@assays$ATAC@counts@Dimnames[[2]], "-[0-9]+") %in% bc))
    ATAC_clean$cell.line.ident[ATAC_clean$group.ident == g &
                                 str_remove(ATAC_clean@assays$ATAC@counts@Dimnames[[2]], "-[0-9]+") %in% bc] <- l
  }
}
print(unique(ATAC_clean$cell.line.ident))
ATAC_clean <- subset(ATAC_clean, cell.line.ident != "unmatched")
save(ATAC_clean, file = "025_ATAC_clean_obj_mapped_human_only.RData")

# QC ####
# compute nucleosome signal score per cell
ATAC_clean <- NucleosomeSignal(object = ATAC_clean)
# compute TSS enrichment score per cell
ATAC_clean <- TSSEnrichment(object = ATAC_clean)

VlnPlot(ATAC_clean, group.by = "group.ident",
        features = c("TSS.enrichment", "nucleosome_signal"))

# combine RNA and ATAC ####
multiomic_obj <- integrated_labeled
for (i in 1:length(groups)){
  g <- groups[i]
  chari <- as.character(i)
  print(g)
  multiomic_obj@assays$RNA@counts@Dimnames[[2]][multiomic_obj$orig.ident == g] <-
    str_replace(multiomic_obj@assays$RNA@counts@Dimnames[[2]][multiomic_obj$orig.ident == g],
                "[0-9]+_[0-9]+", chari)
  multiomic_obj@assays$RNA@data@Dimnames[[2]][multiomic_obj$orig.ident == g] <-
    str_replace(multiomic_obj@assays$RNA@data@Dimnames[[2]][multiomic_obj$orig.ident == g],
                "[0-9]+_[0-9]+", chari)
  multiomic_obj@assays$SCT@counts@Dimnames[[2]][multiomic_obj$orig.ident == g] <-
    str_replace(multiomic_obj@assays$SCT@counts@Dimnames[[2]][multiomic_obj$orig.ident == g],
                "[0-9]+_[0-9]+", chari)
  multiomic_obj@assays$SCT@data@Dimnames[[2]][multiomic_obj$orig.ident == g] <-
    str_replace(multiomic_obj@assays$SCT@data@Dimnames[[2]][multiomic_obj$orig.ident == g],
                "[0-9]+_[0-9]+", chari)
  multiomic_obj@assays$integrated@data@Dimnames[[2]][multiomic_obj$orig.ident == g] <-
    str_replace(multiomic_obj@assays$integrated@data@Dimnames[[2]][multiomic_obj$orig.ident == g],
                "[0-9]+_[0-9]+", chari)
}
multiomic_obj[["ATAC"]] <- CreateAssayObject(GetAssayData(ATAC_clean, slot = "counts"))
multiomic_obj$RNA.clusters <- multiomic_obj$seurat_clusters
multiomic_obj$nucleosome_signal <- ATAC_clean$nucleosome_signal
multiomic_obj$nucleosome_percentile <- ATAC_clean$nucleosome_percentile
multiomic_obj$TSS.enrichment <- ATAC_clean$TSS.enrichment
multiomic_obj$TSS.percentile <- ATAC_clean$TSS.percentile
hist(multiomic_obj$TSS.enrichment)
sum(multiomic_obj$nucleosome_signal < 2 & multiomic_obj$TSS.enrichment > 2)

Idents(multiomic_obj) <- "nucleosome_signal"
multiomic_obj <- subset(multiomic_obj, nucleosome_signal < 2)
Idents(multiomic_obj) <- "TSS.enrichment"
multiomic_obj <- subset(multiomic_obj, TSS.enrichment > 2)

# dimensional reduction ####
ATAC_clean <- subset(ATAC_clean, TSS.enrichment > 2 & nucleosome_signal < 2)
DefaultAssay(ATAC_clean) <- "ATAC"
ATAC_clean <- RunTFIDF(ATAC_clean)
ATAC_clean <- FindTopFeatures(ATAC_clean, min.cutoff = 'q0')
ATAC_clean <- RunSVD(ATAC_clean)
ATAC_clean <- RunUMAP(object = ATAC_clean, reduction = 'lsi', dims = 2:30, seed.use = 42)
ATAC_clean <- FindNeighbors(object = ATAC_clean, reduction = 'lsi', dims = 2:30)
ATAC_clean <- FindClusters(object = ATAC_clean, verbose = FALSE, 
                     algorithm = 3, random.seed = 42)
DimPlot(object = ATAC_clean, cols = DiscretePalette(24, "alphabet2"),
        label = TRUE) + 
  NoLegend() +
  ggtitle("clustering of 18-line ATACseq data")
ATAC_clean$RNA.cell.type <- ""
types <- unique(integrated_labeled$cell.type.forplot)
for (i in 1:length(types)){
  t <- types[i]
  bc <- integrated_labeled@assays$RNA@counts@Dimnames[[2]][integrated_labeled$cell.type.forplot == t]
  bc <- str_remove(bc, "-[0-9]_[0-9]+")
  ATAC_clean$RNA.cell.type[str_remove(ATAC_clean@assays$ATAC@counts@Dimnames[[2]], "-[0-9]+") %in% bc] <- t
}
ATAC_clean$time.ident <- "0hr"
ATAC_clean$time.ident[str_detect(ATAC_clean$group.ident, "-1")] <- "1hr"
ATAC_clean$time.ident[str_detect(ATAC_clean$group.ident, "-6")] <- "6hr"

DimPlot(object = ATAC_clean, label = F, group.by = "group.ident") + 
  ggtitle("by library")
DimPlot(object = ATAC_clean, label = TRUE, group.by = "RNA.cell.type", 
        cols = c("#B33E52", "#E6D2B8", "#CCAA7A", "#54990F")) + 
  NoLegend() +
  ggtitle("labeled by RNA cell types")

types <- unique(ATAC_clean$RNA.cell.type)
times <- unique(ATAC_clean$time.ident)
ATAC_clean$timextype.ident <- "NA"

for (i in 1:length(types)){
  print(types[i])
  for (j in 1:length(times)){
    print(times[j])
    ATAC_clean$timextype.ident[ATAC_clean$time.ident == times[j] &
                                 ATAC_clean$RNA.cell.type == types[i]] <- 
      paste(types[i], times[j], sep = "_")
  }
}
unique(ATAC_clean$timextype.ident)
save(ATAC_clean, file = "025_ATAC_clean_clustered_obj.RData")
