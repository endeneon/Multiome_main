# Chuxuan Li 02/03/2023
# read and match 10x-aggr aggregated 029 ATACseq data with RNAseq data

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
counts <- Read10X_h5("/data/FASTQ/Duan_Project_029/GRCh38_only/Duan_Project_029_aggr/outs/filtered_feature_bc_matrix.h5")
frag <- "/data/FASTQ/Duan_Project_029/GRCh38_only/Duan_Project_029_aggr/outs/atac_fragments.tsv.gz"
# load RNA object
load("../Duan_project_029_RNA/029_RNA_integrated_labeled.RData")
# load annotation
load("../R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")

ATAC <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = frag,
  annotation = ens_use
)
ATAC # 309554 features for 234657 cells
sum(rowSums(ATAC@data))
# 1966629453
sum(rowSums(ATAC@counts))
# 1966629453
save(ATAC, file = "029_ATAC_chromatin_assay.RData")

# clean object ####
# remove rat genes
ATAC_obj <- CreateSeuratObject(counts = ATAC, assay = "ATAC")
use.intv <- ATAC_obj@assays$ATAC@counts@Dimnames[[1]][str_detect(ATAC_obj@assays$ATAC@counts@Dimnames[[1]], "^chr[0-9]+")]
ATAC_counts <- ATAC_obj@assays$ATAC@counts[(which(rownames(ATAC_obj@assays$ATAC@counts) %in% use.intv)), ]
ATAC_obj <- subset(ATAC_obj, features = rownames(ATAC_counts))

# combine with RNA seq object
multiomic_obj <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA")
ATAC_assay <- GetAssayData(object = ATAC_obj[["ATAC"]])
multiomic_obj[["ATAC"]] <- CreateChromatinAssay(counts = ATAC_assay, 
                                                fragments = Fragments(ATAC_obj),
                                                annotation = ens_use)
DefaultAssay(multiomic_obj) <- "ATAC"
rm(ATAC_assay, ATAC_counts)

# assign group and cell line idents to ATACseq obj
multiomic_obj$group.ident <- "unknown"
orig.idents <- as.character(sort(as.numeric(unique(str_extract(multiomic_obj@assays$ATAC@counts@Dimnames[[2]], "[0-9]+$")))))
orig.idents
multiomic_obj$orig.ident <- str_extract(multiomic_obj@assays$ATAC@counts@Dimnames[[2]], "[0-9]+$")
group.idents <- sort(unique(integrated_labeled$orig.ident))

for (i in 1:length(orig.idents)){
  o <- orig.idents[i]
  print(o)
  print(group.idents[i])
  multiomic_obj$group.ident[multiomic_obj$orig.ident == o] <- group.idents[i]
}
unique(multiomic_obj$group.ident)

multiomic_obj$cell.line.ident <- "unmatched"
for (i in 1:length(group.idents)){
  g <- group.idents[i]
  print(g)
  subRNA <- subset(integrated_labeled, orig.ident == g)
  lines <- unique(subRNA$cell.line.ident)
  for (l in lines) {
    linesubRNA <- subset(subRNA, cell.line.ident == l)
    multiomic_obj$cell.line.ident[str_extract(multiomic_obj@assays$ATAC@counts@Dimnames[[2]], 
                                              "^[A-Z]+-") %in%
                                    str_extract(linesubRNA@assays$RNA@counts@Dimnames[[2]],
                                                "^[A-Z]+-") &
                                    multiomic_obj$group.ident == g] <- l
  }
}
unique(multiomic_obj$group.ident)
unique(multiomic_obj$cell.line.ident)
multiomic_obj <- subset(multiomic_obj, cell.line.ident != "unmatched")
save(multiomic_obj, file = "029_clean_multiomic_obj_mapped_to_RNA.RData")
