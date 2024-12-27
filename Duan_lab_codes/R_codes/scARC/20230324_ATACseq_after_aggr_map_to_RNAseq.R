# Chuxuan Li 03/24/2023
# read and match 10x-aggr aggregated 018-029 ATACseq data with RNAseq data

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
options(future.globals.maxSize = 5368709120)

# load data ####
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/018-029_RNA_integrated_labeled_with_harmony.RData")

library(Matrix)
counts1 <- Read10X("./018_to_029_combined/outs/filtered_feature_bc_matrix/part1/gzfiles/")
counts2 <- Read10X("./018_to_029_combined/outs/filtered_feature_bc_matrix/part2/gzfiles/")
counts3 <- Read10X("./018_to_029_combined/outs/filtered_feature_bc_matrix/part3/gzfiles/")
frag <- "./018_to_029_combined/outs/atac_fragments.tsv.gz"
counts <- counts1$Peaks + counts2$Peaks + counts3$Peaks

# load annotation
load("../R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")

ATAC <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = frag,
  annotation = ens_use
)
ATAC # 225934 features for 795150 cells
sum(rowSums(ATAC@data))
# 4185464619
save(ATAC, file = "018-029_ATAC_chromatin_assay.RData")


# clean object ####
# remove rat genes
ATAC_obj <- CreateSeuratObject(counts = ATAC, assay = "ATAC")
use.intv <- ATAC_obj@assays$ATAC@counts@Dimnames[[1]][str_detect(ATAC_obj@assays$ATAC@counts@Dimnames[[1]], "^chr[0-9]+")]
ATAC_counts <- ATAC_obj@assays$ATAC@counts[(which(rownames(ATAC_obj@assays$ATAC@counts) %in% use.intv)), ]
ATAC_obj <- subset(ATAC_obj, features = rownames(ATAC_counts))

# combine with RNA seq object
counts <- counts1$`Gene Expression` + counts2$`Gene Expression` + counts3$`Gene Expression`
multiomic_obj <- CreateSeuratObject(counts = counts, assay = "RNA")
ATAC_assay <- GetAssayData(object = ATAC_obj[["ATAC"]])
multiomic_obj[["ATAC"]] <- CreateChromatinAssay(counts = ATAC_assay, 
                                                fragments = Fragments(ATAC_obj),
                                                annotation = ens_use)
DefaultAssay(multiomic_obj) <- "ATAC"

rm(counts1, counts2, counts3, counts, ATAC_assay, ATAC_counts)
gc()

# assign group and cell line idents to ATACseq obj
multiomic_obj$group.ident <- "unknown"
group_order_csv <- read_csv("Duan_Project_018_to_029_samples.csv")
group_order <- group_order_csv$library_id
orig.idents <- as.character(sort(as.numeric(unique(str_extract(multiomic_obj@assays$ATAC@counts@Dimnames[[2]], "[0-9]+$")))))
orig.idents
multiomic_obj$orig.ident <- str_extract(multiomic_obj@assays$ATAC@counts@Dimnames[[2]], "[0-9]+$")
for (i in 1:length(orig.idents)){
  o <- orig.idents[i]
  print(o)
  print(group_order[i])
  multiomic_obj$group.ident[multiomic_obj$orig.ident == o] <- group_order[i]
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

unique(multiomic_obj$cell.line.ident)
sum(multiomic_obj$cell.line.ident == "unmatched")
multiomic_obj <- subset(multiomic_obj, cell.line.ident != "unmatched")
save(multiomic_obj, file = "018-029_multiomic_obj_mapped_to_human_bc_only.RData")
