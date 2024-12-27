# Chuxuan Li 04/19/2023
# QC, normalize, and cluster aggregated 018-029 ATACseq data
# This code is adapted from 20230203_QC_normalize_cluster_029_ATAC_data.R

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
options(future.seed = TRUE)

load("./018-029_multiomic_obj_mapped_to_human_bc_only.RData")
load("./018-029_RNA_integrated_labeled_with_harmony.RData")

# QC ####
# compute nucleosome signal score per cell
DefaultAssay(multiomic_obj)
multiomic_obj <- NucleosomeSignal(object = multiomic_obj, assay = "ATAC")
# compute TSS enrichment score per cell
multiomic_obj <- TSSEnrichment(object = multiomic_obj)

VlnPlot(multiomic_obj, group.by = "group.ident",
        features = c("TSS.enrichment", "nucleosome_signal"))

hist(multiomic_obj$TSS.enrichment, breaks = 1000, xlim = c(0, 15), 
     main = "029 ATAC-seq data TSS enrichment distribution", xlab = "TSS enrichment")
hist(multiomic_obj$nucleosome_signal, breaks = 1000, xlim = c(0, 2), 
     main = "029 ATAC-seq data nucleosome signal distribution", 
     xlab = "Nucleosome signal")

sum(multiomic_obj$nucleosome_signal < 2 & multiomic_obj$TSS.enrichment > 2)
save(multiomic_obj, file = "018-029_multiomioc_obj_added_TSS_nucsig.RData")

# Idents(multiomic_obj) <- "nucleosome_signal"
# multiomic_obj <- subset(multiomic_obj, nucleosome_signal < 2)
# Idents(multiomic_obj) <- "TSS.enrichment"
# multiomic_obj <- subset(multiomic_obj, TSS.enrichment > 2)

multiomic_obj <- subset(multiomic_obj, nucleosome_signal < 2 & TSS.enrichment > 2)

# dimensional reduction ####
DefaultAssay(multiomic_obj) <- "ATAC"
multiomic_obj <- RunTFIDF(multiomic_obj)
multiomic_obj <- FindTopFeatures(multiomic_obj, min.cutoff = 'q0')
multiomic_obj <- RunSVD(multiomic_obj)
multiomic_obj <- RunUMAP(object = multiomic_obj, reduction = 'lsi', dims = 2:30, seed.use = 42)
multiomic_obj <- FindNeighbors(object = multiomic_obj, reduction = 'lsi', dims = 2:30)
multiomic_obj <- FindClusters(object = multiomic_obj, verbose = FALSE, 
                              algorithm = 3, random.seed = 42, resolution = 0.5)
DimPlot(object = multiomic_obj, cols = DiscretePalette(24, "alphabet2"),
        label = TRUE) + 
  NoLegend() +
  ggtitle("clustering of 029 ATACseq data") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

DimPlot(object = multiomic_obj, label = F, group.by = "group.ident") + 
  ggtitle("by library") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

multiomic_obj$RNA.fine.cell.type <- "unknown"
multiomic_obj$RNA.cell.type <- "unknown"
fine.types <- as.vector(unique(integrated_labeled$fine.cell.type))
types <- unique(integrated_labeled$cell.type)
for (i in 1:length(types)) {
  t <- types[i]
  barcodes <- str_extract(multiomic_obj@assays$ATAC@counts@Dimnames[[2]], 
                          "^[A-Z]+-") %in%
    str_extract(integrated_labeled@assays$RNA@counts@Dimnames[[2]][integrated_labeled$cell.type == t],
                "^[A-Z]+-")
  multiomic_obj$RNA.cell.type[barcodes] <- t
}
for (i in 1:length(fine.types)) {
  t <- fine.types[i]
  print(t)
  barcodes <- str_extract(multiomic_obj@assays$ATAC@counts@Dimnames[[2]], 
                          "^[A-Z]+-") %in%
    str_extract(integrated_labeled@assays$RNA@counts@Dimnames[[2]][integrated_labeled$fine.cell.type == t],
                "^[A-Z]+-")
  multiomic_obj$RNA.fine.cell.type[barcodes] <- t
}
unique(multiomic_obj$RNA.cell.type)
unique(multiomic_obj$RNA.fine.cell.type)

DimPlot(object = multiomic_obj, label = TRUE, group.by = "RNA.cell.type", 
        cols = c("#B33E52", "#E6D2B8", "#CCAA7A", "#54990F")) + 
  NoLegend() +
  ggtitle("029 ATAC-seq data\nlabeled by RNA cell types") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
DimPlot(object = multiomic_obj, label = TRUE, group.by = "RNA.fine.cell.type", 
        cols = c("#f29116", #?glut
                 "#B33E52", #GABA
                 "#E6D2B8", #nmglut
                 "#CCAA7A", #npglut
                 "#347545", #NPC
                 "#4293db"#unknown
        )) + 
  ggtitle("029 ATAC-seq data\nlabeled by RNA sub-cell types") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

times <- unique(integrated_labeled$time.ident)
multiomic_obj$time.ident <- ""
for (t in times) {
  barcodes <- str_extract(multiomic_obj@assays$ATAC@counts@Dimnames[[2]], 
                          "^[A-Z]+-") %in%
    str_extract(integrated_labeled@assays$RNA@counts@Dimnames[[2]][integrated_labeled$time.ident == t],
                "^[A-Z]+-")
  multiomic_obj$time.ident[barcodes] <- t
}
unique(multiomic_obj$time.ident)

multiomic_obj$timextype.ident <- "NA"
for (i in 1:length(types)){
  print(types[i])
  for (j in 1:length(times)){
    print(times[j])
    multiomic_obj$timextype.ident[multiomic_obj$time.ident == times[j] &
                                    multiomic_obj$RNA.cell.type == types[i]] <- 
      paste(types[i], times[j], sep = "_")
  }
}
sort(unique(multiomic_obj$timextype.ident))
save(multiomic_obj, file = "018-029_multiomic_obj_clustered_added_metadata.RData")
