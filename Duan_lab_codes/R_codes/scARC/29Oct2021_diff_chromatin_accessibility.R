# Chuxuan Li 10/29/2021
# Do differential chromatin accessibility analysis on fine cell types by time

# init
library(Seurat)
library(readr)
library(ggplot2)
library(future)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)
library(cowplot)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)

library(stringr)

set.seed(1886)
# set threads and parallelization
plan("multisession", workers = 3)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

# load data (10x aggregated version)
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/aggr_ATACseq_with_all_labels.RData")


# check differential accessibility
# change back to working with peaks instead of gene activities
DefaultAssay(aggr_GABA) <- 'peaks'
DefaultAssay(aggr_NmCo_glut) <- 'peaks'
DefaultAssay(aggr_NpCm_glut) <- 'peaks'

Idents(aggr_GABA) <- "time.ident"

GABA_0v1_peaks <- FindMarkers(
  object = aggr_GABA,
  ident.1 = "1hr",
  ident.2 = "0hr",
  min.pct = 0, 
  logfc.threshold = 0,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)
head(GABA_0v1_peaks)

GABA_0v6_peaks <- FindMarkers(
  object = aggr_GABA,
  ident.1 = "6hr",
  ident.2 = "0hr",
  min.pct = 0, 
  logfc.threshold = 0,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

Idents(aggr_NmCo_glut) <- "time.ident"

NmCp_glut_0v1_peaks <- FindMarkers(
  object = aggr_NmCo_glut,
  ident.1 = "1hr",
  ident.2 = "0hr",
  min.pct = 0, 
  logfc.threshold = 0,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

NmCp_glut_0v6_peaks <- FindMarkers(
  object = aggr_NmCo_glut,
  ident.1 = "6hr",
  ident.2 = "0hr",
  min.pct = 0, 
  logfc.threshold = 0,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

Idents(aggr_NpCm_glut) <- "time.ident"

NpCm_glut_0v1_peaks <- FindMarkers(
  object = aggr_NpCm_glut,
  ident.1 = "1hr",
  ident.2 = "0hr",
  min.pct = 0, 
  logfc.threshold = 0,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

NpCm_glut_0v6_peaks <- FindMarkers(
  object = aggr_NpCm_glut,
  ident.1 = "6hr",
  ident.2 = "0hr",
  min.pct = 0, 
  logfc.threshold = 0,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

save.image(file = "29Oct2021.RData")

write.table(GABA_0v1_peaks, 
            file = "GABA_0v1hr_peaks.txt",
            quote = F,
            sep = "\t",
            row.names = T,
            col.names = T)
write.table(GABA_0v6_peaks, 
            file = "GABA_6v0hr_peaks.txt",
            quote = F,
            sep = "\t",
            row.names = T,
            col.names = T)
write.table(NmCp_glut_0v1_peaks, 
            file = "NmCp_glut_1v0hr_peaks.txt",
            quote = F,
            sep = "\t",
            row.names = T,
            col.names = T)
write.table(NmCp_glut_0v6_peaks, 
            file = "NmCp_glut_6v0hr_peaks.txt",
            quote = F,
            sep = "\t",
            row.names = T,
            col.names = T)
write.table(NpCm_glut_0v1_peaks, 
            file = "NpCm_glut_1v0hr_peaks.txt",
            quote = F,
            sep = "\t",
            row.names = T,
            col.names = T)
write.table(NpCm_glut_0v6_peaks, 
            file = "NpCm_glut_6v0hr_peaks.txt",
            quote = F,
            sep = "\t",
            row.names = T,
            col.names = T)

# get peaks from other cell types too
aggr_signac_unfiltered$fine.cell.type[aggr_signac_unfiltered$RNA.cell.type.ident %in% "SEMA3E+ glut"] <- "SEMA3E+ glut"
SEMA3Ep_glut <- subset(aggr_signac_unfiltered, subset = fine.cell.type %in% "SEMA3E+ glut")

aggr_signac_unfiltered$fine.cell.type[aggr_signac_unfiltered$RNA.cell.type.ident %in% c("MAP2+ NPC", "NPC")] <- "NPC"
NPC <- subset(aggr_signac_unfiltered, subset = fine.cell.type %in% "NPC")

aggr_signac_unfiltered$fine.cell.type[aggr_signac_unfiltered$RNA.cell.type.ident %in% c("immature forebrain glut", 
                                                                          "forebrain NPC")] <- "forebrain"
forebrain <- subset(aggr_signac_unfiltered, subset = fine.cell.type %in% "forebrain")

aggr_signac_unfiltered$fine.cell.type[aggr_signac_unfiltered$RNA.cell.type.ident == "subcerebral immature neuron"] <- "subcerebral immature neuron"
subcerebral <- subset(aggr_signac_unfiltered, subset = fine.cell.type %in% "subcerebral immature neuron")

aggr_signac_unfiltered$fine.cell.type[aggr_signac_unfiltered$RNA.cell.type.ident == "immature neuron"] <- "immature neuron"
immature_neuron <- subset(aggr_signac_unfiltered, subset = fine.cell.type %in% "immature neuron")



cell_types <- c("SEMA3E+ glut", "NPC", "forebrain", "immature neuron", 
                "subcerebral immature neuron")
sep_by_cell_type <- c(SEMA3Ep_glut, NPC, forebrain, immature_neuron, subcerebral)
time_points <- c("1hr", "6hr")

for (i in 1:length(cell_types)){
  ct <- cell_types[i]
  print(ct)
  obj <- sep_by_cell_type[[i]]
  Idents(obj) <- "time.ident"
  for (t in time_points) {
    print(paste("time:", t))
    peaks <- FindMarkers(
      object = obj,
      ident.1 = t,
      ident.2 = "0hr",
      min.pct = 0, 
      logfc.threshold = 0,
      test.use = 'LR',
      latent.vars = 'peak_region_fragments'
    )
    file_name = paste0(ct, str_sub(t, end = -3L), "v0hr_peaks.txt")
    print(file_name)
    write.table(peaks, 
                file = file_name,
                quote = F,
                sep = "\t",
                row.names = T,
                col.names = T)
  }
}

