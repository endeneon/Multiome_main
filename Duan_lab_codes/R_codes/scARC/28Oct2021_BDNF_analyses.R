# Chuxuan Li 10/28/2021
# Check BDNF gene expression changes, any BDNF related analyses

library(Seurat)
library(future)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)
library(patchwork)
library(cowplot)
library(gplots)

# set threads and parallelization
plan("multisession", workers = 6)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

# Plot gene expression changes for BDNF across two types of glut
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/integrated_labeled_updated.RData")
# assign bigger cell types to the RNAseq object first
integrated_renamed_1$cell.type <- Idents(integrated_renamed_1)
integrated_renamed_1$broad.cell.type <- "other"
integrated_renamed_1$broad.cell.type[integrated_renamed_1$cell.type %in% 
                                         c("GAD1+/GAD2+ GABA", "GAD1+ GABA (less mature)", "GAD1+ GABA")] <- "GABA"
integrated_renamed_1$broad.cell.type[integrated_renamed_1$cell.type %in% 
                                         c("NEFM+/CUX2- glut", "NEFM-/CUX2+ glut",
                                           "NEFM+/CUX2-, SST+ glut", "NEFM-/CUX2+, ADCYAP1+ glut",
                                           "NEFM+/CUX2-, ADCYAP1+ glut", "NEFM-/CUX2+, TCERG1L+ glut",
                                           "NEFM+/CUX2+, SST+ glut")] <- "glut"
integrated_renamed_1$broad.cell.type[integrated_renamed_1$cell.type %in% 
                                         c("NPC", "immature neuron", "MAP2+ NPC", "subcerebral immature neuron")] <- "NPC"
integrated_renamed_1$fine.cell.type <- "other"
integrated_renamed_1$fine.cell.type[integrated_renamed_1$cell.type %in% 
                                        c("GAD1+/GAD2+ GABA", "GAD1+ GABA (less mature)", "GAD1+ GABA")] <- "GABA"
integrated_renamed_1$fine.cell.type[integrated_renamed_1$cell.type %in% 
                                        c("NEFM+/CUX2- glut", "NEFM+/CUX2-, SST+ glut", 
                                          "NEFM+/CUX2-, ADCYAP1+ glut")] <- "NEFM+/CUX2- glut"
integrated_renamed_1$fine.cell.type[integrated_renamed_1$cell.type %in% 
                                        c("NEFM-/CUX2+, glut",
                                          "NEFM-/CUX2+, ADCYAP1+ glut",
                                          "NEFM-/CUX2+, TCERG1L+ glut"
                                        )] <- "NEFM-/CUX2+ glut"
integrated_renamed_1$fine.cell.type[integrated_renamed_1$cell.type %in% 
                                        c("NPC", "immature neuron", "MAP2+ NPC", "subcerebral immature neuron")] <- "NPC"

NpCm_glut_RNA <- subset(integrated_renamed_1, subset = fine.cell.type %in% "NEFM+/CUX2- glut")
NmCp_glut_RNA <- subset(integrated_renamed_1, subset = fine.cell.type %in% "NEFM-/CUX2+ glut")
FeaturePlot(NpCm_glut_RNA, 
            features = "BDNF",
            split.by = "time.ident",
            cols = c("lightgrey", "slateblue4"))
FeaturePlot(NmCp_glut_RNA, 
            features = "BDNF",
            split.by = "time.ident",
            pt.size = 0.2,
            cols = c("lightgrey", "slateblue4"))
