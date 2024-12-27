# Siwei Zhang 06 Jan 2024
# plot lib 018 for for identity check
# make a subset object

# Get all libraries from 018 to 030, normalize them separately, then integrate

# init ####
{
  
  library(Seurat)
  library(EnsDb.Hsapiens.v86)
  library(patchwork)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(RColorBrewer)
  library(dplyr)
  library(viridis)
  library(graphics)
  library(readxl)
  library(readr)
  
  library(future)
}

plan("multisession", workers = 12)
options(expressions = 20000)
options(future.globals.maxSize = 161061273600)

load("~/NVME/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/018-030_RNA_integrated_labeled_with_harmony.RData")


# assign active ident
Idents(integrated_labeled) <- "lib.ident"
DefaultAssay(integrated_labeled) <- "integrated"


DimPlot(integrated_labeled,
        cells.highlight = (integrated_labeled$lib.ident %in% c("2-0",
                                                    "2-1",
                                                    "2-6",
                                                    "8-0",
                                                    "8-1",
                                                    "8-6")))

DimPlot(integrated_labeled,
        cells.highlight = colnames(integrated_labeled)[integrated_labeled$lib.ident %in% c("2-0",
                                                                                           "2-1",
                                                                                           "2-6",
                                                                                           "8-0",
                                                                                           "8-1",
                                                                                           "8-6")],
        cols.highlight = "darkred") +
  ggtitle("By Harmony")

DefaultAssay(integrated_labeled) <- "SCT"
DimPlot(integrated_labeled,
        cells.highlight = colnames(integrated_labeled)[integrated_labeled$lib.ident %in% c("2-0",
                                                                                           "2-1",
                                                                                           "2-6",
                                                                                           "8-0",
                                                                                           "8-1",
                                                                                           "8-6")],
        cols.highlight = "darkred") +
  ggtitle("By SCTran")
