# Chuxuan Li 28 Sept 2021
# look at the counts per cell

# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)
library(future)
library(viridis)

combined_gex <- readRDS("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/DESeq2/combined_g2_g8_USE_THIS.rds")

features <- c("BDNF", "PNOC", "GAPDH", "FOS")

# new.cluster.ids <- c("unidentified", "GABAergic", "Glutamatergic", "Glutamatergic", "GABAergic", 
#                      "Glutamatergic", "GABAergic", "Glutamatergic", "Glutamatergic", "Glutamatergic",
#                      "unidentified", "GABAergic", "Glutamatergic", "Glutamatergic", "unidentified",
#                      "unidentified", "Glutamatergic", "Glutamatergic", "Glutamatergic", "unidentified",
#                      "unidentified", "unidentified", "Glutamatergic", "unidentified", "unidentified", "unidentified")
# names(new.cluster.ids) <- levels(combined_gex)
# combined_gex_renamed <- RenameIdents(combined_gex, new.cluster.ids)

# DotPlot(combined_gex, 
#         features = features,
#         cols = "Purples",
#         group.by = "cell.type.ident") + RotatedAxis() 
# 
# DotPlot(combined_gex, 
#         features = features,
#         cols = "Purples",
#         group.by = "time.ident") + RotatedAxis() 
# 
# DotPlot(combined, 
#         features = features,
#         cols = "Purples",
#         group.by = "cell.line.ident") + RotatedAxis() 

combined_clean <- subset(combined_gex, subset = cell.type.ident %in% c("GABA", "Glut"))
# DotPlot(combined_clean, 
#         features = features,
#         cols = "Purples",
#         group.by = "cell.type.ident") + RotatedAxis() 
# 
# DotPlot(combined_clean, 
#         features = features,
#         cols = "Purples",
#         group.by = "time.ident") + RotatedAxis() 

combined_clean <- subset(combined_gex, subset = cell.line.ident %in% c("CD_54", "CD_27", "CD_26", "CD_25", "CD_08"))
combined_clean <- subset(combined_clean, subset = cell.type.ident %in% c("GABA", "Glut"))

DPnew(combined_clean, 
      features = features,
      cols = inferno(100),
      group.by = "cell.type.ident",
      split.by = "time.ident") + RotatedAxis() 

DPnew(combined_clean, 
        features = features,
        cols = "Reds",
        group.by = "cell.type.ident",
        split.by = "time.ident") + RotatedAxis() 

CD_54_clean <- subset(combined_clean, subset = cell.line.ident %in% "CD_54")
CD_27_clean <- subset(combined_clean, subset = cell.line.ident %in% "CD_27")
CD_26_clean <- subset(combined_clean, subset = cell.line.ident %in% "CD_26")
CD_25_clean <- subset(combined_clean, subset = cell.line.ident %in% "CD_25")
CD_08_clean <- subset(combined_clean, subset = cell.line.ident %in% "CD_08")

DPnew(CD_54_clean, 
      features = features,
      cols = "Reds",
      group.by = "cell.type.ident",
      split.by = "time.ident") + 
  RotatedAxis() +
  ggtitle("CD_54")

DPnew(CD_27_clean, 
      features = features,
      cols = "Reds",
      group.by = "cell.type.ident",
      split.by = "time.ident") + 
  RotatedAxis() +
  ggtitle("CD_27")

DPnew(CD_26_clean, 
      features = features,
      cols = "Reds",
      group.by = "cell.type.ident",
      split.by = "time.ident") + 
  RotatedAxis() +
  ggtitle("CD_26")

DPnew(CD_25_clean, 
      features = features,
      cols = "Reds",
      group.by = "cell.type.ident",
      split.by = "time.ident") + 
  RotatedAxis() +
  ggtitle("CD_25")

DPnew(CD_08_clean, 
      features = features,
      cols = "Reds",
      group.by = "cell.type.ident",
      split.by = "time.ident") + 
  RotatedAxis() +
  ggtitle("CD_08")
