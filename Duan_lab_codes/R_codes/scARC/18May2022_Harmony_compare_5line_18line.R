# Chuxuan Li 05/18/2022
# Check the similarity between 18-line and 5-line RNASeq data using Harmony
#embedding values

# init ####
library(Seurat)
library(Signac)
library(edgeR)
library(harmony)

library(FactoMineR)
library(factoextra)

library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(reshape2)


load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v2_combine_5line_18line/5line_18line_combined_labeled_nfeat5000_obj.RData")
unique(integrated_labeled$batch.ident)
harmony_embeddings <- harmony::HarmonyMatrix(
  V, meta_data, 'dataset', do_pca = FALSE, verbose=FALSE
)