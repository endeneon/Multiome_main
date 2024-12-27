# Chuxuan Li 01/05/2022
# Using the Signac-processed 10x-aggregated ATACseq dataset (new_peak_set_after_motif.RData),
#subset the data into only 4 cell types (or 3, combining the 2 gluts), then
#call findmarkers() to find time-point sensitive peaks for each cell type
#which is then used to do enrichment analysis and plot dotplot


# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)

library(readr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(graphics)
library(stringr)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)


set.seed(1105)


#unique(obj_complete$broad.cell.type)
#unique(obj_complete$fine.cell.type)

# four_celltype_obj <- subset(obj_complete, subset = fine.cell.type %in% 
#                               c("GABA", "NPC", "NEFM+/CUX2- glut", "NEFM-/CUX2+ glut"))
# three_celltype_obj <- subset(obj_complete, subset = broad.cell.type %in% 
#                               c("GABA", "NPC", "glut"))
# four_celltype_obj$fine.cell.type[four_celltype_obj$fine.cell.type == "NEFM+/CUX2- glut"] <- "NpCm_glut"
# four_celltype_obj$fine.cell.type[four_celltype_obj$fine.cell.type == "NEFM-/CUX2+ glut"] <- "NmCp_glut"

#unique(three_celltype_obj$broad.cell.type)

# save("four_celltype_obj", file = "four_celltype_obj.RData")
# save("three_celltype_obj", "three_celltype_obj.RData")
setwd("/home/cli/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes")

load("three_celltype_obj.RData")
load("four_celltype_obj.RData")

# first deal with the object with 4 cell types ####
Idents(four_celltype_obj) <- four_celltype_obj$fine.cell.type


for (i in c("GABA", "NPC")) {
  print(i)
  print("1v0")
  # 1hr vs 0hr
  df <- FindMarkers(object = four_celltype_obj,
                        ident.1 = "1hr",
                        group.by = 'time.ident',
                        subset.ident = i,
                        ident.2 = "0hr",
                        min.pct = 0.0,
                        logfc.threshold = 0.0)
  df$Gene_Symbol <- rownames(df)
  file_name <- paste0("./Output/four_cell_types_da_peaks/", 
                      i,
                      "_1v0hr.csv")
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = ",", 
              row.names = F, 
              col.names = F)
  
  print("6v0")
  # 6hr vs 0hr
  df <- FindMarkers(object = four_celltype_obj,
                    ident.1 = "6hr",
                    group.by = 'time.ident',
                    subset.ident = i,
                    ident.2 = "0hr",
                    min.pct = 0.0,
                    logfc.threshold = 0.0)
  df$Gene_Symbol <- rownames(df)
  file_name <- paste0("./Output/four_cell_types_da_peaks/", 
                      i,
                      "_6v0hr.csv"
  )
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = ",", 
              row.names = F, 
              col.names = F)
}

print("finished with four_celltypes")

# three cell types ####
Idents(three_celltype_obj) <- three_celltype_obj$broad.cell.type

for (i in unique(three_celltype_obj$broad.cell.type)) {
  print(i)
  print("1v0")
  # 1hr vs 0hr
  df <- FindMarkers(object = three_celltype_obj,
                    ident.1 = "1hr",
                    group.by = 'time.ident',
                    subset.ident = i,
                    ident.2 = "0hr",
                    min.pct = 0.0,
                    logfc.threshold = 0.0)
  df$peak <- rownames(df)
  file_name <- paste0("./Output/three_cell_types_da_peaks/", 
                      i,
                      "_1v0hr.csv"
  )
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = ",", 
              row.names = F, 
              col.names = F)
  print("6v0")
  # 6hr vs 0hr
  df <- FindMarkers(object = three_celltype_obj,
                    ident.1 = "6hr",
                    group.by = 'time.ident',
                    subset.ident = i,
                    ident.2 = "0hr",
                    min.pct = 0.0,
                    logfc.threshold = 0.0)
  df$peak <- rownames(df)
  file_name <- paste0("./Output/three_cell_types_da_peaks/", 
                      i,
                      "_6v0hr.csv"
  )
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = ",", 
              row.names = F, 
              col.names = F)
}

print("finished with three_celltypes")


