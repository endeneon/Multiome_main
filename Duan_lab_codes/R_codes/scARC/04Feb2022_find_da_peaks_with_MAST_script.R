# Chuxuan Li 12/17/2021
# calculate differentially accessible peaks using findmarkers() for 0v6hr, 1v6hr,
#for GABA and 2 glut groups

# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)
library(future)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)

set.seed(1911)
# set threads and parallelization
plan("multisession", workers = 1)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 11474836480)


# load data
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/four_celltype_obj.RData")

# assign idents
Idents(four_celltype_obj) <- four_celltype_obj$fine.cell.type

# run FindMarkers()
for (i in unique(four_celltype_obj$fine.cell.type)) {
  print(i)
  print("1v0")
  # 1hr vs 0hr
  df <- FindMarkers(object = four_celltype_obj,
                    ident.1 = "1hr",
                    group.by = 'time.ident',
                    subset.ident = i,
                    ident.2 = "0hr",
                    min.pct = 0.0,
                    logfc.threshold = 0.0,
                    test.use = "MAST")
  df$peaks <- rownames(df)
  file_name <- paste0("./Output/four_cell_types_use_MAST/", 
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
                    logfc.threshold = 0.0, 
                    test.use = "MAST")
  df$peaks <- rownames(df)
  file_name <- paste0("./Output/four_cell_types_use_MAST/", 
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


