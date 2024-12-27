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

set.seed(1886)

# using LR test
da_peaks_GABA_6v1hr_LR <- FindMarkers(
  object = obj_complete,
  ident.1 = "GABA_6hr", 
  ident.2 = "GABA_1hr",
  test.use = 'LR', 
  latent.vars = "in_peaks_count",
  logfc.threshold = 0,
  min.pct = 0)

da_peaks_NmCp_glut_6v1hr_LR <- FindMarkers(
  object = obj_complete,
  ident.1 = "NmCp_glut_6hr", 
  ident.2 = "NmCp_glut_1hr",
  test.use = 'LR', 
  latent.vars = "in_peaks_count",
  logfc.threshold = 0,
  min.pct = 0)

da_peaks_NpCm_glut_6v1hr_LR <- FindMarkers(
  object = obj_complete,
  ident.1 = "NpCm_glut_6hr", 
  ident.2 = "NpCm_glut_1hr",
  test.use = 'LR', 
  latent.vars = "in_peaks_count",
  logfc.threshold = 0,
  min.pct = 0)


df_GABA_LR <- data.frame(p_val = da_peaks_GABA_6v1hr_LR$p_val,
                         avg_log2FC = da_peaks_GABA_6v1hr_LR$avg_log2FC,
                         pct.1 = da_peaks_GABA_6v1hr_LR$pct.1,
                         pct.2 = da_peaks_GABA_6v1hr_LR$pct.2,
                         p_val_adj = da_peaks_GABA_6v1hr_LR$p_val_adj, 
                         peak <- rownames(da_peaks_GABA_6v1hr_LR),
                         stringsAsFactors = T)

df_NmCp_LR <- data.frame(p_val = da_peaks_NmCp_glut_6v1hr_LR$p_val,
                         avg_log2FC = da_peaks_NmCp_glut_6v1hr_LR$avg_log2FC,
                         pct.1 = da_peaks_NmCp_glut_6v1hr_LR$pct.1,
                         pct.2 = da_peaks_NmCp_glut_6v1hr_LR$pct.2,
                         p_val_adj = da_peaks_NmCp_glut_6v1hr_LR$p_val_adj, 
                         peak <- rownames(da_peaks_NmCp_glut_6v1hr_LR),
                         stringsAsFactors = T)
df_NpCm_LR <- data.frame(p_val = da_peaks_NpCm_glut_6v1hr_LR$p_val,
                         avg_log2FC = da_peaks_NpCm_glut_6v1hr_LR$avg_log2FC,
                         pct.1 = da_peaks_NpCm_glut_6v1hr_LR$pct.1,
                         pct.2 = da_peaks_NpCm_glut_6v1hr_LR$pct.2,
                         p_val_adj = da_peaks_NpCm_glut_6v1hr_LR$p_val_adj, 
                         peak <- rownames(da_peaks_NpCm_glut_6v1hr_LR),
                         stringsAsFactors = T)

write.table(df_NmCp_LR, 
            file = "NmCp_glut_LR_6V1_all_peaks.csv",
            quote = F, 
            sep = ",", 
            row.names = F, 
            col.names = T)
write.table(df_NpCm_LR, 
            file = "NpCm_glut_LR_6V1_all_peaks.csv",
            quote = F, 
            sep = ",", 
            row.names = F, 
            col.names = T)
write.table(df_GABA_LR, 
            file = "GABA_LR_6V1_all_peaks.csv",
            quote = F, 
            sep = ",", 
            row.names = F, 
            col.names = T)

save.image(file = "da_peaks_complete_6v1hr.RData")
