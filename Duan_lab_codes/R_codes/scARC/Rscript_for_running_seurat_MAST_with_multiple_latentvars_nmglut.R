# Chuxuan Li 10/11/2022
# rscript for bash running single-cell differential expression with findmarkers()
#and MAST, plus co-cultured batch, age, sex, disease, cellular composition as latent.vars

# init ####
library(Seurat)
library(SeuratWrappers)

library(dplyr)
library(tidyr)
library(stringr)

library(ggplot2)
library(RColorBrewer)

load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000_with_covars_from_30Sep2022_single_cell_DE.RData")
nmglut <- subset(obj, cell.type == "NEFM_neg_glut")

# compute deg ####
find_time_DEG <- function(obj, min.pct, type, time){
  DefaultAssay(obj) <- "RNA"
  obj <- ScaleData(obj)
  Idents(obj) <- "time.ident"
  df <- FindMarkers(object = obj, assay = "RNA", slot = "scale.data",
                    ident.1 = time, group.by = 'time.ident', 
                    ident.2 = "0hr",
                    min.pct = min.pct, logfc.threshold = 0.0,
                    test.use = "MAST", 
                    latent.vars = c("orig.ident", "age", "sex", "aff", "cell.type.fraction"), 
                    random.seed = 42)
  df$gene_symbol <- rownames(df)
  df$cell_type <- type
  return(df)
}

#nmglut_res_1v0 <- find_time_DEG(nmglut, 0.01, "nmglut", "1hr")
nmglut_res_6v0 <- find_time_DEG(nmglut, 0.01, "nmglut", "6hr")

#save(nmglut_res_1v0, file = "seurat_MAST_w_more_latentvars_18line_DE_results_nmglut1v0.RData")
save(nmglut_res_6v0, file = "seurat_MAST_w_more_latentvars_18line_DE_results_nmglut6v0.RData")

