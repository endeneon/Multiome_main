# Chuxuan Li 04/11/2022
# use monocle3 to analyze single-cell differential expression of 18-line data

# init ####
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(SingleCellExperiment)

library(dplyr)
library(tidyr)
library(stringr)

library(ggplot2)
library(RColorBrewer)

library(future)

plan("multisession", workers = 1)

load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v4_combine_5_18_20line/5_18_20line_combined_labeled_nfeat5000_obj.RData")
obj <- subset(integrated_labeled, cell.type != "unknown")
rm(integrated_labeled)
# assign affected status
load("../covariates_pooled_together_for_all_linextime_5+18+20.RData")
celllines <- sort(unique(obj$cell.line.ident))
obj$aff <- "control"
for (i in 1:length(celllines)) {
  l <- celllines[i]
  print(l)
  aff <- unique(covar_table$disease[covar_table$cell_line == l])
  obj$aff[obj$cell.line.ident == l] <- aff
}
unique(obj$aff)
sum(obj$aff == "control")
sum(obj$aff == "case")


# same function applied to each cell type ####

monoclePipeline <- function(subtypeobj, typename) {
  print(typename)
  DefaultAssay(subtypeobj) <- "SCT"
  subtypeobj@assays$RNA <- NULL
  subtypeobj@assays$integrated <- NULL
  # filter genes
  pct_exp <- rowSums(subtypeobj@assays$SCT@data == 0, na.rm = T) / ncol(subtypeobj@assays$RNA@counts)
  expressed_genes = names(pct_exp[pct_exp >= 0.01])
  subtypeobj@assays$SCT@data <- 
    subtypeobj@assays$SCT@data[subtypeobj@assays$SCT@data@Dimnames[[1]] %in% expressed_genes, ]
  print(paste0("number of genes left: ", length(expressed_genes)))
  
  times <- c("0hr", "1hr", "6hr")
  for (i in 1:length(times)) {
    print(times[i])
    subtimeobj <- subset(subtypeobj, time.ident == times[i])
    cds <- SeuratWrappers::as.cell_data_set(subtypeobj)
    
    cds <- estimate_size_factors(cds)
    cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(subtypeobj)
    # fit linear model
    fits <- fit_models(cds = cds, 
                       model_formula_str = "~aff + cell.line.ident", 
                       verbose = T)
    # extract coeff table
    fit_coefs <- coefficient_table(fits)
    unique(fit_coefs$term)
    filename <- paste0(typename, "monocle_fit_coefs_", times[i], ".RData")
    print("saving......")
    save(fit_coefs, file = filename)
  }
}

GABA <- subset(obj, cell.type %in% c("GABA", "SST_pos_GABA", "SEMA3E_pos_GABA"))
nmglut <- subset(obj, cell.type == "NEFM_neg_glut")
npglut <- subset(obj, cell.type == "NEFM_pos_glut")

monoclePipeline(GABA, "GABA")
monoclePipeline(nmglut, "nmglut")
monoclePipeline(npglut, "npglut")