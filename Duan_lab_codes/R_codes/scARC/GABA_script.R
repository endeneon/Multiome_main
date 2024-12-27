# Chuxuan Li 06/09/2022
# 5-18-20 line integrated data splitted by cell type - this is GABA
# use Monocle3 to calculate differential gene expression. Run by bash

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

unique(integrated_labeled$cell.line.ident)

unique(integrated_labeled$cell.type)

# GABA #####
GABA <- subset(integrated_labeled, cell.type %in% c("GABA", "SST_pos_GABA", "SEMA3E_pos_GABA"))
pct_exp <- rowSums(GABA@assays$RNA@counts, na.rm = T) / ncol(GABA@assays$RNA@counts)
expressed_genes = names(pct_exp[pct_exp >= 0.01])

GABA@assays$integrated@counts <- 
  GABA@assays$RNA@counts[GABA@assays$RNA@counts@Dimnames[[1]] %in% expressed_genes, ]
GABA@assays$integrated@data <- 
  GABA@assays$RNA@data[GABA@assays$RNA@data@Dimnames[[1]] %in% expressed_genes, ]

DefaultAssay(GABA) <- "integrated"
GABA[["RNA"]] <- NULL
GABA[["SCT"]] <- NULL
GABA_cds <- SeuratWrappers::as.cell_data_set(GABA)
print(paste0("number of genes left: ", length(expressed_genes))) #34268

# need to calculate size factor manually
GABA_cds <- estimate_size_factors(GABA_cds)
# add gene names manually
GABA_cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(GABA)
print("starting fit_models()......")
# fit linear model
# returns a tibble that contains a row for each gene. 
# column contains generalized linear model integrated_labeledects, model formula can be any
#column in the colData table
GABA_gene_fits <- fit_models(cds = GABA_cds, 
                             model_formula_str = "~time.ident + orig.ident", 
                             verbose = T)
print("finished fit_models()!")
# extract coeff table
GABA_fit_coefs <- coefficient_table(GABA_gene_fits)
unique(GABA_fit_coefs$term)
print("saving......")
save(GABA_fit_coefs, file = "GABA_fit_coefs.RData")

