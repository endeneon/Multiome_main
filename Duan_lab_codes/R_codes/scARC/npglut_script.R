# Chuxuan Li 06/09/2022
# 5-18-20 line integrated data splitted by cell type - this is npglut
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

# npglut #####
npglut <- subset(integrated_labeled, cell.type == "NEFM_pos_glut")

npglut_cds <- SeuratWrappers::as.cell_data_set(npglut)

# need to calculate size factor manually
npglut_cds <- estimate_size_factors(npglut_cds)
# add gene names manually
npglut_cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(npglut)
print("starting fit_models()......")
# fit linear model
# returns a tibble that contains a row for each gene. 
# column contains generalized linear model integrated_labeledects, model formula can be any
#column in the colData table
npglut_gene_fits <- fit_models(cds = npglut_cds, 
                               model_formula_str = "~time.ident + orig.ident", 
                               verbose = T)
print("finished fit_models()!")
# extract coeff table
npglut_fit_coefs <- coefficient_table(npglut_gene_fits)
unique(npglut_fit_coefs$term)
print("saving......")
save(npglut_fit_coefs, file = "npglut_fit_coefs.RData")


