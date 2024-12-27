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

load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")

unique(integrated_labeled$cell.line.ident)
obj <- subset(integrated_labeled, cell.line.ident != "unmatched")
rm(integrated_labeled)

unique(obj$cell.type)

# GABA #####
GABA <- subset(obj, cell.type %in% c("GABA", "SST_pos_GABA", "SEMA3E_pos_GABA"))
pct_exp <- rowSums(GABA@assays$RNA@counts == 0, na.rm = T) / ncol(GABA@assays$RNA@counts)
expressed_genes = names(pct_exp[pct_exp >= 0.01])
# hist(rowSums(GABA@assays$RNA@counts[GABA@assays$RNA@counts@Dimnames[[1]] %in% 
#                                       expressed_genes, ]), breaks = 10000, 
#      xlim = c(0, 10))
# sum(rowSums(GABA@assays$RNA@counts[GABA@assays$RNA@counts@Dimnames[[1]] %in% 
#                                  expressed_genes, ]) == 0)
GABA@assays$RNA@counts <- 
  GABA@assays$RNA@counts[GABA@assays$RNA@counts@Dimnames[[1]] %in% expressed_genes, ]
GABA@assays$RNA@data <- 
  GABA@assays$RNA@data[GABA@assays$RNA@data@Dimnames[[1]] %in% expressed_genes, ]

# DefaultAssay(GABA) <- "integrated"
# GABA[["RNA"]] <- NULL
GABA_cds <- SeuratWrappers::as.cell_data_set(GABA)
print(paste0("number of genes left: ", length(expressed_genes))) #34268
# expr_matrix  <- GetAssayData(GABA, assay = 'integrated', slot = 'data')
# # create cell dataset, don't worry about the warning
# GABA_cds <- monocle3::new_cell_data_set(
#   expr_matrix,
#   cell_metadata = GABA@meta.data
# )
# GABA_cds <- preprocess_cds(GABA_cds, method = "PCA", norm_method = "log")
# GABA_cds <- reduce_dimension(GABA_cds, reduction_method = "UMAP", verbose = T)
# GABA_cds <- cluster_cells(GABA_cds, reduction_method = 'UMAP')
# GABA_cds <- learn_graph(GABA_cds, use_partition = TRUE)
# 
# # use the helper function in the monocle3 documentation
# # https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
# # a helper function to identify the root principal points:
# get_earliest_principal_node <- function(cds, time_bin = "0hr"){
#   cell_ids <- which(colData(cds)[, "time.ident"] == time_bin)
#   
#   closest_vertex <-
#     cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   root_pr_nodes <-
#     igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
#                                                               (which.max(table(closest_vertex[cell_ids,]))))]
#   
#   root_pr_nodes
# }
# principal_node <- get_earliest_principal_node(GABA_cds)
# 
# # order cells in pseudotime 
# GABA_cds <- order_cells(GABA_cds, root_pr_nodes = principal_node)

# need to calculate size factor manually
GABA_cds <- estimate_size_factors(GABA_cds)
# add gene names manually
GABA_cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(GABA)

# fit linear model
# returns a tibble that contains a row for each gene. 
# column contains generalized linear model objects, model formula can be any
#column in the colData table
GABA_gene_fits <- fit_models(cds = GABA_cds, 
                             model_formula_str = "~time.ident + orig.ident", 
                             verbose = T)
# extract coeff table
GABA_fit_coefs <- coefficient_table(GABA_gene_fits)
unique(GABA_fit_coefs$term)
save(GABA_fit_coefs, file = "GABA_fit_coefs.RData")

# nmglut #####
nmglut <- subset(obj, cell.type == "NEFM_neg_glut")

nmglut_cds <- SeuratWrappers::as.cell_data_set(nmglut)

# need to calculate size factor manually
nmglut_cds <- estimate_size_factors(nmglut_cds)
# add gene names manually
nmglut_cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(nmglut)

# fit linear model
# returns a tibble that contains a row for each gene. 
# column contains generalized linear model objects, model formula can be any
#column in the colData table
nmglut_gene_fits <- fit_models(cds = nmglut_cds, 
                               model_formula_str = "~time.ident + orig.ident", 
                               verbose = T)
# extract coeff table
nmglut_fit_coefs <- coefficient_table(nmglut_gene_fits)
unique(nmglut_fit_coefs$term)
save(nmglut_fit_coefs, file = "nmglut_fit_coefs.RData")

# npglut #####
npglut <- subset(obj, cell.type == "NEFM_pos_glut")

npglut_cds <- SeuratWrappers::as.cell_data_set(npglut)

# need to calculate size factor manually
npglut_cds <- estimate_size_factors(npglut_cds)
# add gene names manually
npglut_cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(npglut)

# fit linear model
# returns a tibble that contains a row for each gene. 
# column contains generalized linear model objects, model formula can be any
#column in the colData table
npglut_gene_fits <- fit_models(cds = npglut_cds, 
                               model_formula_str = "~time.ident + orig.ident", 
                               verbose = T)
# extract coeff table
npglut_fit_coefs <- coefficient_table(npglut_gene_fits)
unique(npglut_fit_coefs$term)
save(npglut_fit_coefs, file = "npglut_fit_coefs.RData")


