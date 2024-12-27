# Chuxuan 03/09/2022
# use monocle on NormalizeData() normalized 5-line data with only four celltypes

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

load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/27Jan2022_5line_added_cellline_idents.RData")

unique(RNAseq_integrated_labeled$cell.line.ident)
obj <- subset(RNAseq_integrated_labeled, cell.line.ident != "unmatched")
rm(RNAseq_integrated_labeled)

unique(obj$broad.cell.type)

# GABA #####
GABA <- subset(obj, broad.cell.type == "GABA")
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
print(paste0("number of genes left: ", length(expressed_genes)))
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
                        model_formula_str = "~time.ident + lib.ident", 
                        verbose = T)
# extract coeff table
GABA_fit_coefs <- coefficient_table(GABA_gene_fits)
unique(GABA_fit_coefs$term)


# nmglut #####
nmglut <- subset(obj, broad.cell.type == "NEFM_neg_glut")

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
                        model_formula_str = "~time.ident + lib.ident", 
                        verbose = T)
# extract coeff table
nmglut_fit_coefs <- coefficient_table(nmglut_gene_fits)
unique(nmglut_fit_coefs$term)

# npglut #####
npglut <- subset(obj, broad.cell.type == "NEFM_pos_glut")

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
                               model_formula_str = "~time.ident + lib.ident", 
                               verbose = T)
# extract coeff table
npglut_fit_coefs <- coefficient_table(npglut_gene_fits)
unique(npglut_fit_coefs$term)
save(npglut_fit_coefs, file = "npglut_fit_coefs.RData")


# NPC #####
NPC <- subset(obj, broad.cell.type == "NPC")

NPC_cds <- SeuratWrappers::as.cell_data_set(NPC)

# need to calculate size factor manually
NPC_cds <- estimate_size_factors(NPC_cds)
# add gene names manually
NPC_cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(NPC)

# fit linear model
# returns a tibble that contains a row for each gene. 
# column contains generalized linear model objects, model formula can be any
#column in the colData table
NPC_gene_fits <- fit_models(cds = NPC_cds, 
                               model_formula_str = "~time.ident + lib.ident", 
                               verbose = T)
# extract coeff table
npglut_fit_coefs <- coefficient_table(npglut_gene_fits)
unique(npglut_fit_coefs$term)
save(npglut_fit_coefs, file = "npglut_fit_coefs.RData")



save("nmglut_fit_coefs", file = "nmglut_fit_coefs_complete_df.RData")
