# Chuxuan Li 03/14/2022
# script passing NPC subobject to monocle fit_models()

# init ####
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(SingleCellExperiment)
library(BiocGenerics)

library(dplyr)
library(tidyr)
library(stringr)

library(ggplot2)
library(RColorBrewer)

library(future)

plan("multisession", workers = 1)
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/27Jan2022_5line_added_cellline_idents.RData")

print("loaded object")

unique(RNAseq_integrated_labeled$cell.line.ident)
obj <- subset(RNAseq_integrated_labeled, cell.line.ident != "unmatched")
rm(RNAseq_integrated_labeled)

obj <- subset(obj, broad.cell.type == "NPC")
DefaultAssay(obj) <- "RNA"
obj$line.time.ident <- "unmatched"
for (l in unique(obj$cell.line.ident)){
  for (t in unique(obj$time.ident)){
    obj$line.time.ident[obj$cell.line.ident == l & obj$time.ident == t] <- paste0(l, "_", t)
  }
}
unique(obj$line.time.ident)

pct_exp <- rowSums(obj@assays$RNA@counts > 0, na.rm = T) / ncol(obj@assays$RNA@counts)
expressed_genes = names(pct_exp[pct_exp >= 0.01])
print(paste0("number of genes left: ", length(expressed_genes)))

# subset genes
obj@assays$RNA@counts <- 
  obj@assays$RNA@counts[obj@assays$RNA@counts@Dimnames[[1]] %in% expressed_genes, ]
obj@assays$RNA@data <- 
  obj@assays$RNA@data[obj@assays$RNA@data@Dimnames[[1]] %in% expressed_genes, ]


cds <- SeuratWrappers::as.cell_data_set(obj)
fData(cds)$gene_short_name <- rownames(obj)
# need to calculate size factor manually
cds <- estimate_size_factors(cds)
# add gene names manually

print("Start fitting model")
# fit linear model
# returns a tibble that contains a row for each gene.
# column contains generalized linear model objects, model formula can be any
#column in the colData table
gene_fits <- fit_models(cds = cds,
                        model_formula_str = "~time.ident + line.time.ident",
                        verbose = T)
print("finished!")

# extract coeff table
NPC_fit_coefs <- coefficient_table(gene_fits)
save(NPC_fit_coefs, file = "NPC_fit_coefs_filtered.RData")
