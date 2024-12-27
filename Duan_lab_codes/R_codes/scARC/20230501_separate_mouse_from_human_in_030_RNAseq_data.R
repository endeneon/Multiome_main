# Chuxuan Li 05/01/2023
# separate mouse/human data from 030 mapped to hybrid genome data
#by clustering individual libraries, output barcodes

# init ####
library(Seurat)
library(glmGamPoi)
library(sctransform)

library(ggplot2)
library(RColorBrewer)
library(viridis)
library(graphics)
library(ggrepel)

library(dplyr)
library(readr)
library(stringr)
library(readxl)

library(future)
# set threads and parallelization
plan("multisession", workers = 1)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

# load data ####
setwd("/data/FASTQ/Duan_Project_030/")
# note this h5 file contains both atac-seq and gex information
# Read files, separate the ATACseq data from the .h5 matrix
h5list <- list.files(path = ".",
                     pattern = "filtered_feature_bc_matrix.h5",
                     recursive = T,
                     include.dirs = F, 
                     full.names = T)
objlist <- vector(mode = "list", length = length(h5list))
names <- c()
for (i in 1:length(objlist)){
  h5file <- Read10X_h5(filename = h5list[i])
  print(str_extract(string = h5list[i],
                    pattern = "[0-9]+-[0-6]"))
  names <- c(names, as.vector(str_extract(string = h5list[i],
                                          pattern = "[0-9]+-[0-6]")))
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(string = h5list[i],
                                                  pattern = "[0-9]+-[0-6]"))
  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list[i],
                                                  pattern = "-[0-6]"), "-"), "hr")
  print(unique(obj$time.ident))
  obj$group.ident <- str_remove(str_extract(string = h5list[i],
                                            pattern = "[0-9]+-"), "-")
  # check number of genes and cells
  print(paste0(i, " number of genes: ", nrow(obj),
               ", number of cells: ", ncol(obj)))
  objlist[[i]] <- obj
}
setwd("/nvmefs/scARC_Duan_018/Duan_project_030_RNA/")
save(objlist, file = "raw_data_list.RData")


# normalization ####
transformed_lst <- vector(mode = "list", length(h5list))
for (i in 1:length(objlist)){
  print(paste0("now at ", i))
  # proceed with normalization
  obj <- PercentageFeatureSet(objlist[[i]],
                              pattern = c("^MT-"),
                              col.name = "percent.mt")
  obj <- FindVariableFeatures(obj, nfeatures = 3000)
  transformed_lst[[i]] <- SCTransform(obj,
                                      vars.to.regress = "percent.mt",
                                      method = "glmGamPoi", 
                                      return.only.var.genes = F,
                                      variable.features.n = 8000,
                                      seed.use = 2022,
                                      verbose = T)
}

# dimensional reduction and clustering ####
scaled_lst <- vector(mode = "list", length(transformed_lst))
dimreduclust <- function(obj){
  obj <- ScaleData(obj)
  obj <- RunPCA(obj,
                verbose = T,
                seed.use = 2022)
  obj <- RunUMAP(obj,
                 reduction = "pca",
                 dims = 1:30,
                 seed.use = 2022)
  obj <- FindNeighbors(obj,
                       reduction = "pca",
                       dims = 1:30)
  obj <- FindClusters(obj,
                      resolution = 0.5,
                      random.seed = 2022)
  return(obj)
  
}
for (i in 1:length(transformed_lst)){
  scaled_lst[[i]] <- dimreduclust(transformed_lst[[i]])
}

save(transformed_lst, file = "transformed_list.RData")
save(scaled_lst, file = "scaled_clustered_obj_lst.RData")
rm(transformed_lst)

# check and remove rat cells ####
cleanobj_lst <- vector(mode = "list", length = length(scaled_lst))

for (i in 1:length(scaled_lst)) {
  png(paste0("./hybrid_feature_plots/dimplot_by_cluster/", 
             str_extract(h5list[i], "[0-9]+-[0-6]"), "_dimplot_by_cluster.png"))
  p <- DimPlot(scaled_lst[[i]],
               label = T,
               group.by = "seurat_clusters") +
    NoLegend() +
    ggtitle(str_extract(h5list[i], "[0-9]+-[0-6]")) +
    theme(text = element_text(size = 12))
  print(p)
  dev.off()
  
  png(paste0("./hybrid_feature_plots/featplots/", 
             str_extract(h5list[i], "[0-9]+-[0-6]"), "_specific_genes.png"))
  p <- FeaturePlot(scaled_lst[[i]], 
                   features = c("Gfap", "S100b", "GAD1", "SLC17A6"), ncol = 2)
  print(p)
  dev.off()
  
  
  png(paste0("./hybrid_feature_plots/vlnplot_gfap/",
             str_extract(h5list[i], "[0-9]+-[0-6]"), "_Gfap_vln.png"))
  p <- VlnPlot(scaled_lst[[i]], "Gfap") +
    ggtitle(paste0(str_extract(h5list[i], "[0-9]+-[0-6]"), " Gfap expression by cluster"))
  print(p)
  dev.off()
}

# remove clusters that are rat cells, count
rat_clusters <- list(c(11), c(16, 13), c(12), 
                     c(9, 11), c(9), c(8), 
                     c(7, 12), c(11), c(9), 
                     c(7), c(11), c(11))
for (i in 1:length(scaled_lst)) {
  c <- rat_clusters[[i]]
  cat("rat cluster number: ", c)
  scaled_lst[[i]]$mouse.ident <- "human"
  scaled_lst[[i]]$mouse.ident[scaled_lst[[i]]$seurat_clusters %in% c] <- "mouse"
  human_barcodes <- colnames(scaled_lst[[i]])[scaled_lst[[i]]$mouse.ident == "human"]
  write.table(human_barcodes, file = paste0("./human_only_barcodes/", 
                                          str_extract(h5list[i], "[0-9]+-[0-6]"),
                                          "_human_barcodes.txt"), 
              sep = "\t", quote = F, row.names = F, col.names = F)
  cleanobj_lst[[i]] <- subset(scaled_lst[[i]], subset = mouse.ident != "mouse")
}

mouse_cell_count <- matrix(nrow = length(scaled_lst), ncol = 3, 
                         dimnames = list(names, 
                                         c("number of mouse astrocytes", "total number of cells",
                                           "percentage of cells that are astrocytes (%)")))
for (i in 1:length(scaled_lst)){
  mouse_cell_count[i, 1] <- (sum(scaled_lst[[i]]$mouse.ident ==  "mouse"))
  mouse_cell_count[i, 2] <- (length(scaled_lst[[i]]$mouse.ident))
  mouse_cell_count[i, 3] <- mouse_cell_count[i, 1] / mouse_cell_count[i, 2] * 100
}
write.table(mouse_cell_count, file = "mouse_cell_count_by_lib.csv", sep = ",", quote = F)

save(cleanobj_lst, file = "removed_mouse_list_030_mapped_to_hybrid.RData")

