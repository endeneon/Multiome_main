# Chuxuan Li 02/08/2023
# Use data mapped to rat + rabbit, separate rat and rabbit after removing cells
#mapped to human barcodes

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
load("../029_RNA_integrated_labeled.RData")
setwd("/nvmefs/scARC_Duan_029_Oc2_Mm10/")
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
setwd("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/oc2_mm10_analysis/")
save(objlist, file = "mapped_to_oc2_mm10_raw_data_list.RData")

# remove human cells ####
libs <- str_extract(string = h5list, pattern = "[0-9]+-[0-6]")
integrated_labeled <- subset(integrated_labeled, orig.ident %in% libs)
cell_count_by_species <- matrix(nrow = length(objlist), ncol = 6, 
                         dimnames = list(libs, c("total number of cells",
                                           "human (mapped to demuxlet barcodes)", 
                                           "nonhuman (not mapped to barcodes)",
                                           "rabbit cells",
                                           "mouse cells",
                                           "other cells (unknown identity)")))
for (i in 1:length(objlist)) {
  lib <- unique(objlist[[i]]$orig.ident)
  sub_intobj <- subset(integrated_labeled, orig.ident == lib)
  human_bc <- str_remove(sub_intobj@assays$RNA@counts@Dimnames[[2]], "_[0-9]+$")
  obj <- objlist[[i]]
  obj$human.ident <- "nonhuman"
  obj$human.ident[obj@assays$RNA@counts@Dimnames[[2]] %in% human_bc] <- "human"
  cell_count_by_species[i, 1] <- length(obj$human.ident)
  cell_count_by_species[i, 2] <- sum(obj$human.ident == "human")
  cell_count_by_species[i, 3] <- sum(obj$human.ident == "nonhuman")
  obj <- subset(obj, human.ident == "nonhuman")
  objlist[[i]] <- obj
}

save(objlist, file = "removed_human_cells_objlist.RData")

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
save(scaled_lst, file = "scaled_clustered_list.RData")
rm(transformed_lst)

# look at expression ####
for (i in 1:length(scaled_lst)) {
  print(DefaultAssay(scaled_lst[[i]]))
  png(paste0("./dimplot_by_cluster/", 
             str_extract(h5list[i], "[0-9]+-[0-6]"), "_dimplot_by_cluster.png"))
  p <- DimPlot(scaled_lst[[i]],
               label = T,
               group.by = "seurat_clusters") +
    NoLegend() +
    ggtitle(str_extract(h5list[i], "[0-9]+-[0-6]")) +
    theme(text = element_text(size = 12))
  print(p)
  dev.off()
  
  if ("Gfap" %in% rownames(scaled_lst[[i]])) {
    png(paste0("./vlnplot_mouse/", str_extract(h5list[i], "[0-9]+-[0-6]"), "_vln.png"),
        width = 1200, height = 400)
    p <- VlnPlot(scaled_lst[[i]], c("Gfap", "S100b")) 
    print(p)
    dev.off()
 
    png(paste0("./featplot_mouse/", str_extract(h5list[i], "[0-9]+-[0-6]"), "_featplot.png"),
        width = 1000, height = 400)
    p <- FeaturePlot(scaled_lst[[i]], features = c("Gfap", "S100b")) 
    print(p)
    dev.off()
  } 
  
  if ("ENSOCUG00000002862" %in% rownames(scaled_lst[[i]])) {
    png(paste0("./vlnplot_rabbit/",
               str_extract(h5list[i], "[0-9]+-[0-6]"), "_rabbit_vln.png"), width = 1200, height = 400)
    p <- VlnPlot(scaled_lst[[i]], c("ENSOCUG00000002862", "ENSOCUG00000013326")) +
      ggtitle(paste0(str_extract(h5list[i], "[0-9]+-[0-6]"), 
                     " rabbit ribosomal protein + translation initiation factor expression by cluster"))
    print(p)
    dev.off()
    
    png(paste0("./featplot_rabbit/", str_extract(h5list[i], "[0-9]+-[0-6]"), "_featplot.png"),
        width = 1000, height = 400)
    p <- FeaturePlot(scaled_lst[[i]], features = c("ENSOCUG00000002862", "ENSOCUG00000013326")) 
    print(p)
    dev.off()
  } 
  
}

test <- scaled_lst[[1]]
mouse_gene_2 <- sort(rowSums(test)[test$seurat_clusters == 2], decreasing = T)
mouse_gene_4 <- sort(rowSums(test)[test$seurat_clusters == 4], decreasing = T)
mouse_gene_6 <- sort(rowSums(test)[test$seurat_clusters == 6], decreasing = T)
mouse_gene_7 <- sort(rowSums(test)[test$seurat_clusters == 7], decreasing = T)

rabbit_clusters <- list(c(0, 1, 3, 5), c(0, 1, 2), c(0, 1, 2, 4), #20088
                     c(0, 1, 2, 6), c(0, 3, 5), c(0, 1, 2, 4)) #70179
mouse_clusters <- list(c(2, 4, 6, 7), c(3, 4, 5, 6), c(3, 5, 6, 7), #20088
                        c(3, 4, 5), c(1, 2), c(3, 5 ,6)) #70179
for (i in 1:length(scaled_lst)) {
  rabbit <- rabbit_clusters[[i]]
  mouse <- mouse_clusters[[i]]
  obj <- scaled_lst[[i]]
  obj$species.ident <- "unknown"
  obj$species.ident[scaled_lst[[i]]$seurat_clusters %in% rabbit] <- "rabbit"
  obj$species.ident[scaled_lst[[i]]$seurat_clusters %in% mouse] <- "mouse"
  cell_count_by_species[i, 4] <- sum(obj$species.ident == "rabbit")
  cell_count_by_species[i, 5] <- sum(obj$species.ident == "mouse")
  cell_count_by_species[i, 6] <- sum(obj$species.ident == "unknown")
  scaled_lst[[i]] <- obj
}

save(scaled_lst, file = "labeled_by_species_obj_list.RData")
write.table(cell_count_by_species, 
            file = "cell_count_by_species_separated_by_lib_analysis.csv", 
            quote = F, sep = ",", row.names = T, col.names = T)
save(cell_count_by_species, file = "cell_count_by_species_separated_by_lib_analysis.RData")


### extract mouse barcodes
load("~/NVME/scARC_Duan_018/Duan_project_029_RNA/oc2_mm10_analysis/removed_human_cells_objlist.RData")
objlist[[1]]$
  load("~/NVME/scARC_Duan_018/Duan_project_029_RNA/oc2_mm10_analysis/rabbit_mouse_labeled_cells_20230420.RData")
scaled_lst[[1]]$mouse_or_rabbit
