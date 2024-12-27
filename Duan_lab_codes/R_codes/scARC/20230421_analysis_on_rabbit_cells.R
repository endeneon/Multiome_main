# Chuxuan Li 04/21/2023
# Use the mouse-rabbit separation from 20230420 to do analysis on rabbit cells ####

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

load("rabbit_mouse_labeled_cells_20230420.RData")
load("removed_human_cells_objlist.RData")

library(ensembldb)

# extract rabbit cells only ####
for (i in 1:length(objlist)) {
  if (i == 1) {
    sample_libs <- c(as.character(unique(objlist[[i]]$orig.ident)))
    query_libs <- c(as.character(unique(scaled_lst[[i]]$orig.ident)))
  } else {
    sample_libs <- c(sample_libs, as.character(unique(objlist[[i]]$orig.ident)))
    query_libs <- c(query_libs, as.character(unique(scaled_lst[[i]]$orig.ident)))
  }
}
names(objlist) <- sample_libs
for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  queryobj <- scaled_lst[[i]]
  rabbit_bc <- queryobj@assays$RNA@counts@Dimnames[[2]][queryobj$mouse_or_rabbit == "rabbit"]
  obj$mouse_or_rabbit <- "mouse"
  obj$mouse_or_rabbit[obj@assays$RNA@counts@Dimnames[[2]] %in% rabbit_bc] <- "rabbit"
  obj <- subset(obj, mouse_or_rabbit == "rabbit")
  objlist[[i]] <- obj
  print(unique(obj$mouse_or_rabbit))
}

save(objlist, file = "rabbit_only_objlist_from_20230421_analysis_before_normalization.Rdata")

# normalize, dim reduction, clustering ####
dimreduclust <- function(obj){
  obj <- ScaleData(obj)
  obj <- RunPCA(obj,
                verbose = T,
                seed.use = 2023)
  obj <- RunUMAP(obj,
                 reduction = "pca",
                 dims = 1:30,
                 seed.use = 2023)
  obj <- FindNeighbors(obj,
                       reduction = "pca",
                       dims = 1:30)
  obj <- FindClusters(obj,
                      resolution = 0.5,
                      random.seed = 2023)
  return(obj)
  
}

processed_lst <- vector("list", length(objlist))
for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  print(paste0("now at ", i))
  # proceed with normalization
  obj <- PercentageFeatureSet(objlist[[i]], pattern = c("^MT-"), col.name = "percent.mt")
  obj <- FindVariableFeatures(obj, nfeatures = 3000)
  obj <- SCTransform(obj, vars.to.regress = "percent.mt", method = "glmGamPoi", 
                     return.only.var.genes = F, variable.features.n = 8000, 
                     seed.use = 2022, verbose = T)
  processed_lst[[i]] <- dimreduclust(obj)
}
names(processed_lst) <- sample_libs
save(processed_lst, file = "rabbit_only_objlist_clustered.RData")

# plot ####
for (i in 1:length(processed_lst)) {
  png(filename = paste0("./rabbit_only_dimplots/rabbit_dimplots_", 
                        names(processed_lst)[i], ".png"), width = 600, height = 400)
  p <- DimPlot(processed_lst[[i]], label = T, group.by = "seurat_clusters") +
    NoLegend() +
    ggtitle(paste0("Clustering of rabbit cells - ", names(processed_lst)[i])) +
    theme(text = element_text(size = 10), axis.text = element_text(size = 10))
  print(p)
  dev.off()
  
  png(filename = paste0("./rabbit_only_featplots/rabbit_GAD1_GAD2_SLC17A6_featplot_", 
                        names(processed_lst)[i], ".png"), width = 600, height = 400)
  p <- FeaturePlot(processed_lst[[i]], features = c("ENSOCUG00000012644", 
                                                    "ENSOCUG00000017480",
                                                    "ENSOCUG00000016389")) 
  p$patches$plots[[1]] <- p$patches$plots[[1]] + ggtitle("GAD1")
  p$patches$plots[[2]] <- p$patches$plots[[2]] + ggtitle("GAD2")
  p$labels$title <- "SLC17A6"
  print(p)
  dev.off()
  
  genes <- c("VIM", "MAP2", "PAX6", "GFAP")
  ensids <- c("ENSOCUG00000009222", "ENSOCUG00000005163", "ENSOCUG00000011075",
              "ENSOCUG00000006895")
  png(filename = paste0("./rabbit_only_featplots/rabbit_VIM_MAP2_PAX6_GFAP_featplot_", 
                        names(processed_lst)[i], ".png"), width = 600, height = 400)
  p <- FeaturePlot(processed_lst[[i]], features = ensids)
  p$patches$plots[[1]] <- p$patches$plots[[1]] + ggtitle("VIM")
  p$patches$plots[[2]] <- p$patches$plots[[2]] + ggtitle("MAP2")
  p$patches$plots[[3]] <- p$patches$plots[[3]] + ggtitle("PAX6")
  p$labels$title <- "GFAP"
  print(p)
  dev.off() 
}
 
