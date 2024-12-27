# Chuxuan Li 05/03/2023
# separate human from mouse in 025-17 & 46 libraries data (rat + human)

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

# set parallelization plan to multisession w 1 worker 
#creates multiple R processes that work in parallel
plan("multisession", workers = 1)
# plan("sequential")
# plan()

#Set max number of expressions to evaluate in a single future. 
#defines the size of each processing 'chunk' 
# high values improve performance but consume memory 
options(expressions = 20000)

options(future.globals.maxSize = 207374182400)

# This code is for loading data from a specific directory containing 
#multiple '.h5' files with ATAC-seq and gene expression (gex) information.
# setwd("/data/FASTQ/Duan_Project_025_17_46_hg38_rn6/")
setwd("/data/FASTQ/Duan_Project_025_17_46/hg38_rn6_hybrid/")
# note this h5 file contains both atac-seq and gex information
# Read files, separate the ATACseq data from the .h5 matrix
h5list <- list.files(path = ".",
                     pattern = "filtered_feature_bc_matrix.h5",
                     recursive = T, #searing for files in sub directories 
                     include.dirs = F, #only include files not directories 
                     full.names = T) # this should only extract 3 h5 files from 40201 
h5list
#create an empty list to store dat objects 
objlist <- vector(mode = "list", length = length(h5list))
#create empty vector to store the names of the data objects 
names <- c()

# This code iterates through each '.h5' file in the 'objlist', reads 
#the data from the files, and creates Seurat objects for further analysis.
for (i in 1:length(objlist)){
  #read data frim the i-th .h5 file 
  h5file <- Read10X_h5(filename = h5list[i])
  # extract and print the info from the file name 
  print(str_extract(string = h5list[i],
                    pattern = "[0-9]+-[0-6]"))
  # Extract unique identifier from the file name and add it 
  #to the 'names' vector
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

setwd("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part1_hybrid_genome/")
save(objlist, file = "raw_data_list_mapped_to_hybrid_genome_17_46_new_libs.RData")


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

save(transformed_lst, file = "transformed_list_hybrid_genome_17_46_new_libs.RData")
save(scaled_lst, file = "scaled_clustered_obj_lst_hybrid_genome_17_46_new_libs.RData")
rm(transformed_lst)

# check and remove rat cells ####
cleanobj_lst <- vector(mode = "list", length = length(scaled_lst))

for (i in 1:length(scaled_lst)) {
  png(paste0("./hybrid_feature_plots_17_46/dimplot_by_cluster/", 
             str_extract(h5list[i], "[0-9]+-[0-6]"), "_dimplot_by_cluster.png"))
  p <- DimPlot(scaled_lst[[i]],
               label = T,
               group.by = "seurat_clusters") +
    NoLegend() +
    ggtitle(str_extract(h5list[i], "[0-9]+-[0-6]")) +
    theme(text = element_text(size = 12))
  print(p)
  dev.off()
  
  png(paste0("./hybrid_feature_plots_17_46/featplots/", 
             str_extract(h5list[i], "[0-9]+-[0-6]"), "_specific_genes.png"))
  p <- FeaturePlot(scaled_lst[[i]], 
                   features = c("Gfap", "S100b", 
                                "GAD1", "SLC17A6"), 
                   ncol = 2)
  print(p)
  dev.off()
  
  
  png(paste0("./hybrid_feature_plots_17_46/vlnplot_gfap/",
             str_extract(h5list[i], "[0-9]+-[0-6]"), "_Gfap_vln.png"))
  p <- VlnPlot(scaled_lst[[i]], "Gfap") +
    ggtitle(paste0(str_extract(h5list[i], 
                               "[0-9]+-[0-6]"), 
                   " Gfap expression by cluster"))
  print(p)
  dev.off()
}

# remove clusters that are rat cells, count
rat_clusters <- list(c(5), c(7), c(8),
                     c(8, 9), c(8), c(9))
for (i in 1:length(scaled_lst)) {
  c <- rat_clusters[[i]]
  cat("rat cluster number: ", c)
  scaled_lst[[i]]$rat.ident <- "human"
  scaled_lst[[i]]$rat.ident[scaled_lst[[i]]$seurat_clusters %in% c] <- "rat"
  human_barcodes <- colnames(scaled_lst[[i]])[scaled_lst[[i]]$rat.ident == "human"]
  write.table(human_barcodes, file = paste0("./025_17_46_human_only_barcodes/", 
                                            str_extract(h5list[i], "[0-9]+-[0-6]"),
                                            "_human_barcodes.txt"), 
              sep = "\t", quote = F, row.names = F, col.names = F)
  cleanobj_lst[[i]] <- subset(scaled_lst[[i]], subset = rat.ident != "rat")
}

rat_cell_count <- matrix(nrow = length(scaled_lst), ncol = 3, 
                           dimnames = list(names, 
                                           c("number of rat astrocytes", "total number of cells",
                                             "percentage of cells that are astrocytes (%)")))
for (i in 1:length(scaled_lst)){
  rat_cell_count[i, 1] <- (sum(scaled_lst[[i]]$rat.ident ==  "rat"))
  rat_cell_count[i, 2] <- (length(scaled_lst[[i]]$rat.ident))
  rat_cell_count[i, 3] <- rat_cell_count[i, 1] / rat_cell_count[i, 2] * 100
}
write.table(rat_cell_count, file = "rat_cell_count_by_lib_025_17_46.csv", sep = ",", quote = F)

save(cleanobj_lst, file = "removed_rat_list_025_17_46_mapped_to_hybrid.RData")

