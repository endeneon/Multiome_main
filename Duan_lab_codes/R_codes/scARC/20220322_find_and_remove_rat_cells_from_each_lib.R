# Chuxuan Li 03/22/2022
# Extract rat cells from 024 data by clustering individual libraries

# init ####
{
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
}
# set threads and parallelization
plan("multisession", workers = 4)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

# load data ####
setwd("/data/FASTQ/Duan_Project_024/hybrid_output/")
# note this h5 file contains both atac-seq and gex information
# Read files, separate the ATACseq data from the .h5 matrix
h5list <- list.files(path = ".",
                     pattern = "filtered_feature_bc_matrix.h5",
                     recursive = T,
                     include.dirs = F, 
                     full.names = T)
objlist <- vector(mode = "list", length = length(h5list))
transformed_lst <- vector(mode = "list", length(h5list))

for (i in 1:length(objlist)){
  
  h5file <- Read10X_h5(filename = h5list[i])
  print(str_extract(string = h5list[i],
                    pattern = "[0-9]+-[0-6]"))
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(string = h5list[i],
                                                  pattern = "[0-9]+-[0-6]"))
  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list[i],
                               pattern = "-[0-6]"), "-"), "hr")
  print(unique(obj$time.ident))
  # check number of genes and cells
  print(paste0(i, " number of genes: ", nrow(obj),
               ", number of cells: ", ncol(obj)))
  objlist[[i]] <- obj
}
setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate")
save(objlist, file = "raw_data_list.RData")


# normalization ####
for (i in 1:length(objlist)){
  print(paste0("now at ", i))
  # proceed with normalization
  objlist[[i]] <-
    PercentageFeatureSet(objlist[[i]],
                         pattern = c("^MT-"),
                         col.name = "percent.mt")
  
  objlist[[i]] <- FindVariableFeatures(objlist[[i]], nfeatures = 3000)
  transformed_lst[[i]] <- SCTransform(objlist[[i]],
                                      vars.to.regress = "percent.mt",
                                      method = "glmGamPoi", 
                                      return.only.var.genes = F,
                                      variable.features.n = 8000,
                                      seed.use = 322,
                                      verbose = T)
}

# dimensional reduction and clustering ####
dimreduclust <- function(obj){
  obj <- ScaleData(obj)
  obj <- RunPCA(obj,
                verbose = T,
                seed.use = 11)
  obj <- RunUMAP(obj,
                 reduction = "pca",
                 dims = 1:30,
                 seed.use = 11)
  obj <- FindNeighbors(obj,
                       reduction = "pca",
                       dims = 1:30)
  obj <- FindClusters(obj,
                      resolution = 0.5,
                      random.seed = 11)
  return(obj)
  
}
for (i in 1:length(transformed_lst)){
  objlist[[i]] <- dimreduclust(transformed_lst[[i]])
}

save(transformed_lst, file = "transformed_list.RData")
rm(transformed_lst)

# check and remove rat cells ####
cleanobj_lst <- vector(mode = "list", length = length(objlist))

setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/QC_for_select_rat_plots/")
jpeg(paste0(str_extract(h5list[15], "[0-9]+-[0-6]"), "_dimplot_by_cluster.jpeg"))
p <- DimPlot(objlist[[15]],
        label = T,
        group.by = "seurat_clusters") +
  NoLegend() +
  ggtitle(str_extract(h5list[15], "[0-9]+-[0-6]")) +
  theme(text = element_text(size = 12))
print(p)
dev.off()

jpeg(paste0(str_extract(h5list[15], "[0-9]+-[0-6]"), "_Gfap.jpeg"))
p <- FeaturePlot(objlist[[15]], features = "Gfap") +
  ggtitle(str_extract(h5list[15], "[0-9]+-[0-6]"))
print(p)
dev.off()

jpeg(paste0(str_extract(h5list[15], "[0-9]+-[0-6]"), "_GAD1.jpeg"))
p <- FeaturePlot(objlist[[15]], features = "GAD1") +
  ggtitle(str_extract(h5list[15], "[0-9]+-[0-6]"))
print(p)
dev.off()

jpeg(paste0(str_extract(h5list[15], "[0-9]+-[0-6]"), "_SLC17A6.jpeg"))
p <- FeaturePlot(objlist[[15]], features = "SLC17A6") +
  ggtitle(str_extract(h5list[15], "[0-9]+-[0-6]"))
print(p)
dev.off()

jpeg(paste0(str_extract(h5list[15], "[0-9]+-[0-6]"), "_Gfap_vln.jpeg"))
p <- VlnPlot(objlist[[15]], "Gfap") +
  ggtitle(paste0(str_extract(h5list[15], "[0-9]+-[0-6]"), " Gfap expression by cluster"))
print(p)
dev.off()

FeaturePlot(objlist[[15]], features = "S100b")
# remove clusters that are rat cells, count
objlist[[5]]$rat.ident <- "human"
objlist[[5]]$rat.ident[objlist[[5]]$seurat_clusters %in% c(6)] <- "rat"

cleanobj_lst[[5]] <- subset(objlist[[5]], subset = rat.ident != "rat")
length(unique(cleanobj_lst[[5]]$seurat_clusters))

for (i in objlist){
  print(sum(i$rat.ident ==  "rat"))
}
sum(objlist[[5]]$rat.ident ==  "rat")
save(cleanobj_lst, file = "removed_rat_list.RData")


# map to demuxed barcodes ####
setwd("/data/FASTQ/Duan_Project_024/hybrid_output")
pathlist <- list.files(path = ".", full.names = T, pattern = ".best.tsv", recursive = T)
barcode_list <- vector(mode = "list", length = length(pathlist))
for (i in 1:length(pathlist)){
  barcode_list[[i]] <- read.delim(pathlist[i], header = F, row.names = NULL)
  colnames(barcode_list[[i]]) <- c("barcode", "line")
}

# sep lines in each group's barcode list, then assign line by barcode
for (i in 1:length(barcode_list)){
  lines <- unique(barcode_list[[i]]$line)
  print(lines)
  cleanobj_lst[[i]]$cell.line.ident <- "unmatched"
  for (j in 1:length(lines)){
    line_spec_barcodes <- barcode_list[[i]]$barcode[barcode_list[[i]]$line == lines[j]]
    cleanobj_lst[[i]]$cell.line.ident[cleanobj_lst[[i]]@assays$RNA@counts@Dimnames[[2]] %in% 
                                        line_spec_barcodes] <- lines[j]
  }
}

setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate")
save(cleanobj_lst, file = "demux_mapped_list.RData")



# check library 5
DefaultAssay(transformed_lst[[3]]) <- "SCT"
FeaturePlot(transformed_lst[[1]], features = "SLC17A7") + ggtitle("5-0")
FeaturePlot(transformed_lst[[2]], features = "SLC17A7") + ggtitle("5-1")
FeaturePlot(transformed_lst[[3]], features = "SLC17A7") + ggtitle("5-6")

# remake 5-0 plots
seurat_obj_2_plot <- dimreduclust(transformed_lst[[1]])
DimPlot(seurat_obj_2_plot,
        label = T,
        repel = T, 
        alpha = 0.8,
        raster = F) +
  NoLegend() +
  # theme(legend.location = "none") +
  ggtitle("5-0")

seurat_obj_2_plot$

FeaturePlot(seurat_obj_2_plot,
            features = "Gfap",
            alpha = 0.8)
FeaturePlot(seurat_obj_2_plot,
            features = "GAD1",
            alpha = 0.8)
FeaturePlot(seurat_obj_2_plot,
            features = "SLC17A6",
            alpha = 0.8)

# add barcodes, default will be "rat_astrocytes"

pathlist <- list.files(path = "/nvmefs/scARC_Duan_024_raw_fastq/GRCh38_output/group_05/common_barcodes_05", 
                       full.names = T, 
                       pattern = ".best.tsv", 
                       recursive = T)
barcode_list <- 
  read.delim(pathlist[1], 
             header = F, 
             row.names = NULL)
colnames(barcode_list) <- c("barcode", "line")
indiv_lines <-
  unique(barcode_list$line)

colnames(seurat_obj_2_plot)

seurat_obj_2_plot$indiv_line <-
  "Rat_astrocyte"
for (i in indiv_lines) {
  seurat_obj_2_plot$indiv_line[colnames(seurat_obj_2_plot) %in%
                                 barcode_list$barcode[barcode_list$line == i]] <-
    i
}

Idents(seurat_obj_2_plot) <-
  "indiv_line"
seurat_obj_2_plot$indiv_line <-
  factor(seurat_obj_2_plot$indiv_line,
         levels = c("CD_02",
                    "CD_05",
                    "CD_06",
                    "Rat_astrocyte"))
DimPlot(seurat_obj_2_plot,
        label = F,
        repel = F, 
        alpha = 0.5,
        cols = brewer.pal(n = 5,
                          name = "Set1"),
        raster = F)


Idents(seurat_obj_2_plot) <-
  "nFeature_RNA"
Idents(seurat_obj_2_plot)
Idents(seurat_obj_2_plot) <-
  "indiv_line"
Idents(seurat_obj_2_plot) <-
  "orig.ident"
DefaultAssay(seurat_obj_2_plot)
VlnPlot(seurat_obj_2_plot,
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3,
        pt.size = 0.1,
        alpha = 0.05) #+
  theme(axis.text.x = element_text(angle = 0))
unique(seurat_obj_2_plot$orig.ident)
