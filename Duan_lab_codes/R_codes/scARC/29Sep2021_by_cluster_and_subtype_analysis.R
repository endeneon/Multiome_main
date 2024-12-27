# Chuxuan 9/29/2021
# Look at individual clusters and check for subtypes; 0, 1, 6 and groups not combined


rm(list=ls())
# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)
library(future)
library(RColorBrewer)
library(dplyr)

# set threads and parallelization
plan("multisession", workers = 3)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)


# load human only data, create Seurat objects
# load data use read10x_h5
# note this h5 file contains both atac-seq and gex information
g_2_0_read <- 
  Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_2_0/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_2_0_read_gex <- 
  CreateSeuratObject(counts = g_2_0_read$`Gene Expression`,
                     project = "g_2_0_read_gex")

g_2_1_read <- 
  Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_2_1/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_2_1_read_gex <- 
  CreateSeuratObject(counts = g_2_1_read$`Gene Expression`,
                     project = "g_2_1_read_gex")

g_2_6_read <- 
  Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_2_6/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_2_6_read_gex <- 
  CreateSeuratObject(counts = g_2_6_read$`Gene Expression`,
                     project = "g_2_6_read_gex")


g_8_0_read <- 
  Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_8_0/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_8_0_read_gex <- 
  CreateSeuratObject(counts = g_8_0_read$`Gene Expression`,
                     project = "g_8_0_read_gex")

g_8_1_read <- 
  Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_8_1/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_8_1_read_gex <- 
  CreateSeuratObject(counts = g_8_1_read$`Gene Expression`,
                     project = "g_8_1_read_gex")

g_8_6_read <- 
  Read10X_h5(filename = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/libraries_8_6/outs/filtered_feature_bc_matrix.h5")
# load gex information first
g_8_6_read_gex <- 
  CreateSeuratObject(counts = g_8_6_read$`Gene Expression`,
                     project = "g_8_6_read_gex")

# make a separate seurat object for subsequent operation
# Note: all genes with total counts = 0 have been pre-removed
g_2_0 <- g_2_0_read_gex
g_2_1 <- g_2_1_read_gex
g_2_6 <- g_2_6_read_gex
g_8_0 <- g_8_0_read_gex
g_8_1 <- g_8_1_read_gex
g_8_6 <- g_8_6_read_gex

# read barcode files
g_2_0_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_2_0_exclude_barcodes.txt",
                             header = F)
g_2_1_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_2_1_exclude_barcodes.txt",
                             header = F)
g_2_6_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_2_6_exclude_barcodes.txt",
                             header = F)
g_8_0_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_8_0_exclude_barcodes.txt",
                             header = F)
g_8_1_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_8_1_exclude_barcodes.txt",
                             header = F)
g_8_6_barcodes <- read.delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/QC_barcodes/g_8_6_exclude_barcodes.txt",
                             header = F)

# test6 <- c(unlist(g_2_6_barcodes, use.names = F), unlist(g_8_6_barcodes, use.names = F))

# select the data with barcodes not in the barcode files
# the subset() function does not take counts@dimnames, has to give cells an ident first
g_2_0@meta.data$exclude.status <- "include"
g_2_0@meta.data$exclude.status[colnames(g_2_0@assays$RNA@counts) %in% g_2_0_barcodes] <- "exclude"
g_2_0 <- subset(x = g_2_0, subset = exclude.status %in% "include")

g_2_1@meta.data$exclude.status <- "include"
g_2_1@meta.data$exclude.status[colnames(g_2_1@assays$RNA@counts) %in% g_2_1_barcodes] <- "exclude"
g_2_1 <- subset(x = g_2_1, subset = exclude.status %in% "include")

g_2_6@meta.data$exclude.status <- "include"
g_2_6@meta.data$exclude.status[colnames(g_2_6@assays$RNA@counts) %in% g_2_0_barcodes] <- "exclude"
g_2_6 <- subset(x = g_2_6, subset = exclude.status %in% "include")

g_8_0@meta.data$exclude.status <- "include"
g_8_0@meta.data$exclude.status[colnames(g_8_0@assays$RNA@counts) %in% g_8_0_barcodes] <- "exclude"
g_8_0 <- subset(x = g_8_0, subset = exclude.status %in% "include")

g_8_1@meta.data$exclude.status <- "include"
g_8_1@meta.data$exclude.status[colnames(g_8_1@assays$RNA@counts) %in% g_8_1_barcodes] <- "exclude"
g_8_1 <- subset(x = g_8_1, subset = exclude.status %in% "include")

g_8_6@meta.data$exclude.status <- "include"
g_8_6@meta.data$exclude.status[colnames(g_8_6@assays$RNA@counts) %in% g_8_0_barcodes] <- "exclude"
g_8_6 <- subset(x = g_8_6, subset = exclude.status %in% "include")

#save.image("pre-SCT_QCed_g2_g8.RData")


# assign time point ident
g_2_0@meta.data$time.ident <- "0hr"
g_2_1@meta.data$time.ident <- "1hr"
g_2_6@meta.data$time.ident <- "6hr"

g_8_0@meta.data$time.ident <- "0hr"
g_8_1@meta.data$time.ident <- "1hr"
g_8_6@meta.data$time.ident <- "6hr"

# assign group ident
g_2_0@meta.data$group.ident <- "group_2"
g_2_1@meta.data$group.ident <- "group_2"
g_2_6@meta.data$group.ident <- "group_2"

g_8_0@meta.data$group.ident <- "group_8"
g_8_1@meta.data$group.ident <- "group_8"
g_8_6@meta.data$group.ident <- "group_8"

save(list = c("g_2_0", "g_2_1", "g_2_6", "g_8_0", "g_8_1", "g_8_6"), file = "raw_group_2_group_8_seurat_objects.RData")


# create a list of objects
obj_lst <- c(g_2_0, g_2_1, g_2_6, g_8_0, g_8_1, g_8_6)

# use this function to do the steps before clustering
process_data_full <- function(seurat_obj){
  # store mitochondrial percentage in object meta data
  # human starts with "MT-", rat starts with "Mt-"
  seurat_obj <- PercentageFeatureSet(seurat_obj,
                                     pattern = c("^MT-"),
                                     col.name = "percent.mt",
                                     assay = 'RNA')
  #use SCTransform()
  seurat_obj <-
    SCTransform(seurat_obj,
                vars.to.regress = "percent.mt",
                method = "glmGamPoi",
                variable.features.n = 8000,
                verbose = T,
                seed.use = 42)
  seurat_obj <- RunPCA(seurat_obj,
                       seed.use = 42,
                       verbose = T)
  seurat_obj <- RunUMAP(seurat_obj,
                        dims = 1:30,
                        seed.use = 42,
                        verbose = T)
  seurat_obj <- FindNeighbors(seurat_obj,
                              dims = 1:30,
                              verbose = T)
  
  return(seurat_obj)
}

# initialize a list to store resulting objects from the first few steps 
processed_obj_lst <- vector(mode = "list", length = length(obj_lst))

# use a loop to apply the previous function to every object
for (i in 1:length(obj_lst)){
  new_obj <- process_data_full(obj_lst[[i]])
  processed_obj_lst[[i]] <- new_obj
}

# pick the first one and test the optimal resolution for this object, use as reference
g_2_0_processed <- processed_obj_lst[[1]]

g_2_0_processed <- FindClusters(g_2_0_processed,
                             verbose = T, 
                             random.seed = 42,
                             resolution = 0.14) # resolution of 0.14 returns 10 clusters

length(levels(g_2_0_processed@meta.data$seurat_clusters)) # check the number of clusters

# set the target cluster number 
target <- 10

# store the optimal resolution values
r_lst <- vector(mode = "list", length = length(obj_lst))
r_lst[[1]] <- 0.14

# find optimal resolution for each object to get the same number of clusters
lambda <- 0.95 # scaling factor for d
d <- 0.05 # amount changed each loop to adjust
maxsteps <- 15 # max number of rounds in the loop

# do the loop, find optimal resolution value for each object to reach 10 clusters
for (i in 2:length(processed_obj_lst)){
  obj <- processed_obj_lst[[i]]
  obj <- FindClusters(obj,
                      verbose = T, 
                      resolution = 0.14)
  n_cluster = length(levels(obj@meta.data$seurat_clusters))
  r = 0.14
  nstep = 0
  while ((n_cluster != target) && nstep <= maxsteps){
    if (n_cluster > target){
      r = r - d
      d = d * lambda
    } else {
      r = r + d
      d = d * lambda
    }
    obj <- FindClusters(obj,
                        verbose = T, 
                        resolution = r)
    n_cluster = length(levels(obj@meta.data$seurat_clusters))
    nstep = nstep + 1
  }
  r_lst[[i]] <- r
}

# initialize list to store the final clustering results
clustered_obj_lst <- vector(mode = "list", length = length(obj_lst))

# use the optimal clustering resolution to cluster the objects again, store the clustering information
for (i in 1:length(processed_obj_lst)){
  obj <- processed_obj_lst[[i]] # extract the object from processed list
  r <- r_lst[[i]] # get the resolution value
  # store the resulting object from clustering
  clustered_obj_lst[[i]] <- FindClusters(obj,
                                         verbose = T, 
                                         resolution = r,
                                         random.seed = 42)
}

# get the clustered objects
g_2_0_p <- clustered_obj_lst[[1]]
g_2_1_p <- clustered_obj_lst[[2]]
g_2_6_p <- clustered_obj_lst[[3]]
g_8_0_p <- clustered_obj_lst[[4]]
g_8_1_p <- clustered_obj_lst[[5]]
g_8_6_p <- clustered_obj_lst[[6]]



# plot DimPlots, check if all have 10 clusters. If not, manually adjust
DimPlot(g_2_0_p, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of group 2 at 0hr")

DimPlot(g_2_1_p, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of group 2 at 1hr")


# adjust resolution
g_2_6_p <- FindClusters(processed_obj_lst[[3]],
                        verbose = T, 
                        resolution = 0.205,
                        random.seed = 42)
DimPlot(g_2_6_p, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of group 2 at 6hr")


g_8_0_p <- FindClusters(processed_obj_lst[[4]],
                        verbose = T, 
                        resolution = 0.1,
                        random.seed = 42)

DimPlot(g_8_0_p, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of group 8 at 0hr")

g_8_1_p <- FindClusters(processed_obj_lst[[5]],
                        verbose = T, 
                        resolution = 0.11,
                        random.seed = 42)

DimPlot(g_8_1_p, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of group 8 at 1hr")

g_8_6_p <- FindClusters(processed_obj_lst[[6]],
                        verbose = T, 
                        resolution = 0.0707,
                        random.seed = 42)

DimPlot(g_8_6_p, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of group 8 at 6hr")


# plot feature plots with a set of marker genes
genes <- c("GAD1", "SLC17A6", "FOS", "BDNF", "VIM", "SOX2")
FeaturePlot(g_2_0_p, 
            features = genes, 
            pt.size = 0.2,
            ncol = 6)
FeaturePlot(g_2_1_p, 
            features = genes, 
            pt.size = 0.2,
            ncol = 6)
FeaturePlot(g_2_6_p, 
            features = genes, 
            pt.size = 0.2,
            ncol = 6)

FeaturePlot(g_8_0_p, 
            features = genes, 
            pt.size = 0.2,
            ncol = 6)
FeaturePlot(g_8_1_p, 
            features = genes, 
            pt.size = 0.2,
            ncol = 6)
FeaturePlot(g_8_6_p, 
            features = genes, 
            pt.size = 0.2,
            ncol = 6)

# check subtypes using the combined object
# load the combined object
combined <- readRDS("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/DESeq2/combined_g2_g8_USE_THIS.rds")

markers_for_vln <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
                     "EBF1", "SERTAD4", # striatal
                     "FOXG1", # forebrain 
                     "TBR1", # pallial glutamatergic
                     "POU3F2", # cortical
                     "NPY", "SST", "DLX2", "DLX5", # inhibitory
                     "SATB2", "CBLN2", "CUX2", "NEFM", # excitatory
                     "VIM", "SOX2", "NES" #NPC
                     )
StackedVlnPlot(obj = combined, features = markers_for_vln) +
  coord_flip()

# plot combined dimplot again with time point identities (He Xin's request)
trasparent_combined_dimplot <- DimPlot(combined,
                                       cols = c("steelblue3", "mediumpurple4", "indianred3"),
                                       # cols = c("cyan", "magenta", "yellow", "black"),
                                       # cols = "Set2",
                                       # cols = g_2_6_gex@meta.data$cell.line.ident,
                                       group.by = "time.ident",
                                       label = F) +
  ggtitle("clustering timed by 0, 1, 6hr")
trasparent_combined_dimplot[[1]]$layers[[1]]$aes_params$alpha = .1
trasparent_combined_dimplot

# plot 3 dimplots with only one time point colored
DimPlot(combined,
        cols = c("steelblue3", "transparent", "transparent"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("clustering of merged object, 0hr")

DimPlot(combined,
        cols = c("transparent", "mediumpurple4", "transparent"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("clustering of merged object, 1hr")

DimPlot(combined,
        cols = c("transparent", "transparent", "indianred3"),
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("clustering of merged object, 6hr")


# plot combined marker gene expressions (He Xin's request)
FeaturePlot(combined, 
            features = genes, 
            pt.size = 0.2,
            ncol = 3, 
            cols = c("lightgrey", "#0083cf"), 
            keep.scale = "all")

# other marker genes 
subtype_markers <- c("EBF1", "SERTAD4", # striatal
                     "FOXG1", # forebrain 
                     "TBR1", "FOXP2", "TLE4", # pallial glutamatergic
                     "FEZF2", "ADCYAP1", "TCERG1L", "SEMA3E", # subcerebral
                     "POU3F2", "CUX1", "BCL11B", "FPXP2", # cortical
                     "LHX2", "EOMES", "PLXND1", # general cortex
                     "NPY", "SST", "DLX2", "DLX5", # inhibitory
                     "SATB2", "CBLN2", "CUX2", "NEFM" # excitatory
                     )
cortical_markers <- c("POU3F2", "CUX1", "BCL11B", "FOXP2")

exc_markers <- c("BCL11B","TBR1", "SATB2", "POU3F2", "CUX1")
inh_markers <- c("SST", "CALB1", "CALB2", "PVALB")
FeaturePlot(g_2_0_p, 
            features = cortical_markers, 
            pt.size = 0.2,
            ncol = 4, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(g_2_0_p, 
            features = exc_markers, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

FeaturePlot(g_2_0_p, 
            features = inh_markers, 
            pt.size = 0.2,
            ncol = 5, 
            cols = c("lightgrey", "#0083cf"))

# look at combined object again
clusters <- levels(combined$seurat_clusters)
cluster_cell_count_by_time <- data.frame(group_2_0hr = 1:27,
                                         group_2_1hr = 1:27,
                                         group_2_6hr = 1:27,
                                         group_8_0hr = 1:27,
                                         group_8_1hr = 1:27,
                                         group_8_6hr = 1:27,
                                         stringsAsFactors = F)
for (i in 1:length(clusters)){
  c <- clusters[i]

  cluster_cell_count_by_time$group_2_0hr[i] <- sum((combined$time.ident == "0hr") & 
                                                  (combined$group.ident == "group_2") & 
                                                  (combined$seurat_clusters == c))
  cluster_cell_count_by_time$group_2_1hr[i] <- sum((combined$time.ident == "1hr") & 
                                                  (combined$group.ident == "group_2") & 
                                                  (combined$seurat_clusters == c))
  cluster_cell_count_by_time$group_2_6hr[i] <- sum((combined$time.ident == "6hr") & 
                                                  (combined$group.ident == "group_2") & 
                                                  (combined$seurat_clusters == c))
  cluster_cell_count_by_time$group_8_0hr[i] <- sum((combined$time.ident == "0hr") & 
                                                  (combined$group.ident == "group_8") & 
                                                  (combined$seurat_clusters == c))
  cluster_cell_count_by_time$group_8_1hr[i] <- sum((combined$time.ident == "1hr") & 
                                                  (combined$group.ident == "group_8") & 
                                                  (combined$seurat_clusters == c))
  cluster_cell_count_by_time$group_8_6hr[i] <- sum((combined$time.ident == "6hr") & 
                                                  (combined$group.ident == "group_8") & 
                                                  (combined$seurat_clusters == c))
}

cluster_cell_count_by_time_percentage <-
  as.matrix(cluster_cell_count_by_time)/rowSums(cluster_cell_count_by_time)
cluster_cell_count_by_time_percentage <- 
  as.data.frame(cluster_cell_count_by_time_percentage)
cluster_cell_count_by_time_percentage$total_cells <- rowSums(cluster_cell_count_by_time)

rownames(cluster_cell_count_by_time) <- 0:26


# # now do not analyze the time points separately; combined times within the two groups 
# group_2 <- merge(g_2_0, c(g_2_1, g_2_6))
# group_8 <- merge(g_8_0, c(g_8_1, g_8_6))

# # cluster within the two groups
# group_2_p <- process_data(group_2)
# group_8_p <- process_data(group_8)
# 
# group_2_p <- FindClusters(group_2_p,
#                         verbose = T, 
#                         resolution = 0.501,
#                         random.seed = 42)
# group_8_p <- FindClusters(group_8_p,
#                         verbose = T, 
#                         resolution = 0.5,
#                         random.seed = 42)

# # get the number of cells in each time, divide by total number of cell 
# clusters_g2 <- levels(group_2_p$seurat_clusters)
# clusters_g8 <- levels(group_8_p$seurat_clusters)
# 
# 
# for (c in clusters)
# sum(group_2$time.ident == "0hr")
# sum(group_2$time.ident == "0hr")

GABA <- subset(combined, subset = cell.type.ident == "GABA")
GABA <- RunPCA(GABA,
               seed.use = 42,
               verbose = T)
GABA <- RunUMAP(GABA,
                dims = 1:30,
                seed.use = 42,
                verbose = T)
GABA <- FindNeighbors(GABA,
                      dims = 1:30,
                      verbose = T)
GABA <- FindClusters(GABA,
                     verbose = T, 
                     resolution = 0.5,
                     random.seed = 42)

clusters <- levels(GABA$seurat_clusters)
times <- unique(GABA$time.ident)
lines <- unique(GABA$cell.line.ident)
lines <- lines[!(is.na(lines))]

nc <- length(clusters)
nt <- length(times)
nl <- length(lines)


df_4_plot <- data.frame(value = rep_len(0, length.out = (nc * nt * nl)),
                        cluster = rep_len(NA, length.out = (nc * nt * nl)),
                        time = rep_len(NA, length.out = (nc * nt * nl)),
                        cell_line = rep_len(NA, length.out = (nc * nt * nl)),
                        stringsAsFactors = F)

i <- 1
for (c in clusters){
  for (t in times){
    for (l in lines){
      # print(l)
      s <- sum((GABA$time.ident == t) & (GABA$cell.line.ident %in% l) & (GABA$seurat_clusters == c))
      df_4_plot$value[i] <- s
      df_4_plot$cluster[i] <- c
      df_4_plot$time[i] <- t
      df_4_plot$cell_line[i] <- l
      i = i + 1
    }
  }
}

df_4_plot$cluster <- factor(df_4_plot$cluster, 
                            levels = as.character(0:nc))
#c("0",  "1",  "2",  "3",  "4",  5  6  7  8  9 10 
#11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26))
df_4_plot$time <- factor(df_4_plot$time, 
                         levels = c("6hr", "1hr", "0hr"))
df_4_plot$cell_line <- factor(df_4_plot$cell_line)

hist(df_4_plot$value)
ggplot(data = df_4_plot, aes(x = cell_line, 
                             y = cluster,
                             fill = time,
                             size = value,
                             group = time)) +
  scale_fill_manual(values = brewer.pal(3, "Dark2")) +
  geom_point(shape = 21,
             na.rm = TRUE,
             position = position_dodge(width = 0.9)) +
  scale_size_area() +
  theme_bw() +
  coord_flip() +
  ggtitle("GABAergic cells only")

# get the top markers for time = 0hr, try to use these as anchors
time_0 <- merge(g_2_0, g_8_0)

time_0_p <- process_data(time_0)

time_0_p <- FindClusters(time_0_p,
                         verbose = T, 
                         resolution = 0.5,
                         random.seed = 42)
DimPlot(time_0_p, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of merged data at 0hr")


time_0_markers <- FindAllMarkers(time_0_p)
top_markers <- time_0_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.table(top_markers, sep = "," , file = "top_markers_0hr.csv")


# test the Seurat Find anchors function
test_anchors <- FindIntegrationAnchors(c(g_2_0, g_2_1, g_2_6),
                       anchor.features = 100)
# test integration using the Seurat integrate function
test_combined <- IntegrateData(anchorset = test_anchors)

test_combined_p <- process_data(test_combined)

DimPlot(test_combined_p, 
        label = F) +
  ggtitle("cluster by time for integrated data for group 2")

test_combined_p <- FindClusters(test_combined_p,
                                verbose = T, 
                                resolution = 0.5,
                                random.seed = 42)
DimPlot(test_combined_p, 
        label = T) +
  NoLegend() +
  ggtitle("clustering for integrated data for group 2")
