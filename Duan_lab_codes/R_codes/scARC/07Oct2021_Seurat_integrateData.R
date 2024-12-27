# Chuxuan Li 10/07/2021
# use the IntegrateData function in Seurat to integrate three datasets


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
library(viridis)
library(graphics)
library(patchwork)
library(cowplot)
library(gplots)

# set threads and parallelization
plan("multisession", workers = 6)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

# start with individual RNAseq objects that are separated by group and time
load("/home/cli/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/raw_group_2_group_8_seurat_objects.RData")

# merge and get the objects for each time point
time_0 <- merge(g_2_0, g_8_0)
time_1 <- merge(g_2_1, g_8_1)
time_6 <- merge(g_2_6, g_8_6)

# prepare the datasets
# process_data <- function(seurat_obj){
#   seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
#   SCTransform(seurat_obj,
#               vars.to.regress = "percent.mt",
#               method = "glmGamPoi",
#               variable.features.n = 8000,
#               verbose = T,
#               seed.use = 42)
#   return(seurat_obj)
# }

obj_lst <- c(time_0, time_1, time_6)

time_0_p <- PercentageFeatureSet(time_0, pattern = "^MT-", col.name = "percent.mt")
time_0_p <- SCTransform(time_0_p,
                        vars.to.regress = "percent.mt",
                        method = "glmGamPoi",
                        variable.features.n = 8000,
                        verbose = T,
                        seed.use = 42)

time_1_p <- PercentageFeatureSet(time_1, pattern = "^MT-", col.name = "percent.mt")
time_1_p <- SCTransform(time_1_p,
                        vars.to.regress = "percent.mt",
                        method = "glmGamPoi",
                        variable.features.n = 8000,
                        verbose = T,
                        seed.use = 42)

time_6_p <- PercentageFeatureSet(time_6, pattern = "^MT-", col.name = "percent.mt")
time_6_p <- SCTransform(time_6_p,
                        vars.to.regress = "percent.mt",
                        method = "glmGamPoi",
                        variable.features.n = 8000,
                        verbose = T,
                        seed.use = 42)

p_lst <- c(time_0_p, time_1_p, time_6_p)

# find anchors
features <- SelectIntegrationFeatures(object.list = p_lst)
p_lst_prepped <- PrepSCTIntegration(object.list = p_lst, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = p_lst_prepped, 
                                  normalization.method = "SCT",
                                  anchor.features = features)

# integrate
integrated <- IntegrateData(anchorset = anchors, 
                            normalization.method = "SCT")

# downstream analysis

DefaultAssay(integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
integrated <- RunPCA(integrated, 
                     assay = "integrated",
                     verbose = T, 
                     seed.use = 42)
integrated <- RunUMAP(integrated, 
                      reduction = "pca", 
                      dims = 1:30,
                      seed.use = 42)
integrated <- FindNeighbors(integrated, 
                            reduction = "pca", 
                            dims = 1:30)
integrated <- FindClusters(integrated, 
                           resolution = 0.5,
                           random.seed = 42)

save("integrated", file = "integrated.RData")
#load("../07Oct2021.RData")

DimPlot(integrated,
        repel = F, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of Seurat-combined data")

integrated_dimplot <- DimPlot(integrated,
                              cols = c("steelblue3", "mediumpurple4", "indianred3"),
                              # cols = c("cyan", "magenta", "yellow", "black"),
                              # cols = "Set2",
                              # cols = g_2_6_gex@meta.data$cell.line.ident,
                              group.by = "time.ident",
                              label = F) +
  ggtitle("integrated data by 0, 1, 6hr")

integrated_dimplot[[1]]$layers[[1]]$aes_params$alpha = .1
integrated_dimplot

DimPlot(integrated,
        # cols = c("cyan", "magenta", "yellow", "black"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        split.by = "time.ident",
        label = F) +
  ggtitle("integrated data separated by 0, 1, 6hr")




integrated_by_Seurat <- readRDS("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/integrated_by_Seurat.rds")


# plot stacked vlnplot
subtype_markers <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
                     "EBF1", "SERTAD4", # striatal
                     "FOXG1", # forebrain 
                     "TBR1", "FOXP2", "TLE4", # pallial glutamatergic
                     "FEZF2", "ADCYAP1", "TCERG1L", "SEMA3E", # subcerebral
                     "POU3F2", "CUX1", "BCL11B",  # cortical
                     "LHX2", # general cortex
                     #"EOMES", # hindbrain
                     "NPY", "SST", "DLX2", "DLX5", # inhibitory
                     "SATB2", "CBLN2", "CUX2", "NEFM", # excitatory
                     "VIM", "SOX2", "NES", #NPC
)
DefaultAssay(integrated) <- "RNA"

StackedVlnPlot(obj = integrated_by_Seurat, features = subtype_markers) +
  coord_flip()

# more stacked vlnplots for late/early response genes
new_markers <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
                 "EBF1", "SERTAD4", # striatal
                 "FOXG1", # forebrain 
                 "FOXP2", "TLE4", # pallial glutamatergic
                 "ADCYAP1", "TCERG1L", "SEMA3E", # subcerebral
                 "POU3F2", "CUX1", "BCL11B",  # cortical
                 "LHX2", # general cortex
                 "SST", # inhibitory
                 "SATB2", "CUX2", "NEFM", # excitatory
                 "VIM", "SOX2", "NES", #NPC
                 "MAP2", "DCX"
)
StackedVlnPlot(obj = integrated_by_Seurat, features = new_markers) +
  coord_flip()

early_late_markers <- c("EGR1", "FOS", "FOSB", "NPAS4", # early response
                        "VGF", "BDNF", "NRN1", "PNOC" # late response
)
StackedVlnPlot(obj = integrated_by_Seurat, features = early_late_markers) +
  coord_flip()


# try a dendrogram plot
integrated <- BuildClusterTree(object = integrated_by_Seurat, 
                               assay = "SCT",
                               features = subtype_markers,
                               verbose = T)

# integrated <- BuildClusterTree(object = integrated, 
#                                assay = "integrated",
#                                features = integrated@assays$integrated@data@Dimnames[[1]],
#                                verbose = T)

PlotClusterTree(integrated_by_Seurat, 
                direction = "rightwards",
                use.edge.length = F)

# saveRDS(integrated, file = "integrated_by_Seurat.rds")

# violin plot for striatal markers
VlnPlot(integrated, 
        features = c("EBF1", "SERTAD4"),
        split.by = "time.ident",
        cols = viridis(3),
        assay = "RNA",
        ncol = 1,
        group.by = "seurat_clusters",
        pt.size = 0) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1) +
  ylim(c(-1,2.5)) 

# vlnplot for forebrain marker
VlnPlot(integrated, 
        features = c("FOXG1"),
        split.by = "time.ident",
        cols = viridis(3),
        assay = "RNA",
        ncol = 1,
        group.by = "seurat_clusters",
        pt.size = 0) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1) +
  ylim(c(-1,2.5)) 

# violin plot for pallial glutamatergic markers
VlnPlot(integrated, 
        features = c("TBR1", "FOXP2", "TLE4"),
        split.by = "time.ident",
        cols = viridis(3),
        assay = "RNA",
        ncol = 1,
        group.by = "seurat_clusters",
        pt.size = 0) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1) +
  ylim(c(-1,2.5)) 

# violin plot for subcerebral markers
VlnPlot(integrated, 
        features = c("FEZF2", "ADCYAP1", "TCERG1L", "SEMA3E"),
        split.by = "time.ident",
        cols = viridis(3),
        assay = "RNA",
        ncol = 1,
        group.by = "seurat_clusters",
        pt.size = 0) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1) +
  ylim(c(-1,2.5)) 

# violin plot for cortical markers
VlnPlot(integrated, 
        features = c("POU3F2", "CUX1", "BCL11B"),
        split.by = "time.ident",
        cols = viridis(3),
        assay = "RNA",
        ncol = 1,
        group.by = "seurat_clusters",
        pt.size = 0) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1) +
  ylim(c(-1,2.5)) 

# violin plot for general cortex markers
VlnPlot(integrated, 
        features = c("LHX2", "EOMES", "PLXND1"),
        split.by = "time.ident",
        cols = viridis(3),
        assay = "RNA",
        ncol = 1,
        group.by = "seurat_clusters",
        pt.size = 0) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1) +
  ylim(c(-1,2.5)) 

# violin plot for inhibitory markers
VlnPlot(integrated, 
        features = c("NPY", "SST", "DLX2", "DLX5"),
        split.by = "time.ident",
        cols = viridis(3),
        assay = "RNA",
        ncol = 1,
        group.by = "seurat_clusters",
        pt.size = 0) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1) +
  ylim(c(-1,2.5)) 

# violin plot for excitatory markers
VlnPlot(integrated, 
        features = c("SATB2", "CBLN2", "CUX2", "NEFM"),
        split.by = "time.ident",
        cols = viridis(3),
        assay = "RNA",
        ncol = 1,
        group.by = "seurat_clusters",
        pt.size = 0) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.1) +
  ylim(c(-1,2.5)) 


# cell type identification
DefaultAssay(integrated) <- "RNA"
# create a list to store all markers
num_clust <- length(levels(integrated$seurat_clusters))
lst_of_markers <- vector(mode = "list", length = num_clust)

# loop and calculate all markers
for (i in 0:(num_clust-1)){
  ms <- FindConservedMarkers(integrated, ident.1 = i, grouping.var = "time.ident", verbose = FALSE)
  
  lst_of_markers[i] <- ms
}

# test FindConservedMarkers
# cluster_0_markers <- FindConservedMarkers(integrated, ident.1 = 0, grouping.var = "time.ident", verbose = FALSE)
# head(cluster_0_markers)

for (m in subtype_markers){
  if (m %in% rownames(cluster_0_markers)){
    print(paste(m, which(rownames(cluster_0_markers) == m)))
  }
}



# assign cluster identities based on the violin plots
new.cluster.ids <- 
  c("GAD1+/GAD2+ GABA", "NEFM+/CUX2- glut", "NEFM-/CUX2+ glut", "forebrain NPC", "GAD1+ GABA (less mature)", 
    "immature neuron", "NEFM+/CUX2-, SST+ glut", "NEFM-/CUX2+, ADCYAP1+ glut", "NEFM+/CUX2- glut", "SEMA3E+ glut", 
    "subcerebral immature neuron", "GAD1+ GABA", "MAP2+ NPC", "NEFM+/CUX2- glut", "NEFM-/CUX2+ glut", 
    "NEFM+/CUX2-, ADCYAP1+ glut", "NEFM+/CUX2-, SST+ glut", "NPC", "SEMA3E+ glut", "GAD1+/GAD2+ GABA", 
    "immature forebrain glut", "NEFM-/CUX2+, TCERG1L+ glut", "immature forebrain glut", "NEFM+/CUX2+, SST+ glut", "NEFM-/CUX2+, ADCYAP1+ glut", 
    "NEFM-/CUX2+, TCERG1L+ glut")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated$seurat_clusters))
names(new.cluster.ids) <- levels(integrated)

integrated_renamed_1 <- RenameIdents(integrated, new.cluster.ids)


labeled_dimplot <- DimPlot(integrated_renamed_1, 
                           reduction = "umap", 
                           label = TRUE, 
                           repel = T, 
                           pt.size = 0.3) +
  ggtitle("clusters labeled by cell type")
# set transparency
labeled_dimplot[[1]]$layers[[1]]$aes_params$alpha = .1
labeled_dimplot

save("integrated_renamed", file = "integrated_renamed.RData")


# get the count of each marker gene for each cluster
# initialize a matrix to store. Rows are genes, columns are clusters
marker_counts <- as.data.frame(matrix(data = NA, 
                                      nrow = length(subtype_markers), 
                                      ncol = length(levels(integrated$seurat_clusters))))
rownames(marker_counts) <- subtype_markers

hist_res <- vector(mode = "list", length = length(subtype_markers) * length(levels(integrated$seurat_clusters)))

for (i in 1:length(subtype_markers)){
  for (j in 1:length(levels(integrated$seurat_clusters))){
    print(paste("i: ", i, "j: ", j))
    marker <- subtype_markers[i]
    cluster_num <- levels(integrated$seurat_clusters[j])
    counts <- integrated@assays$RNA@counts[(integrated@assays$RNA@counts@Dimnames[[1]] %in% marker), 
                                           integrated$seurat_clusters %in% cluster_num]
    # print(hist(counts, 
    #      breaks = 10,
    #      plot = F))
    res <- hist(counts,
                breaks = 10,
                plot = F)$density
    hist_res[i*j] <- res
  }
}


# check cell line composition of each cluster
# first read in all the cell line barcode files
g_2_0_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_0_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
g_2_0_CD_27_barcodes <- unlist(g_2_0_CD_27_barcodes)
g_2_0_CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_0_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
g_2_0_CD_54_barcodes <- unlist(g_2_0_CD_54_barcodes)


g_2_1_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_1_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
g_2_1_CD_27_barcodes <- unlist(g_2_1_CD_27_barcodes)
g_2_1_CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_1_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
g_2_1_CD_54_barcodes <- unlist(g_2_1_CD_54_barcodes)

g_2_6_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_6_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
g_2_6_CD_27_barcodes <- unlist(g_2_6_CD_27_barcodes)
g_2_6_CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_6_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
g_2_6_CD_54_barcodes <- unlist(g_2_6_CD_54_barcodes)


# import group 8 cell line barcodes
g_8_0_CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_0_common_CD08_CD25_CD26.best.CD_08", 
           col_names = FALSE)
g_8_0_CD_08_barcodes <- unlist(g_8_0_CD_08_barcodes)

g_8_0_CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_0_common_CD08_CD25_CD26.best.CD_25", 
           col_names = FALSE)
g_8_0_CD_25_barcodes <- unlist(g_8_0_CD_25_barcodes)

g_8_0_CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_0_common_CD08_CD25_CD26.best.CD_26", 
           col_names = FALSE)
g_8_0_CD_26_barcodes <- unlist(g_8_0_CD_26_barcodes)


g_8_1_CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_1_common_CD08_CD25_CD26.best.CD_08", 
           col_names = FALSE)
g_8_1_CD_08_barcodes <- unlist(g_8_1_CD_08_barcodes)

g_8_1_CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_1_common_CD08_CD25_CD26.best.CD_25", 
           col_names = FALSE)
g_8_1_CD_25_barcodes <- unlist(g_8_1_CD_25_barcodes)

g_8_1_CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_1_common_CD08_CD25_CD26.best.CD_26", 
           col_names = FALSE)
g_8_1_CD_26_barcodes <- unlist(g_8_1_CD_26_barcodes)


g_8_6_CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_6_common_CD08_CD25_CD26.best.CD_08", 
           col_names = FALSE)
g_8_6_CD_08_barcodes <- unlist(g_8_6_CD_08_barcodes)

g_8_6_CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_6_common_CD08_CD25_CD26.best.CD_25", 
           col_names = FALSE)
g_8_6_CD_25_barcodes <- unlist(g_8_6_CD_25_barcodes)

g_8_6_CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_6_common_CD08_CD25_CD26.best.CD_26", 
           col_names = FALSE)
g_8_6_CD_26_barcodes <- unlist(g_8_6_CD_26_barcodes)

# assign cell line identity
integrated$cell.line.ident <- NA
integrated$cell.line.ident[str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                             c(g_2_0_CD_27_barcodes,
                               g_2_1_CD_27_barcodes,
                               g_2_6_CD_27_barcodes)] <- "CD_27"
integrated$cell.line.ident[str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                             c(g_2_0_CD_54_barcodes,
                               g_2_1_CD_54_barcodes,
                               g_2_6_CD_54_barcodes)] <- "CD_54"
integrated$cell.line.ident[str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                             c(g_8_0_CD_08_barcodes,
                               g_8_1_CD_08_barcodes,
                               g_8_6_CD_08_barcodes)] <- "CD_08"
integrated$cell.line.ident[str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                             c(g_8_0_CD_25_barcodes,
                               g_8_1_CD_25_barcodes,
                               g_8_6_CD_25_barcodes)] <- "CD_25"
integrated$cell.line.ident[str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                             c(g_8_0_CD_26_barcodes,
                               g_8_1_CD_26_barcodes,
                               g_8_6_CD_26_barcodes)] <- "CD_26"
# check if anything remains unassigned
sum(is.na(integrated$cell.line.ident))
# get the number of clusters, time points, lines
clusters <- levels(integrated$seurat_clusters)
times <- unique(integrated$time.ident)
lines <- unique(integrated$cell.line.ident)
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
      s <- sum((integrated$time.ident == t) & (integrated$cell.line.ident %in% l) & (integrated$seurat_clusters == c))
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
  coord_flip()

DefaultAssay(integrated) <- "RNA"
# differential expression 
DotPlot(integrated, 
        features = early_late_markers, 
        cols = c("steelblue3", "mediumpurple3", "indianred3"), 
        dot.scale = 5, 
        split.by = "time.ident") +
  theme_bw()+
  # RotatedAxis() +
  coord_flip() 


DefaultAssay(integrated_renamed_1) <- "SCT"
DPnew2(integrated_renamed_1, 
       features = early_late_markers, 
       cols = c("royalblue3", "mediumpurple3", "red3"), 
       dot.scale = 6, 
       split.by = "time.ident") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 60, 
                                   hjust = 1), validate = TRUE) +
  coord_flip() 


new.cluster.ids <- 
  c("GAD1+/GAD2+ GABA", "NEFM+/CUX2- glut", "NEFM-/CUX2+ glut", "forebrain NPC", "GAD1+ GABA (less mature)", 
    "immature neuron", "NEFM+/CUX2-, SST+ glut", "NEFM-/CUX2+, ADCYAP1+ glut", "NEFM+/CUX2- glut", "SEMA3E+ glut", 
    "subcerebral immature neuron", "GAD1+ GABA", "MAP2+ NPC", "NEFM+/CUX2- glut", "NEFM-/CUX2+ glut", 
    "NEFM+/CUX2-, ADCYAP1+ glut", "NEFM+/CUX2-, SST+ glut", "NPC", "SEMA3E+ glut", "GAD1+/GAD2+ GABA", 
    "immature forebrain glut", "NEFM-/CUX2+, TCERG1L+ glut", "immature forebrain glut", "NEFM+/CUX2+, SST+ glut", "NEFM-/CUX2+, ADCYAP1+ glut", 
    "NEFM-/CUX2+, TCERG1L+ glut")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated$seurat_clusters))
names(new.cluster.ids) <- levels(integrated)

integrated_renamed_1 <- RenameIdents(integrated, new.cluster.ids)

save("integrated_renamed_1", file = "integrated_labeled_updated.RData")

# identify differentially expressed genes across conditions with scatterplot
theme_set(theme_cowplot())
Idents(integrated_renamed_1)
integrated_renamed_1$time.ident.mod <- NA
integrated_renamed_1$time.ident.mod[integrated_renamed_1$time.ident == "0hr"] <- "zero_hour"
integrated_renamed_1$time.ident.mod[integrated_renamed_1$time.ident == "1hr"] <- "one_hour"
integrated_renamed_1$time.ident.mod[integrated_renamed_1$time.ident == "6hr"] <- "six_hour"


# write a loop to generate plots repeatedly
for (id in new.cluster.ids){
  cells <- subset(integrated_renamed_1, idents = id) # choose the cell type
  cells <- NormalizeData(object = cells, assay = "RNA")
  Idents(cells) <- "time.ident.mod" # set default ident to be time
  avg.cells <- as.data.frame(log1p(AverageExpression(cells, verbose = FALSE)$RNA)) # calculate
  avg.cells$gene <- rownames(avg.cells) 
  
  # save plot as pdf: prepare file name
  clean_id <- gsub(pattern = "/", replacement = " ", x = id)
  print(clean_id)
  file.name <- paste0("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/g_2_g_8_combined_plots/Using IntegrateData/diff_gex_plots/", 
                      clean_id, 
                      "_diff_gex_scatterplot.pdf")
  pdf(file = file.name,
      width = 10.5, # The width of the plot in inches
      height = 5.27) # The height of the plot in inches
  
  # plot title is cell type + xvy hr
  plot1.title <- paste(id, "cells\n 0v1 hr")
  plot2.title <- paste(id, "cells\n 0v6 hr")
  p1 <- ggplot(avg.cells, 
               aes(x = zero_hour, y = one_hour)) + 
    geom_point(size = 0.5, color = "steelblue3", alpha = 0.6) + 
    geom_abline(slope = 1, intercept = 0, color = "red3") +
    xlim(0, 6) +
    ylim(0, 6) +
    ggtitle(plot1.title)
  p1 <- LabelPoints(plot = p1, 
                    points = early_late_markers, 
                    box.padding = 0,
                    min.segment.length = 0,
                    max.overlaps = 10000,
                    repel = TRUE)
  p2 <- ggplot(avg.cells, 
               aes(x = zero_hour, y = six_hour)) + 
    geom_point(size = 0.5, color = "steelblue3", alpha = 0.6) + 
    geom_abline(slope = 1, intercept = 0, color = "red3") +
    xlim(0, 6) +
    ylim(0, 6) +
    ggtitle(plot2.title)
  p2 <- LabelPoints(plot = p2, 
                    points = early_late_markers, 
                    box.padding = 0,
                    min.segment.length = 0,
                    max.overlaps = 10000,
                    repel = TRUE)
  print(p1 + p2)
  
  dev.off()
}

FeaturePlot(object = integrated_renamed_1, 
            features = c("FOS", "NPAS4", "BDNF", "VGF"),
            split.by = "time.ident")

# count number of cells in each cell type
cell_types <- unique(new.cluster.ids)
count_lst <- vector(mode = "list", length = length(cell_types))
names(count_lst) <- cell_types
for (i in 1:length(count_lst)){
  t <- cell_types[i]
  sum <- sum(integrated_renamed_1@active.ident %in% t)
  print(paste(t, ": ", sum))
  count_lst[i] <- sum
}

# group the subtypes into broad categories
integrated_renamed_1$celltype <- Idents(integrated_renamed_1)
integrated_renamed_1$broad.cell.type <- NA
integrated_renamed_1$broad.cell.type[integrated_renamed_1$celltype %in% 
                                       c("GAD1+/GAD2+ GABA", "GAD1+ GABA (less mature)", "GAD1+ GABA")] <- "GABA"
integrated_renamed_1$broad.cell.type[integrated_renamed_1$celltype %in% 
                                       c("NEFM+/CUX2- glut", "NEFM-/CUX2+ glut",
                                         "NEFM+/CUX2-, SST+ glut", "NEFM-/CUX2+, ADCYAP1+ glut",
                                         "NEFM+/CUX2-, ADCYAP1+ glut", "NEFM-/CUX2+, TCERG1L+ glut",
                                         "NEFM+/CUX2+, SST+ glut")] <- "glut"
integrated_renamed_1$broad.cell.type[integrated_renamed_1$celltype %in% 
                                       c("NPC", "immature neuron", "MAP2+ NPC", "subcerebral immature neuron")] <- "NPC"
# divide the glut
integrated_renamed_1$fine.cell.type <- NA
integrated_renamed_1$fine.cell.type[integrated_renamed_1$celltype %in% 
                                       c("GAD1+/GAD2+ GABA", "GAD1+ GABA (less mature)", "GAD1+ GABA")] <- "GABA"
integrated_renamed_1$fine.cell.type[integrated_renamed_1$celltype %in% 
                                       c("NEFM+/CUX2- glut", "NEFM+/CUX2-, SST+ glut", 
                                         "NEFM+/CUX2-, ADCYAP1+ glut")] <- "NEFM+/CUX2- glut"
integrated_renamed_1$fine.cell.type[integrated_renamed_1$celltype %in% 
                                       c("NEFM-/CUX2+, glut",
                                         "NEFM-/CUX2+, ADCYAP1+ glut",
                                         "NEFM-/CUX2+, TCERG1L+ glut"
                                         )] <- "NEFM-/CUX2+ glut"
integrated_renamed_1$fine.cell.type[integrated_renamed_1$celltype %in% 
                                       c("NPC", "immature neuron", "MAP2+ NPC", "subcerebral immature neuron")] <- "NPC"



# find variable markers for each time point
GABA_subset <- subset(integrated_renamed_1, subset = broad.cell.type == "GABA")
glut_subset <- subset(integrated_renamed_1, subset = broad.cell.type == "glut")
NPC_subset <- subset(integrated_renamed_1, subset = broad.cell.type == "NPC")
type_lst <- c(GABA_subset, glut_subset, NPC_subset)

for (t in type_lst){
  t_name <- unique(t$broad.cell.type)
  print(t_name)
  
  # differentially expressed genes comparing one big cell type with all others
  Idents(t) <- "time.ident"
  print(DefaultAssay(t))
  zero_one <- FindMarkers(t, 
                          ident.1 = "1hr", 
                          logfc.threshold = 0, 
                          min.pct = 0.1,
                          ident.2 = "0hr", 
                          #group.by = "time.ident",
                          test.use = "LR",
                          verbose = TRUE)
  zero_six <- FindMarkers(t, 
                          ident.1 = "6hr", 
                          logfc.threshold = 0,
                          min.pct = 0.1,
                          ident.2 = "0hr", 
                          #group.by = "time.ident",
                          test.use = "LR",
                          verbose = TRUE)
  
  write.table(zero_one, 
            file = paste0(t_name, "_DE_markers_1v0hr.csv"),
            col.names = T,
            row.names = F,
            sep = ",",
            quote = F)
  print(nrow(zero_one))
  write.table(zero_six, 
              file = paste0(t_name, "_DE_markers_6v0hr.csv"),
              col.names = T,
              row.names = F,
              sep = ",",
              quote = F)
}

markers = FindMarkers(seu, 
                      group.by = "my_clustering", 
                      ident.1 = "my_cluster1", 
                      ident.2 = "my_cluster2", 
                      min.cells.group = 1, 
                      min.cells.feature = 1,
                      min.pct = 0,
                      logfc.threshold = 0,
                      only.pos = FALSE
)


# check differentially expressed genes between cell types
# time_name <- unique(t$time.ident)
# print(time_name)
# Idents(t) <- "fine.cell.type"
# between_glut_diff <- FindMarkers(t, 
#                                  ident.1 = "NEFM+/CUX2- glut", 
#                                  ident.2 = "NEFM-/CUX2+ glut", 
#                                  verbose = TRUE)
# 
# # differentially expressed genes comparing one big cell type with all others
# Idents(t) <- "broad.cell.type"
# glut_diff <- FindMarkers(t, 
#                          ident.1 = "glut", 
#                          ident.2 = c("GABA", "NPC"), 
#                          verbose = TRUE)
# GABA_diff <- FindMarkers(t, 
#                          ident.1 = "GABA", 
#                          ident.2 = c("glut", "NPC"), 
#                          verbose = TRUE)
# NPC_diff <- FindMarkers(t, 
#                         ident.1 = "NPC", 
#                         ident.2 = c("GABA", "glut"), 
#                         verbose = TRUE)
# 
# write.csv(glut_diff, 
#           paste0(time_name, "_DE_markers_glut.csv"),
#           row.names = T)
# write.csv(GABA_diff, 
#           paste0(time_name, "_DE_markers_GABA.csv"), 
#           row.names = T)
# write.csv(NPC_diff, 
#           paste0(time_name, "_DE_markers_NPC.csv"), 
#           row.names = T)
# write.csv(between_glut_diff, 
#           paste0(time_name, "_DE_markers_DE_markers_within_glut.csv"),
#           row.names = T)



# plot dotplot split by time but only ext/inh
DefaultAssay(integrated_renamed_1) <- "SCT"
grant_markers <- c("FOS", "BDNF", "NPAS4", "SP4", "TCF4")

two_celltypes_only <- subset(x = integrated_renamed_1,
                               subset = broad.cell.type %in% c("glut", "GABA"))
DPnew2(two_celltypes_only, 
       features = grant_markers, 
       cols = rev(brewer.pal(9, "Blues")),
       dot.scale = 10, 
       split.by = "time.ident",
       group.by = "broad.cell.type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                  hjust = 1), validate = TRUE) +
  coord_flip() 




# by cell line dotplots
library(stringr)
two_celltypes_only@meta.data$cell.line.ident <- NA
two_celltypes_only@meta.data$cell.line.ident[str_sub(two_celltypes_only@assays$SCT@counts@Dimnames[[2]], end = -5L) %in% 
                                         c(g_2_0_CD_27_barcodes, g_2_1_CD_27_barcodes, g_2_6_CD_27_barcodes)] <- "CD_27"
two_celltypes_only@meta.data$cell.line.ident[str_sub(two_celltypes_only@assays$SCT@counts@Dimnames[[2]], end = -5L) %in% 
                                         c(g_8_0_CD_25_barcodes, g_8_1_CD_25_barcodes, g_8_6_CD_25_barcodes)] <- "CD_25"
two_celltypes_only@meta.data$cell.line.ident[str_sub(two_celltypes_only@assays$SCT@counts@Dimnames[[2]], end = -5L) %in% 
                                         c(g_8_0_CD_26_barcodes, g_8_1_CD_26_barcodes, g_8_6_CD_26_barcodes)] <- "CD_26"
two_celltypes_only@meta.data$cell.line.ident[str_sub(two_celltypes_only@assays$SCT@counts@Dimnames[[2]], end = -5L) %in% 
                                         c(g_8_0_CD_08_barcodes, g_8_1_CD_08_barcodes, g_8_6_CD_08_barcodes)] <- "CD_08"
two_celltypes_only@meta.data$cell.line.ident[str_sub(two_celltypes_only@assays$SCT@counts@Dimnames[[2]], end = -5L) %in% 
                                         c(g_2_0_CD_54_barcodes, g_2_1_CD_54_barcodes, g_2_6_CD_54_barcodes)] <- "CD_54"

CD27 <- subset(x = two_celltypes_only, subset = cell.line.ident == "CD_27")
CD54 <- subset(x = two_celltypes_only, subset = cell.line.ident == "CD_54")
CD25 <- subset(x = two_celltypes_only, subset = cell.line.ident == "CD_25")
CD26 <- subset(x = two_celltypes_only, subset = cell.line.ident == "CD_26")
CD08 <- subset(x = two_celltypes_only, subset = cell.line.ident == "CD_08")

# sum(is.na(integrated_renamed_1$cell.line.ident))
# sum((integrated_renamed_1$cell.line.ident) %in% "CD_54")
DPnew2(CD27, 
       features = heatmap_markers, 
       cols = "BrBG",
       dot.scale = 10, 
       split.by = "time.ident",
       group.by = "broad.cell.type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1), validate = TRUE) +
  coord_flip() +
  ggtitle("CD_26")

DPnew2(CD27, 
       features = grant_markers, 
       cols = rev(brewer.pal(9, "Blues")),
       dot.scale = 10, 
       split.by = "time.ident",
       group.by = "broad.cell.type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1), validate = TRUE) +
  coord_flip() +
  ggtitle("CD_27")

DPnew2(CD25, 
       features = heatmap_markers, 
       cols = rev(brewer.pal(9, "Blues")),
       dot.scale = 10, 
       split.by = "time.ident",
       group.by = "broad.cell.type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1), validate = TRUE) +
  coord_flip() +
  ggtitle("CD_25")

DPnew2(CD26, 
       features = heatmap_markers, 
       cols = rev(brewer.pal(9, "Blues")),
       dot.scale = 10, 
       split.by = "time.ident",
       group.by = "broad.cell.type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1), validate = TRUE) +
  coord_flip() +
  ggtitle("CD_26")

DPnew2(CD54, 
       features = grant_markers, 
       cols = rev(brewer.pal(9, "Blues")),
       dot.scale = 10, 
       split.by = "time.ident",
       group.by = "broad.cell.type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1), validate = TRUE) +
  coord_flip() +
  ggtitle("CD_54")

DPnew2(CD08, 
       features = grant_markers, 
       cols = rev(brewer.pal(9, "Blues")),
       dot.scale = 10, 
       split.by = "time.ident",
       group.by = "broad.cell.type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1), validate = TRUE) +
  coord_flip() +
  ggtitle("CD_08")

heatmap_markers <- c("JUNB", "NR4A1", "FOS", "BTG2", "ATF3", "NPAS4", "NR4A3", "DUSP1", "EGR1",
                     "VGF", "PCSK1", "DUSP4", "ATP1B1", "SCG2", "SLC7A5", "CREM", "BDNF", "NPTX1")
exp_df <- GetAvgExpScaled(CD26, 
                 features = heatmap_markers, 
                 cols = rev(brewer.pal(9, "Blues")),
                 dot.scale = 10, 
                 split.by = "time.ident",
                 group.by = "broad.cell.type") 
df_toplot <- rbind(exp_df[3, ] - exp_df[1, ], exp_df[5, ] - exp_df[1, ],
                   exp_df[4, ] - exp_df[2, ], exp_df[6, ] - exp_df[2, ])

df_4_heatmap <- as.data.frame(t(df_toplot))
df_4_heatmap <- df_4_heatmap[, c(2, 1, 4, 3)]


heatmap.2(as.matrix(df_4_heatmap),
          Rowv = F,
          Colv = F,
          dendrogram = "none",
          scale = "none",
          col = colorRampPalette(colors = c("steelblue3", "steelblue1", "white", "red3"))(100),
          trace = "none",
          margins = c(10,10),
          cellnote = round(as.matrix(df_4_heatmap), digits = 2),
          notecol = "black")
dev.off()

cell_lines <- c(CD08, CD25, CD26, CD27, CD54)
for (cl in cell_lines){
  raw_df <- GetAvgExpScaled(cl, 
                            features = heatmap_markers, 
                            cols = rev(brewer.pal(9, "Blues")),
                            dot.scale = 10, 
                            split.by = "time.ident",
                            group.by = "broad.cell.type") 
  rearr_df <- rbind(raw_df[3, ] - raw_df[1, ], raw_df[5, ] - raw_df[1, ],
                     raw_df[4, ] - raw_df[2, ], raw_df[6, ] - raw_df[2, ])
  
  heat_df <- as.data.frame(t(rearr_df))
  heat_df <- heat_df[, c(2, 1, 4, 3)]
  
  file_name <- (unique(cl$cell.line.ident))
  print(file_name)
  pdf(file = paste0(file_name, "_AvgExpScaled_heatmap.pdf"))
  heatmap.2(as.matrix(heat_df),
            Rowv = F,
            Colv = F,
            dendrogram = "none",
            scale = "none",
            col = colorRampPalette(colors = c("steelblue3", "steelblue1", "white", "red3"))(100),
            trace = "none",
            margins = c(10,10),
            cellnote = round(as.matrix(heat_df), digits = 2),
            notecol = "black")
  dev.off()
}

# rearrange order of columns for CD26
raw_df <- GetAvgExpScaled(CD26, 
                          features = heatmap_markers, 
                          cols = rev(brewer.pal(9, "Blues")),
                          dot.scale = 10, 
                          split.by = "time.ident",
                          group.by = "broad.cell.type") 
rearr_df <- rbind(raw_df[3, ] - raw_df[1, ], raw_df[5, ] - raw_df[1, ],
                  raw_df[4, ] - raw_df[2, ], raw_df[6, ] - raw_df[2, ])

heat_df <- as.data.frame(t(rearr_df))
heat_df <- heat_df[, c(2, 1, 3, 4)]

file_name <- (unique(CD26$cell.line.ident))
print(file_name)
pdf(file = paste0(file_name, "_AvgExpScaled_heatmap.pdf"))
heatmap.2(as.matrix(heat_df),
          Rowv = F,
          Colv = F,
          dendrogram = "none",
          scale = "none",
          col = colorRampPalette(colors = c("steelblue3", "steelblue1", "white", "red3"))(100),
          trace = "none",
          margins = c(10,10),
          cellnote = round(as.matrix(heat_df), digits = 2),
          notecol = "black")
dev.off()


# assign cell line identities to the main object
# separate into cell lines
# import cell ident barcodes
times <- c("0", "1", "6")
lines_g2 <- c("CD_27", "CD_54")
lines_g8 <- c("CD_08", "CD_25", "CD_26")
g_2_0_bc_lst <- vector(mode = "list", length = 2)
g_2_1_bc_lst <- vector(mode = "list", length = 2)
g_2_6_bc_lst <- vector(mode = "list", length = 2)
g_8_0_bc_lst <- vector(mode = "list", length = 3)
g_8_1_bc_lst <- vector(mode = "list", length = 3)
g_8_6_bc_lst <- vector(mode = "list", length = 3)

for (t in times){
  for (i in 1:length(lines_g2)){
    l <- lines_g2[i]
    file_name <- paste0("/nvmefs/scARC_Duan_018/GRCh38_mapped_only/common_barcodes/g_2_",
                        t,
                        "_common_CD27_CD54.best.",
                        l)
    print(file_name)
    bc <- read_csv(file_name, col_names = FALSE)
    bc <- unlist(bc)
    bc <- str_sub(bc, end = -3L)
    if (t %in% "0"){
      g_2_0_bc_lst[[i]] <- bc
    } else if (t %in% "1"){
      g_2_1_bc_lst[[i]] <- bc
    } else {
      g_2_6_bc_lst[[i]] <- bc
    }
  }
}

for (t in times){
  for (i in 1:length(lines_g8)){
    l <- lines_g8[i]
    file_name <- paste0("/nvmefs/scARC_Duan_018/GRCh38_mapped_only/common_barcodes/g_8_",
                        t,
                        "_common_CD08_CD25_CD26.best.",
                        l)
    print(file_name)
    bc <- read_csv(file_name, col_names = FALSE)
    bc <- unlist(bc)
    bc <- str_sub(bc, end = -3L)
    if (t %in% "0"){
      g_8_0_bc_lst[[i]] <- bc
    } else if (t %in% "1"){
      g_8_1_bc_lst[[i]] <- bc
    } else {
      g_8_6_bc_lst[[i]] <- bc
    }
  }
}

# assign cell line identities to the object
integrated_renamed_1$cell.line.ident <- NA
integrated_renamed_1$cell.line.ident[str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]], end = -7L) %in% 
                                         c(g_2_0_bc_lst[[1]], g_2_1_bc_lst[[1]], g_2_6_bc_lst[[1]])] <- "CD_27"
integrated_renamed_1$cell.line.ident[str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]], end = -7L) %in% 
                                         c(g_2_0_bc_lst[[2]], g_2_1_bc_lst[[2]], g_2_6_bc_lst[[2]])] <- "CD_54"

integrated_renamed_1$cell.line.ident[str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]], end = -7L) %in% 
                                         c(g_8_0_bc_lst[[1]], g_8_1_bc_lst[[1]], g_8_6_bc_lst[[1]])] <- "CD_08"
integrated_renamed_1$cell.line.ident[str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]], end = -7L) %in% 
                                         c(g_8_0_bc_lst[[2]], g_8_1_bc_lst[[2]], g_8_6_bc_lst[[2]])] <- "CD_25"
integrated_renamed_1$cell.line.ident[str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]], end = -7L) %in% 
                                         c(g_8_0_bc_lst[[3]], g_8_1_bc_lst[[3]], g_8_6_bc_lst[[3]])] <- "CD_26"
sum(is.na(integrated_renamed_1$cell.line.ident))

integrated_renamed_1$fine.cell.type[integrated_renamed_1$cell.type == "SEMA3E+ glut"] <- "SEMA3E+ glut"
integrated_renamed_1$fine.cell.type[integrated_renamed_1$cell.type %in% c("MAP2+ NPC", "NPC")] <- "NPC"
integrated_renamed_1$fine.cell.type[integrated_renamed_1$cell.type %in% c("immature forebrain glut", 
                                                                          "forebrain NPC")] <- "forebrain"
integrated_renamed_1$fine.cell.type[integrated_renamed_1$cell.type == "subcerebral immature neuron"] <- "subcerebral immature neuron"
integrated_renamed_1$fine.cell.type[integrated_renamed_1$cell.type == "immature neuron"] <- "immature neuron"

save("integrated_renamed_1", file = "integrated_with_all_labels.RData")

# separate cell lines
sep_by_line <- vector(mode = "list", length = 5)
lines = c("CD_27", "CD_54", "CD_08", "CD_25", "CD_26")
for (i in 1:length(lines)){
  l <- lines[i]
  sep_by_line[[i]] <- subset(integrated_renamed_1, subset = cell.line.ident %in% l)
}


sep_by_time <- vector(mode = "list", length = 5*3)
time_lst <- c("0hr", "1hr", "6hr")
count = 0
for (i in 1:length(time_lst)){
  t <- time_lst[i]
  print(t)
  for (j in 1:length(sep_by_line)) {
    count <- count + 1
    print(unique(sep_by_line[[j]]$cell.line.ident))
    sep_by_time[[count]] <- subset(sep_by_line[[j]], subset = time.ident == t)
  }
}

cell_types <- c("SEMA3E+ glut", "NPC", "forebrain", "immature neuron", 
                "subcerebral immature neuron", 
                "NEFM+/CUX2- glut", "NEFM-/CUX2+ glut", "GABA")
sep_by_cell_type <- vector(mode = "list", length = 5*3*length(cell_types))
order_lst <- vector(length = 5*3*length(cell_types))

count <- 0
for (i in 1:length(cell_types)){
  ct <- cell_types[i]
  print(ct)
  for (j in 1:length(sep_by_time)) {
    print(unique(sep_by_time[[j]]$fine.cell.type))
    count <- count + 1
    sep_by_cell_type[[count]] <- subset(sep_by_time[[j]], subset = fine.cell.type == ct)
    order_lst[count] <- paste(unique(sep_by_cell_type[[count]]$cell.line.ident)[1], 
                              unique(sep_by_cell_type[[count]]$time.ident[1]),
                              ct,
                              sep = "_")
  }
}

# strip all the unallowed symbols in the name of output files
orders <- gsub(pattern = ' ',
               replacement = "_",
               x = order_lst)

orders <- gsub(pattern = '\\+',
               replacement = "p",
               x = orders)

orders <- gsub(pattern = '-',
               replacement = "m",
               x = orders)

orders <- gsub(pattern = "\\/",
               replacement = "_",
               x = orders)
unique(orders)

# unique(order_lst)

for (i in 1:length(sep_by_cell_type)) {
  file_name <- paste0(orders[i], ".txt")
  write.table(x = gsub(pattern = "^([A-Z]+\\-1)_[0-9]_[0-9]", 
                       replacement = "\\1",
                       x = sep_by_cell_type[[i]]@assays$RNA@counts@Dimnames[[2]]),
              file = paste0("subsetted_barcodes_by_subtype_time_line/",
                            file_name), 
              quote = F,
              row.names = F, 
              col.names = F,
              sep = "\t")
}


test <- gsub(pattern = "^([A-Z]+\\-1)_[0-9]_[0-9]", 
             replacement = "\\1", 
             x = sep_by_cell_type[[1]]@assays$RNA@counts@Dimnames[[2]])

length(strsplit(test, ""))
nchar(test)
