# Chuxuan Li 12/23/2021
# RNAseq data analysis by SCT-normalizing 6 libraries individually first, then
# use IntegrateData() in Seurat to integrate 6 libraries


# init ####
library(Seurat)
library(future)
library(readr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)
library(patchwork)
library(cowplot)
library(gplots)

# set threads and parallelization
plan("multisession", workers = 1)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 107374182400)


# load data ####
# note this h5 file contains both atac-seq and gex information
mat_lst <- vector(mode = "list", length = 6L)
mat_lst[[1]] <- 
  Read10X_h5(filename = "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_0/outs/filtered_feature_bc_matrix.h5")
mat_lst[[2]] <- 
  Read10X_h5(filename = "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_1/outs/filtered_feature_bc_matrix.h5")
mat_lst[[3]] <- 
  Read10X_h5(filename = "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_6/outs/filtered_feature_bc_matrix.h5")
mat_lst[[4]] <- 
  Read10X_h5(filename = "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_0/outs/filtered_feature_bc_matrix.h5")
mat_lst[[5]] <- 
  Read10X_h5(filename = "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_1/outs/filtered_feature_bc_matrix.h5")
mat_lst[[6]] <- 
  Read10X_h5(filename = "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_6/outs/filtered_feature_bc_matrix.h5")

# load gex information and create objects
obj_lst <- vector(mode = "list", length = 6L)
normed_lst <- vector(mode = "list", length = 6L)
for (i in 1:length(obj_lst)){
  obj <- CreateSeuratObject(counts = mat_lst[[i]]$`Gene Expression`,
                            project = "pilot")
  print(paste0(i, " number of genes: ", nrow(obj),
               ", number of cells: ", ncol(obj)))
  obj_lst[[i]] <-
    PercentageFeatureSet(obj,
                         pattern = c("^MT-"),
                         col.name = "percent.mt")
  obj_lst[[i]]$lib.ident <- i
  normed_lst[[i]] <- NormalizeData(obj_lst[[i]],
                                normalization.method = "CLR")
  normed_lst[[i]] <- FindVariableFeatures(normed_lst[[i]],
                                          selection.method = "vst",
                                          nfeatures = 5000)
  # all.genes <- rownames(normed_lst[[i]])
  # normed_lst[[i]] <- ScaleData(normed_lst[[i]], features = all.genes)
}
scaled_lst <- normed_lst
for (i in 1:length(normed_lst)){
  print(sum(colSums(normed_lst[[i]]@assays$RNA@data))/ncol(normed_lst[[i]]@assays$RNA@data))
}

# integration and init downstream analysis ####
# find anchors
features <- SelectIntegrationFeatures(object.list = scaled_lst, 
                                      nfeatures = 5000,
                                      fvf.nfeatures = 5000)
anchors <- FindIntegrationAnchors(object.list = scaled_lst, 
                                  anchor.features = features)
# Retained 17601 anchors

# integrate
integrated <- IntegrateData(anchorset = anchors, 
                            normalization.method = "LogNormalize")

DefaultAssay(integrated)

# downstream analysis
integrated$group.ident <- NA
integrated$group.ident[integrated$lib.ident %in% c(1, 2, 3)] <- "group_2"
integrated$group.ident[integrated$lib.ident %in% c(4, 5, 6)] <- "group_8"
unique(integrated$group.ident)

integrated$time.ident <- NA
integrated$time.ident[integrated$lib.ident %in% c(1, 4)] <- "0hr"
integrated$time.ident[integrated$lib.ident %in% c(2, 5)] <- "1hr"
integrated$time.ident[integrated$lib.ident %in% c(3, 6)] <- "6hr"
unique(integrated$time.ident)

sum(integrated$lib.ident == "1") #g20: 14572
sum(integrated$lib.ident == "2") #g21: 13560
sum(integrated$lib.ident == "3") #g26: 14721
sum(integrated$lib.ident == "4") #g80: 13261
sum(integrated$lib.ident == "5") #g81: 12534
sum(integrated$lib.ident == "6") #g86: 14445

# cell line composition ####
times <- c("0", "1", "6")
lines_g2 <- c("CD_27", "CD_54")
lines_g8 <- c("CD_08", "CD_25", "CD_26")
g_2_0_bc_lst <- vector(mode = "list", length = 2)
g_2_1_bc_lst <- vector(mode = "list", length = 2)
g_2_6_bc_lst <- vector(mode = "list", length = 2)
g_8_0_bc_lst <- vector(mode = "list", length = 3)
g_8_1_bc_lst <- vector(mode = "list", length = 3)
g_8_6_bc_lst <- vector(mode = "list", length = 3)

# import cell ident barcodes
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
integrated$cell.line.ident <- "unmatched"
integrated$cell.line.ident[str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                                            c(g_2_0_bc_lst[[1]], g_2_1_bc_lst[[1]], g_2_6_bc_lst[[1]])] <- "CD_27"
integrated$cell.line.ident[str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                                            c(g_2_0_bc_lst[[2]], g_2_1_bc_lst[[2]], g_2_6_bc_lst[[2]])] <- "CD_54"

integrated$cell.line.ident[str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                                            c(g_8_0_bc_lst[[1]], g_8_1_bc_lst[[1]], g_8_6_bc_lst[[1]])] <- "CD_08"
integrated$cell.line.ident[str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                                            c(g_8_0_bc_lst[[2]], g_8_1_bc_lst[[2]], g_8_6_bc_lst[[2]])] <- "CD_25"
integrated$cell.line.ident[str_sub(integrated@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                                            c(g_8_0_bc_lst[[3]], g_8_1_bc_lst[[3]], g_8_6_bc_lst[[3]])] <- "CD_26"
sum(integrated$cell.line.ident == "unmatched") # 25173

human_only <- subset(integrated, subset = cell.line.ident != "unmatched")
# look at the cell line and time point composition of the dataset
# get the number of clusters, time points, lines
clusters <- levels(human_only$seurat_clusters)
times <- unique(human_only$time.ident)
lines <- unique(human_only$cell.line.ident)
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
      s <- sum((human_only$time.ident == t) & (human_only$cell.line.ident %in% l) & (human_only$seurat_clusters == c))
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


# clustering ####
human_only <- ScaleData(human_only)
human_only <- RunPCA(human_only, 
                     verbose = T, 
                     seed.use = 11)
human_only <- RunUMAP(human_only, 
                      reduction = "pca", 
                      dims = 1:30,
                      seed.use = 11)
human_only <- FindNeighbors(human_only, 
                            reduction = "pca", 
                            dims = 1:30)
human_only <- FindClusters(human_only, 
                           resolution = 0.5,
                           random.seed = 11)

DimPlot(human_only, 
        label = T) +
  theme(text = element_text(size = 12)) +
  NoLegend()

# check time point distribution 
human_only_dimplot <- DimPlot(human_only,
                              cols = c("cyan", "magenta", "yellow"),
                              group.by = "time.ident",
                              label = F) +
  ggtitle("colored by time") +
  theme(text = element_text(size = 12))
human_only_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
human_only_dimplot

# check library distribution
by_lib_dimplot <- DimPlot(human_only,
                          cols = "Set1",
                          group.by = "lib.ident",
                          label = F) +
  ggtitle("colored by library") +
  theme(text = element_text(size = 12))
by_lib_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
by_lib_dimplot


# assign cell types ####
type_markers <- c("GAD1", "GAD2", "SLC17A6",
                     "EBF1", # striatal
                     "ADCYAP1", "TCERG1L", "SEMA3E", # subcerebral
                     "CUX1", "BCL11B",  # cortical
                     "SST", # inhibitory
                     "SATB2", "CUX2", "NEFM", # excitatory
                     "VIM", "SOX2", "NES", #NPC
                     "MAP2", "DCX",
                     "SLC17A7", "SERTAD4", "FOXG1", # forebrain
                     "POU3F2", "LHX2" # general cortex
)


StackedVlnPlot(obj = human_only, features = type_markers) +
  coord_flip()

FeaturePlot(human_only, features = c("GAD1", "GAD2"))
FeaturePlot(human_only, features = c("NEFM"))
FeaturePlot(human_only, features = c("SEMA3E"))
FeaturePlot(human_only, features = c("FOXG1", "NES"))
FeaturePlot(human_only, features = c("VIM", "SOX2"))
DefaultAssay(human_only) <- "RNA"
FeaturePlot(human_only, features = c("SLC17A6", "SLC17A7"))
FeaturePlot(human_only, features = c("CUX2"))

# assign cluster identities based on the violin plots
new.cluster.ids <- 
  c("NEFM_pos_glut", "NEFM_neg_glut", "GABA", "GABA", "NEFM_pos_glut", 
    "NEFM_pos_glut", "NEFM_pos_glut", "GABA", "NEFM_neg_glut", "SEMA3E_pos_glut", 
    "subcerebral_neuron", "NPC", "NEFM_neg_glut", "NEFM_neg_glut", "GABA", 
    "SEMA3E_pos_glut", "NEFM_neg_glut", "forebrain_NPC", "NEFM_neg_glut", "NPC", 
    "SOX2+ NPC", "unknown", "unknown", "forebrain_NPC", "subcerebral_neuron", 
    "SST+ neuron", "SEMA3E_pos_glut", "cortical_NPC", "subcerebral_neuron")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(human_only$seurat_clusters))

names(new.cluster.ids) <- levels(human_only)
human_only <- RenameIdents(human_only, new.cluster.ids)

human_only$cell.type <- human_only@active.ident
unique(human_only$cell.type)

human_only$cell.type.for.plot <- str_replace_all(human_only$cell.type, 
                                             "_",
                                             " ")
human_only$cell.type.for.plot <- str_replace(human_only$cell.type.for.plot, 
                                             " pos",
                                             "+")
human_only$cell.type.for.plot <- str_replace(human_only$cell.type.for.plot, 
                                             " neg",
                                             "-")
unique(human_only$cell.type.for.plot)
Idents(human_only) <- "cell.type.for.plot"
DimPlot(human_only, 
        reduction = "umap", 
        label = TRUE, 
        repel = T, 
        pt.size = 0.3) +
  ggtitle("clusters labeled by cell type") +
  theme(text = element_text(size = 11)) +
  NoLegend()

sum(human_only@active.ident == "NPC") #3009
sum(human_only@active.ident == "NEFM- glut") #14145
sum(human_only@active.ident == "NEFM+ glut") #16322
sum(human_only@active.ident == "GABA") #13583
sum(human_only@active.ident == "SEMA3E+ glut") #4337

filtered_obj <- subset(human_only, 
                       cell.type %in% c("NEFM_pos_glut", "NEFM_neg_glut", "GABA", "NPC"))


# check early late response gene markers ####
zero <- subset(human_only, time.ident == "0hr")
one <- subset(human_only, time.ident == "1hr")
six <- subset(human_only, time.ident == "6hr")

early_late_markers <- c("FOS", "VGF", "BDNF", "GAPDH")

FeaturePlot(zero, features = early_late_markers)
FeaturePlot(one, features = early_late_markers)
FeaturePlot(six, features = early_late_markers)

filtered_obj$type.time.ident <- "others"
filtered_obj$type.time.for.plot <- "others"
for (ti in unique(filtered_obj$time.ident)){
  for (ty in unique(filtered_obj$cell.type)){
    fpty <- unique(filtered_obj$cell.type.for.plot[filtered_obj$time.ident == ti &
                                                       filtered_obj$cell.type == ty])
    fpname <- paste0(fpty, "-", ti)
    name <- paste0(ty, "_", ti)
    filtered_obj$type.time.ident[filtered_obj$time.ident == ti &
                                   filtered_obj$cell.type == ty] <- name
    filtered_obj$type.time.for.plot[filtered_obj$time.ident == ti &
                                   filtered_obj$cell.type == ty] <- fpname
  }
}
unique(filtered_obj$type.time.ident)
unique(filtered_obj$type.time.for.plot)

# differentially expressed genes ####
