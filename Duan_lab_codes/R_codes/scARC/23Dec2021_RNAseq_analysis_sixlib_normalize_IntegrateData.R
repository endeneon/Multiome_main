# Chuxuan Li 12/23/2021
# RNAseq data analysis by SCT-normalizing 6 libraries individually first, then
# use IntegrateData() in Seurat to integrate 6 libraries


# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)
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
transformed_lst <- vector(mode = "list", length = 6L)
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
  transformed_lst[[i]] <- SCTransform(obj_lst[[i]], 
                                      vars.to.regress = "percent.mt",
                                      method = "glmGamPoi",
                                      variable.features.n = 8000,
                                      seed.use = 115,
                                      verbose = T)
}


# integration and init downstream analysis ####
# find anchors
features <- SelectIntegrationFeatures(object.list = transformed_lst, 
                                      nfeatures = 5000,
                                      fvf.nfeatures = 5000)
transformed_lst_prepped <- PrepSCTIntegration(object.list = transformed_lst, 
                                              anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = transformed_lst_prepped, 
                                  normalization.method = "SCT",
                                  anchor.features = features)
# Retained 22442 anchors

# integrate
integrated <- IntegrateData(anchorset = anchors, 
                            normalization.method = "SCT")

DefaultAssay(integrated) <- "integrated"

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

# integrated <- SCTransform(integrated, 
#                           vars.to.regress = "percent.mt",
#                           method = "glmGamPoi",
#                           variable.features.n = 8000,
#                           seed.use = 115,
#                           verbose = T)
#
# do not do SCT again: causes time points to separate


# Run the standard workflow for visualization and clustering
integrated <- RunPCA(integrated, 
                     verbose = T, 
                     seed.use = 11)
integrated <- RunUMAP(integrated, 
                      reduction = "pca", 
                      dims = 1:30,
                      seed.use = 11)
integrated <- FindNeighbors(integrated, 
                            reduction = "pca", 
                            dims = 1:30)
integrated <- FindClusters(integrated, 
                           resolution = 0.5,
                           random.seed = 11)


DimPlot(integrated, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of independently normalized, \n Seurat-integrated data")

# check time point distribution 
integrated_dimplot <- DimPlot(integrated,
                              #cols = c("steelblue3", "mediumpurple4", "indianred3"),
                              cols = c("cyan", "magenta", "yellow"),
                              # cols = "Set2",
                              # cols = g_2_6_gex@meta.data$cell.line.ident,
                              group.by = "time.ident",
                              label = F) +
  ggtitle("independently normalized RNAseq data \n integrated by Seurat \n colored by time")
integrated_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
integrated_dimplot

# check library distribution
by_lib_dimplot <- DimPlot(integrated,
                          #cols = c("steelblue3", "mediumpurple4", "indianred3"),
                          cols = "Set1",
                          #cols = c("transparent", "transparent", "transparent", "transparent", "transparent", "steelblue"),
                          # cols = g_2_6_gex@meta.data$cell.line.ident,
                          group.by = "lib.ident",
                          label = F) +
  ggtitle("independently normalized RNAseq data \n integrated by Seurat \n colored by library")
by_lib_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
by_lib_dimplot

sum(integrated$lib.ident == "1") #g20: 14572
sum(integrated$lib.ident == "2") #g21: 13560
sum(integrated$lib.ident == "3") #g26: 14721
sum(integrated$lib.ident == "4") #g80: 13261
sum(integrated$lib.ident == "5") #g81: 12534
sum(integrated$lib.ident == "6") #g86: 14445

save("integrated", file = "normalize_by_6libs_integrated.RData")


# assign cell types ####

DefaultAssay(integrated) <- "SCT"
# plot stacked vlnplot
# subtype_markers <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
#                      "EBF1", "SERTAD4", # striatal
#                      "FOXG1", # forebrain 
#                      "TBR1", "FOXP2", "TLE4", # pallial glutamatergic
#                      "FEZF2", "ADCYAP1", "TCERG1L", "SEMA3E", # subcerebral
#                      "POU3F2", "CUX1", "BCL11B",  # cortical
#                      "LHX2", # general cortex
#                      #"EOMES", # hindbrain
#                      "NPY", "SST", "DLX2", "DLX5", # inhibitory
#                      "SATB2", "CBLN2", "CUX2", "NEFM", # excitatory
#                      "VIM", "SOX2", "NES" #NPC
# )
# 
# StackedVlnPlot(obj = integrated, 
#                features = subtype_markers) +
#   coord_flip()

trimmed_markers <- c("GAD1", "GAD2", "SLC17A6",
                     "EBF1", # striatal
                     "ADCYAP1", "TCERG1L", "SEMA3E", # subcerebral
                     "CUX1", "BCL11B",  # cortical
                     "SST", # inhibitory
                     "SATB2", "CUX2", "NEFM", # excitatory
                     "VIM", "SOX2", "NES", #NPC
                     "MAP2", "DCX"
)

other_markers <- c("SLC17A7", "SERTAD4", "FOXG1", # forebrain
                   "POU3F2", "LHX2" # general cortex
)

StackedVlnPlot(obj = integrated, features = trimmed_markers) +
  coord_flip()
StackedVlnPlot(obj = integrated, features = other_markers) +
  coord_flip()

StackedVlnPlot(obj = integrated, features = "NEFM") +
  coord_flip() 


FeaturePlot(integrated, features = c("GAD1", "GAD2"))
FeaturePlot(integrated, features = c("SLC17A6", "SLC17A7"), max.cutoff = 15)
FeaturePlot(integrated, features = c("CUX2"), max.cutoff = 15)
FeaturePlot(integrated, features = c("NEFM"))
FeaturePlot(integrated, features = c("SEMA3E"), max.cutoff = 20)
FeaturePlot(integrated, features = c("FOXG1", "NES"), max.cutoff = 10)
FeaturePlot(integrated, features = c("VIM", "SOX2"))


early_late_markers <- c("EGR1", "FOS", "FOSB", "NPAS4", # early response
                        "VGF", "BDNF", "NRN1", "PNOC" # late response
)
StackedVlnPlot(obj = integrated, features = early_late_markers) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


# try a dendrogram plot
integrated <- BuildClusterTree(object = integrated, 
                               assay = "SCT",
                               features = trimmed_markers,
                               verbose = T)

PlotClusterTree(integrated, 
                direction = "rightwards",
                use.edge.length = F)

# assign cluster identities based on the violin plots
new.cluster.ids <- 
  c("NEFM_pos_glut", "forebrain_NPC", "GABA", "GABA", "NEFM_neg_glut", 
    "NEFM_neg_glut", "NEFM_pos_glut", "unknown", "NEFM_pos_glut", "unknown", 
    "NEFM_pos_glut", "GABA", "SEMA3E_pos_glut", "NPC", "unknown", 
    "NEFM_neg_glut", "GABA", "SEMA3E_pos_glut", "NPC", "NEFM_neg_glut", 
    "forebrain_NPC", "GABA", "unknown", "NEFM_pos_glut", "forebrain_NPC", 
    "NEFM_pos_glut", "NEFM_neg_glut")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated$seurat_clusters))


names(new.cluster.ids) <- levels(integrated)

RNAseq_integrated_labeled <- RenameIdents(integrated, new.cluster.ids)


DimPlot(RNAseq_integrated_labeled, 
        reduction = "umap", 
        label = TRUE, 
        repel = T, 
        pt.size = 0.3) +
  ggtitle("independently normalized RNAseq data \n integrated by Seurat, labeled by cell type")

sum(RNAseq_integrated_labeled@active.ident == "NPC") #5068
sum(RNAseq_integrated_labeled@active.ident == "forebrain_NPC") #7767
sum(RNAseq_integrated_labeled@active.ident == "NEFM_neg_glut") #14948
sum(RNAseq_integrated_labeled@active.ident == "NEFM_pos_glut") #19519
sum(RNAseq_integrated_labeled@active.ident == "GABA") #18752
sum(RNAseq_integrated_labeled@active.ident == "SEMA3E_pos_glut") #5510
sum(RNAseq_integrated_labeled@active.ident == "unknown") #11529


# differentially expressed genes ####

RNAseq_integrated_labeled$broad.cell.type <- NA
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled@active.ident == "SEMA3E_pos_glut"] <- "SEMA3E_pos_glut"
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled@active.ident == "forebrain_NPC"] <- "forebrain_NPC"
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled@active.ident == "NEFM_pos_glut"] <- "NEFM_pos_glut"
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled@active.ident == "NEFM_neg_glut"] <- "NEFM_neg_glut"
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled@active.ident == "NPC"] <- "NPC"
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled@active.ident == "GABA"] <- "GABA"
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled@active.ident == "unknown"] <- "unknown"
unique(RNAseq_integrated_labeled$broad.cell.type)

RNAseq_integrated_labeled$spec.cell.type <- "others"
RNAseq_integrated_labeled$spec.cell.type[RNAseq_integrated_labeled$broad.cell.type == "NEFM_pos_glut"] <- "NEFM_pos_glut"
RNAseq_integrated_labeled$spec.cell.type[RNAseq_integrated_labeled$broad.cell.type == "NEFM_neg_glut"] <- "NEFM_neg_glut"
RNAseq_integrated_labeled$spec.cell.type[RNAseq_integrated_labeled$broad.cell.type == "NPC"] <- "NPC"
RNAseq_integrated_labeled$spec.cell.type[RNAseq_integrated_labeled$broad.cell.type == "GABA"] <- "GABA"
unique(RNAseq_integrated_labeled$spec.cell.type)

save("RNAseq_integrated_labeled", file = "nfeature5000_labeled_obj.RData")


# use new method to plot dotplot

DefaultAssay(RNAseq_integrated_labeled) <- "integrated"
unique(Idents(RNAseq_integrated_labeled))
Idents(RNAseq_integrated_labeled) <- "spec.cell.type"
unique(Idents(RNAseq_integrated_labeled))


j <- 0
# calc hour 1 vs 0
for (i in sort(unique(RNAseq_integrated_labeled$spec.cell.type))) {
  print(i)
  j = j + 1
  if (j == 1) {
    print("first cell type")
    df_1vs0 <- FindMarkers(object = RNAseq_integrated_labeled,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0,
                           test.use = "MAST")
    df_1vs0$gene_symbol <- rownames(df_1vs0)
    df_1vs0$cell_type <- i
  } else {
    print("others")
    df_to_append <- FindMarkers(object = RNAseq_integrated_labeled,
                                ident.1 = "1hr",
                                group.by = 'time.ident',
                                subset.ident = i,
                                ident.2 = "0hr",
                                min.pct = 0.0,
                                logfc.threshold = 0.0,
                                test.use = "MAST")
    df_to_append$gene_symbol <- rownames(df_to_append)
    df_to_append$cell_type <- i
    df_1vs0 <- rbind(df_1vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
    }
}

j <- 0
# calculate 6hr vs 0hr
for (i in sort(unique(RNAseq_integrated_labeled$spec.cell.type))) {
  print(i)
  j = j + 1
  if (j == 1) {
    print("first cell type")
    df_6vs0 <- FindMarkers(object = RNAseq_integrated_labeled,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0,
                           test.use = "MAST")
    df_6vs0$gene_symbol <- rownames(df_6vs0)
    df_6vs0$cell_type <- i
  } else {
    print("others")
    df_to_append <- FindMarkers(object = RNAseq_integrated_labeled,
                                ident.1 = "6hr",
                                group.by = 'time.ident',
                                subset.ident = i,
                                ident.2 = "0hr",
                                min.pct = 0.0,
                                logfc.threshold = 0.0,
                                test.use = "MAST")
    df_to_append$gene_symbol <- rownames(df_to_append)
    df_to_append$cell_type <- i
    df_6vs0 <- rbind(df_6vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
  }
}


# make a bubble plot
gene_list <- c("BDNF", "FOS", "NPAS4", "EGR1")
main_types <- main_types <- c("GABA", "NEFM_pos_glut", "NEFM_neg_glut")
df_to_plot <- rbind(df_1vs0[df_1vs0$gene_symbol %in% gene_list & 
                              df_1vs0$cell_type %in% main_types, ],
                    df_6vs0[df_6vs0$gene_symbol %in% gene_list & 
                              df_1vs0$cell_type %in% main_types, ])
df_to_plot$source <- c(rep_len("1v0hr", length.out = 12),
                       rep_len("6v0hr", length.out = 12))


ggplot(df_to_plot,
       aes(x = as.factor(source),
           y = as.factor(cell_type), 
           size = pct.1 * 100,
           fill = avg_log2FC * pct.1 / pct.2)) +
  geom_point(shape = 21) +
  facet_grid(cols = vars(gene_symbol)) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  theme_classic() #+
#coord_flip()

# find variable markers for each time point
GABA_subset <- subset(RNAseq_integrated_labeled, subset = broad.cell.type == "GABA")
NPCM_pos_glut_subset <- subset(RNAseq_integrated_labeled, subset = broad.cell.type == "NEFM_pos_glut")
NPCM_neg_glut_subset <- subset(RNAseq_integrated_labeled, subset = broad.cell.type == "NEFM_neg_glut")
NPC_subset <- subset(RNAseq_integrated_labeled, subset = broad.cell.type == "NPC")
forebrain_NPC_subset <- subset(RNAseq_integrated_labeled, subset = broad.cell.type == "forebrain_NPC")


for (t in c(GABA_subset, NPCM_pos_glut_subset, NPCM_neg_glut_subset, NPC_subset, forebrain_NPC_subset)){
  t_name <- unique(t$broad.cell.type)
  print(t_name)
  
  # differentially expressed genes comparing one big cell type with all others
  Idents(t) <- "time.ident"
  print(DefaultAssay(t))
  zero_one <- FindMarkers(t, 
                          ident.1 = "1hr", 
                          logfc.threshold = 0, 
                          min.pct = 0,
                          ident.2 = "0hr", 
                          group.by = "time.ident",
                          test.use = "MAST",
                          verbose = TRUE)
  zero_six <- FindMarkers(t, 
                          ident.1 = "6hr", 
                          logfc.threshold = 0,
                          min.pct = 0,
                          ident.2 = "0hr", 
                          group.by = "time.ident",
                          test.use = "MAST",
                          verbose = TRUE)
  
  # output the list of variable markers 
  write.table(zero_one, 
              file = paste0(t_name, "_DE_markers_1v0hr.csv"),
              col.names = T,
              row.names = T,
              sep = ",",
              quote = F)
  print(nrow(zero_one))
  write.table(zero_six, 
              file = paste0(t_name, "_DE_markers_6v0hr.csv"),
              col.names = T,
              row.names = T,
              sep = ",",
              quote = F)
  print(nrow(zero_six))
  
}


# examine cell line composition ####

# assign cell line identities
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
RNAseq_integrated_labeled$cell.line.ident <- NA
RNAseq_integrated_labeled$cell.line.ident[str_sub(RNAseq_integrated_labeled@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                                            c(g_2_0_bc_lst[[1]], g_2_1_bc_lst[[1]], g_2_6_bc_lst[[1]])] <- "CD_27"
RNAseq_integrated_labeled$cell.line.ident[str_sub(RNAseq_integrated_labeled@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                                            c(g_2_0_bc_lst[[2]], g_2_1_bc_lst[[2]], g_2_6_bc_lst[[2]])] <- "CD_54"

RNAseq_integrated_labeled$cell.line.ident[str_sub(RNAseq_integrated_labeled@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                                            c(g_8_0_bc_lst[[1]], g_8_1_bc_lst[[1]], g_8_6_bc_lst[[1]])] <- "CD_08"
RNAseq_integrated_labeled$cell.line.ident[str_sub(RNAseq_integrated_labeled@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                                            c(g_8_0_bc_lst[[2]], g_8_1_bc_lst[[2]], g_8_6_bc_lst[[2]])] <- "CD_25"
RNAseq_integrated_labeled$cell.line.ident[str_sub(RNAseq_integrated_labeled@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                                            c(g_8_0_bc_lst[[3]], g_8_1_bc_lst[[3]], g_8_6_bc_lst[[3]])] <- "CD_26"
sum(is.na(RNAseq_integrated_labeled$cell.line.ident))


# look at the cell line and time point composition of the dataset
# get the number of clusters, time points, lines
clusters <- levels(RNAseq_integrated_labeled$seurat_clusters)
times <- unique(RNAseq_integrated_labeled$time.ident)
lines <- unique(RNAseq_integrated_labeled$cell.line.ident)
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
      s <- sum((RNAseq_integrated_labeled$time.ident == t) & (RNAseq_integrated_labeled$cell.line.ident %in% l) & (RNAseq_integrated_labeled$seurat_clusters == c))
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




# Pseudobulk de analysis ####
RNAseq_integrated_labeled$cell.line.ident[is.na(RNAseq_integrated_labeled$cell.line.ident)] <- "unmatched"
celllines <- unique(RNAseq_integrated_labeled$cell.line.ident)
celllines <- celllines[celllines != "unmatched"]

# first add idents indicating both cell line and time
RNAseq_integrated_labeled$cell.line_time.point <- NA
for (cl in celllines){
  print(cl)
  RNAseq_integrated_labeled$cell.line_time.point[RNAseq_integrated_labeled$time.ident == "0hr" &
                                                   RNAseq_integrated_labeled$cell.line.ident == cl] <- paste0(cl, "_0hr")
  RNAseq_integrated_labeled$cell.line_time.point[RNAseq_integrated_labeled$time.ident == "1hr" &
                                                   RNAseq_integrated_labeled$cell.line.ident == cl] <- paste0(cl, "_1hr")
  RNAseq_integrated_labeled$cell.line_time.point[RNAseq_integrated_labeled$time.ident == "6hr" &
                                                   RNAseq_integrated_labeled$cell.line.ident == cl] <- paste0(cl, "_6hr")
}
unique(RNAseq_integrated_labeled$cell.line_time.point)

line_time <- sort(unique(RNAseq_integrated_labeled$cell.line_time.point))

GABA <- subset(RNAseq_integrated_labeled, subset = broad.cell.type == "GABA")
glut <- subset(RNAseq_integrated_labeled, broad.cell.type %in% c("NEFM_pos_glut", 
                                                                 "NEFM_neg_glut"))

genelist <- rownames(GABA@assays$RNA@data)
filler <- rep(0, length(genelist))
GABA_pseudobulk <- data.frame(CD_08_0hr = filler, CD_08_1hr = filler, CD_08_6hr = filler,
                            CD_25_0hr = filler, CD_25_1hr = filler, CD_25_6hr = filler,
                            CD_26_0hr = filler, CD_26_1hr = filler, CD_26_6hr = filler,
                            CD_27_0hr = filler, CD_27_1hr = filler, CD_27_6hr = filler,
                            CD_54_0hr = filler, CD_54_1hr = filler, CD_54_6hr = filler, row.names = genelist)
genelist <- rownames(glut@assays$RNA@data)
filler <- rep(0, length(genelist))
glut_pseudobulk <- data.frame(CD_08_0hr = filler, CD_08_1hr = filler, CD_08_6hr = filler,
                              CD_25_0hr = filler, CD_25_1hr = filler, CD_25_6hr = filler,
                              CD_26_0hr = filler, CD_26_1hr = filler, CD_26_6hr = filler,
                              CD_27_0hr = filler, CD_27_1hr = filler, CD_27_6hr = filler,
                              CD_54_0hr = filler, CD_54_1hr = filler, CD_54_6hr = filler, row.names = genelist)


for (i in 1:length(line_time)){
  lt <- line_time[i]
  subobj <- subset(GABA, subset = cell.line_time.point == lt)
  GABA_pseudobulk[, i] <- as.array(rowSums(subobj@assays$RNA@data))
  subobj <- subset(glut, subset = cell.line_time.point == lt)
  glut_pseudobulk[, i] <- as.array(rowSums(subobj@assays$RNA@data))
}

save("GABA_pseudobulk", "glut_pseudobulk", file = "data_frames_for_pseudobulk_new.RData")
save("RNAseq_integrated_labeled", file = "27Jan2022_5line_added_cellline_idents")

# other analysis ####
DotPlot(RNAseq_integrated_labeled, 
        features = early_late_markers, 
        cols = c("steelblue3", "mediumpurple3", "indianred3"), 
        dot.scale = 5, 
        split.by = "time.ident") +
  theme_bw()+
  # RotatedAxis() +
  coord_flip() 

# plot dotplot split by time but only ext/inh
DefaultAssay(RNAseq_integrated_labeled) <- "SCT"
grant_markers <- c("FOS", "BDNF", "NPAS4", "SP4", "TCF4")

two_celltypes_only <- subset(x = RNAseq_integrated_labeled,
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

# feature plot for some genes 
RNAseq_integrated_labeled$cell.line.ident[is.na(RNAseq_integrated_labeled$cell.line.ident)] <- "NA"
RNAseq_integrated_labeled <- subset(RNAseq_integrated_labeled, cell.line.ident != "NA")
DefaultAssay(RNAseq_integrated_labeled) <- "SCT"
RNAseq_integrated_labeled <- ScaleData(RNAseq_integrated_labeled)
FeaturePlot(subset(RNAseq_integrated_labeled, broad.cell.type == "GABA"), 
            features = c("ACADM", "HSD17B4", "ACAT1", "CNR1"),
            cols = c("transparent", "blue"), max.cutoff = 7.5,
            split.by = "time.ident", keep.scale = "all") &
  theme(legend.position = c(0.75, 0.5))
FeaturePlot(subset(RNAseq_integrated_labeled, broad.cell.type == "GABA"), slot = "scale.data",
            features = c("ACADM", "HSD17B4", "ACAT1", "CNR1"),
            cols = c("transparent", "darkblue"), max.cutoff = 7.5,
            split.by = "time.ident", keep.scale = "all") & 
  theme(legend.position = c(0.75, 0.5))
