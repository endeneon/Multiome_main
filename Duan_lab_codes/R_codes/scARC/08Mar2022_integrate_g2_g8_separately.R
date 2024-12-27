# Chuxuan Li 03/08/2022
# Check if only integrate group 2 and group 8 within themselves, without integrating
#these two groups, would expression levels for marker genes and differential expression
#show expected patterns

# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)

library(ggplot2)
library(RColorBrewer)
library(graphics)

library(stringr)
library(readr)

library(future)

plan("multisession", workers = 1)

# load data ####
# note this h5 file contains both atac-seq and gex information
h5list <- list.files("/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only", 
                     pattern = "filtered_feature_bc_matrix.h5", 
                     recursive = T) 
h5list <- paste0("/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/", h5list)
objlist <- vector(mode = "list", length = length(h5list))
for (i in 1:length(h5list)){
  h5file <- Read10X_h5(filename = h5list[[i]])
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(h5list[[i]], "[0-8]_[0-6]"))
  obj$lib.ident <- paste0(str_extract(h5list[[i]], "[2|8]_[0|1|6]"))
  print(paste0(str_extract(h5list[[i]], "[2|8]_[0|1|6]"), 
               "number of genes: ", nrow(obj), 
               ", number of cells: ", ncol(obj)))
  obj <- PercentageFeatureSet(obj,
                         pattern = c("^MT-"),
                         col.name = "percent.mt")
  objlist[[i]] <- SCTransform(obj, 
                              vars.to.regress = "percent.mt",
                              seed.use = 10,
                              verbose = T)
}

# integrate g2 ####
features_g2 <- SelectIntegrationFeatures(object.list = objlist[1:3], 
                                      nfeatures = 6000,
                                      assay = rep("SCT", 3))
list_g2 <- PrepSCTIntegration(object.list = objlist[1:3], 
                              anchor.features = features_g2,
                              assay = "SCT")
anchors_g2 <- FindIntegrationAnchors(object.list = list_g2, 
                                  normalization.method = "SCT",
                                  anchor.features = features_g2,
                                  dims = 1:20,
                                  reduction = "cca", 
                                  verbose = T)
# Retained 22596 anchors
integrated_g2 <- IntegrateData(anchorset = anchors_g2)

# integrate g8 ####
features_g8 <- SelectIntegrationFeatures(object.list = objlist[4:6], 
                                         nfeatures = 6000,
                                         assay = rep("SCT", 3))
list_g8 <- PrepSCTIntegration(object.list = objlist[4:6], 
                              anchor.features = features_g8,
                              assay = "SCT")
anchors_g8 <- FindIntegrationAnchors(object.list = list_g8, 
                                     normalization.method = "SCT",
                                     anchor.features = features_g8,
                                     dims = 1:20,
                                     reduction = "cca", 
                                     verbose = T)
# Retained _ anchors
integrated_g8 <- IntegrateData(anchorset = anchors_g8)

save("integrated_g2", "integrated_g8", file = "integrated_objects_g2_g8.RData")

# assign time idents
unique(integrated_g2$lib.ident)
integrated_g2$time.ident <- ""
for (i in unique(integrated_g2$lib.ident)){
  integrated_g2$time.ident[integrated_g2$lib.ident == i] <- paste0(str_extract(i, "[0|1|6]"),
                                                                   "hr")
}
unique(integrated_g2$time.ident)

integrated_g8$time.ident <- ""
for (i in unique(integrated_g8$lib.ident)){
  integrated_g8$time.ident[integrated_g8$lib.ident == i] <- paste0(str_extract(i, "[0|1|6]"),
                                                                   "hr")
}
unique(integrated_g8$time.ident)

# dimensional reduction and clustering ####
standard_clustering <- function(obj, seed = 100){
  DefaultAssay(obj) <- "integrated"
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, 
                verbose = T, 
                seed.use = seed)
  obj <- RunUMAP(obj, 
                 reduction = "pca", 
                 dims = 1:30,
                 seed.use = seed)
  obj <- FindNeighbors(obj, 
                       reduction = "pca", 
                       dims = 1:30)
  obj <- FindClusters(obj, 
                      resolution = 0.5,
                      random.seed = seed)
  return(obj)
}
integrated_g2 <- standard_clustering(integrated_g2)
integrated_g8 <- standard_clustering(integrated_g8)

# g2 ####
DimPlot(integrated_g2, 
        label = T) +
  NoLegend() +
  ggtitle("group 2")

# check time point distribution 
time_dimplot <- DimPlot(integrated_g2,
                        cols = c("cyan", "magenta", "yellow"),
                        group.by = "time.ident",
                        label = F) +
  ggtitle("group 2 colored by time")
time_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
time_dimplot

# check library distribution
by_lib_dimplot <- DimPlot(integrated_g2,
                          cols = "Set1",
                          group.by = "lib.ident",
                          label = F) +
  ggtitle("group 2 colored by library")
by_lib_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
by_lib_dimplot

# cell type assignment
DefaultAssay(integrated_g2) <- "SCT"
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
                     "VIM", "SOX2", "NES" #NPC
)

StackedVlnPlot(obj = integrated_g2,
               features = subtype_markers) +
  coord_flip()

FeaturePlot(integrated_g2, features = c("GAD1", "GAD2"))
FeaturePlot(integrated_g2, features = c("SLC17A6", "SLC17A7"), max.cutoff = 15)
FeaturePlot(integrated_g2, features = c("NEFM"))
FeaturePlot(integrated_g2, features = c("FOXG1", "NES"), max.cutoff = 10)
FeaturePlot(integrated_g2, features = c("VIM", "SOX2"))

FeaturePlot(integrated_g2, features = c("FOS", "BDNF"), split.by = "time.ident")

# assign cluster identities based on the violin plots
new.cluster.ids <- 
  c("NEFM_neg_glut", "NEFM_pos_glut", "forebrain_NPC", "GABA", "NEFM_neg_glut", 
    "SEMA3E_pos_glut", "GABA", "NEFM_pos_glut", "unknown", "GABA", 
    "unknown", "unknown", "NPC", "GABA", "SEMA3E_pos_glut", 
    "SST_pos_glut", "NEFM_neg_glut", "NPC", "cortical_NPC", "cortical_NPC", 
    "NPC", "GABA", "unknown", "unknown")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated_g2$seurat_clusters))
names(new.cluster.ids) <- levels(integrated_g2)
labeled_g2 <- RenameIdents(integrated_g2, new.cluster.ids)

DimPlot(labeled_g2, 
        reduction = "umap", 
        label = TRUE, 
        repel = T, 
        pt.size = 0.3) +
  ggtitle("group 2 labeled by cell type")

sum(labeled_g2@active.ident == "NPC") #2394
sum(labeled_g2@active.ident == "NEFM_neg_glut") #9761
sum(labeled_g2@active.ident == "NEFM_pos_glut") #6687
sum(labeled_g2@active.ident == "GABA") #9263

labeled_g2$cell.type <- labeled_g2@active.ident
unique(labeled_g2$cell.type)

# g8 ####
DimPlot(integrated_g8, 
        label = T) +
  NoLegend() +
  ggtitle("group 8")

# check time point distribution 
time_dimplot <- DimPlot(integrated_g8,
                        cols = c("cyan", "magenta", "yellow"),
                        group.by = "time.ident",
                        label = F) +
  ggtitle("group 8 colored by time")
time_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
time_dimplot

# check library distribution
by_lib_dimplot <- DimPlot(integrated_g8,
                          cols = "Set1",
                          group.by = "lib.ident",
                          label = F) +
  ggtitle("group 8 colored by library")
by_lib_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
by_lib_dimplot

# cell type assignment
DefaultAssay(integrated_g8) <- "SCT"
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
                     "VIM", "SOX2", "NES" #NPC
)

StackedVlnPlot(obj = integrated_g8,
               features = subtype_markers) +
  coord_flip()

FeaturePlot(integrated_g8, features = c("GAD1", "GAD2"))
FeaturePlot(integrated_g8, features = c("SLC17A6", "SLC17A7"), max.cutoff = 15)
FeaturePlot(integrated_g8, features = c("NEFM"))
FeaturePlot(integrated_g8, features = c("FOXG1", "NES"), max.cutoff = 10)
FeaturePlot(integrated_g8, features = c("VIM", "SOX2"))

# assign cluster identities based on the violin plots
new.cluster.ids <- 
  c("NEFM_pos_glut", "NEFM_neg_glut", "GABA", "GABA", "GABA", 
    "SST_pos_glut", "NPC", "SEMA3E_pos_glut", "unknown", "NEFM_pos_glut", 
    "NEFM_pos_glut", "unknown", "cortical_NPC", "unknown", "NEFM_neg_glut", 
    "NEFM_pos_glut", "NEFM_neg_glut", "unknown", "cortical_NPC", "GABA", 
    "unknown", "unknown")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated_g8$seurat_clusters))
names(new.cluster.ids) <- levels(integrated_g8)
labeled_g8 <- RenameIdents(integrated_g8, new.cluster.ids)

DimPlot(labeled_g8, 
        reduction = "umap", 
        label = TRUE, 
        repel = T, 
        pt.size = 0.3) +
  ggtitle("group 8 labeled by cell type")

sum(labeled_g8@active.ident == "NPC") #2071
sum(labeled_g8@active.ident == "NEFM_neg_glut") #5996
sum(labeled_g8@active.ident == "NEFM_pos_glut") #10640
sum(labeled_g8@active.ident == "GABA") #9035

labeled_g8$cell.type <- labeled_g8@active.ident
unique(labeled_g8$cell.type)

save("labeled_g2", "labeled_g8", file = "labeled_objs_g2_g8.RData")


# only keep 4 cell types in de analysis
fc_g2 <- subset(labeled_g2, cell.type %in% c("NEFM_neg_glut", "NEFM_pos_glut",
                                             "GABA", "NPC"))
fc_g8 <- subset(labeled_g8, cell.type %in% c("NEFM_neg_glut", "NEFM_pos_glut",
                                             "GABA", "NPC"))
# assign cell line ####
# assign cell line identitiesfc_g2
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
fc_g2$cell.line.ident <- "unmatched"
fc_g2$cell.line.ident[str_sub(fc_g2@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                        c(g_2_0_bc_lst[[1]], g_2_1_bc_lst[[1]], g_2_6_bc_lst[[1]])] <- "CD_27"
fc_g2$cell.line.ident[str_sub(fc_g2@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                        c(g_2_0_bc_lst[[2]], g_2_1_bc_lst[[2]], g_2_6_bc_lst[[2]])] <- "CD_54"

fc_g8$cell.line.ident <- "unmatched"
fc_g8$cell.line.ident[str_sub(fc_g8@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                        c(g_8_0_bc_lst[[1]], g_8_1_bc_lst[[1]], g_8_6_bc_lst[[1]])] <- "CD_08"
fc_g8$cell.line.ident[str_sub(fc_g8@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                        c(g_8_0_bc_lst[[2]], g_8_1_bc_lst[[2]], g_8_6_bc_lst[[2]])] <- "CD_25"
fc_g8$cell.line.ident[str_sub(fc_g8@assays$RNA@counts@Dimnames[[2]], end = -5L) %in% 
                        c(g_8_0_bc_lst[[3]], g_8_1_bc_lst[[3]], g_8_6_bc_lst[[3]])] <- "CD_26"

sum(fc_g2$cell.line.ident == "unmatched")
sum(fc_g8$cell.line.ident == "unmatched")
unique(fc_g2$cell.line.ident)
unique(fc_g8$cell.line.ident)

fc_g2 <- subset(fc_g2, cell.line.ident != "unmatched")
fc_g8 <- subset(fc_g8, cell.line.ident != "unmatched")
save("fc_g2", "fc_g8", file = "labeled_four_celltype_objs.RData")

CD_27 <- fc_g2 <- subset(fc_g2, cell.line.ident == "CD_27")
CD_54 <- fc_g2 <- subset(fc_g2, cell.line.ident == "CD_54")
CD_08 <- fc_g8 <- subset(fc_g8, cell.line.ident != "CD_08")
CD_25 <- fc_g8 <- subset(fc_g8, cell.line.ident != "CD_25")
CD_26 <- fc_g8 <- subset(fc_g8, cell.line.ident != "CD_26")

# differential expression analysis ####
Idents(fc_g2)
DefaultAssay(fc_g2) <- "integrated"
DotPlot(fc_g2, 
        features = c("BDNF", "VGF", "NPAS4", "FOS"), 
        cols = c("steelblue3", "steelblue3", "steelblue3"), 
        dot.scale = 5, 
        split.by = "time.ident") +
  theme_bw() +
  theme(text = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
  coord_flip() +
  ggtitle("group 2")


DefaultAssay(fc_g2) <- "RNA"
Idents(fc_g2)
fc_g2 <- ScaleData(fc_g2)
j <- 0
# calc hour 1 vs 0
for (i in sort(unique(fc_g2$cell.type))) {
  print(i)
  j = j + 1
  if (j == 1) {
    print("first cell type")
    df_1vs0 <- FindMarkers(object = fc_g2, ident.1 = "1hr", ident.2 = "0hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           min.pct = 0.0, logfc.threshold = 0.0,
                           features = c("BDNF", "VGF", "NPAS4", "FOS"),
                           test.use = "MAST")
    df_1vs0$gene_symbol <- rownames(df_1vs0)
    df_1vs0$cell_type <- i
  } else {
    print("others")
    df_to_append <- FindMarkers(object = fc_g2, ident.1 = "1hr",ident.2 = "0hr",
                                group.by = 'time.ident',
                                subset.ident = i,
                                min.pct = 0.0, logfc.threshold = 0.0,
                                features = c("BDNF", "VGF", "NPAS4", "FOS"),
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
for (i in sort(unique(fc_g2$cell.type))) {
  print(i)
  j = j + 1
  if (j == 1) {
    print("first cell type")
    df_6vs0 <- FindMarkers(object = fc_g2,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0,
                           features = c("BDNF", "VGF", "NPAS4", "FOS"),
                           test.use = "MAST")
    df_6vs0$gene_symbol <- rownames(df_6vs0)
    df_6vs0$cell_type <- i
  } else {
    print("others")
    df_to_append <- FindMarkers(object = fc_g2,
                                ident.1 = "6hr",
                                group.by = 'time.ident',
                                subset.ident = i,
                                ident.2 = "0hr",
                                min.pct = 0.0,
                                logfc.threshold = 0.0,
                                features = c("BDNF", "VGF", "NPAS4", "FOS"),
                                test.use = "MAST")
    df_to_append$gene_symbol <- rownames(df_to_append)
    df_to_append$cell_type <- i
    df_6vs0 <- rbind(df_6vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
  }
}
df_to_plot <- rbind(df_1vs0,
                    df_6vs0)
df_to_plot$source <- c(rep_len("1v0hr", length.out = 16),
                       rep_len("6v0hr", length.out = 16))

ggplot(df_to_plot,
       aes(x = as.factor(source),
           y = as.factor(cell_type), 
           size = pct.1 * 100,
           fill = avg_log2FC)) +
  geom_point(shape = 21) +
  facet_grid(cols = vars(gene_symbol)) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  theme_classic() 
