# Chuxuan Li 003/07/2022
# use RNAseq analysis v3 on pilot data to test the method

# init ####
library(Seurat)
library(sctransform)

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)
library(patchwork)
library(cowplot)
library(gplots)

library(future)
library(readr)
library(stringr)

# load data ####
h5list <- list.files(path = "/home/cli/NVME/scRNASeq_1st_test/Duan_Project_017-GEX", 
           pattern = "filtered_feature_bc_matrix.h5", recursive = T)
h5list <- paste0("/home/cli/NVME/scRNASeq_1st_test/Duan_Project_017-GEX/", h5list)
objlist <- vector(mode = "list", length = length(h5list))
transformed_lst <- vector(mode = "list", length = length(h5list))

for (i in 1:length(objlist)){
  lib <- str_extract(h5list[[i]], "SID[0-6]")
  print(lib)
  obj <- Read10X_h5(h5list[[i]])
  obj <- CreateSeuratObject(counts = obj,
                            project = lib)
  print(paste0(i, " number of genes: ", nrow(obj), 
               ", number of cells: ", ncol(obj)))
  objlist[[i]] <- 
    PercentageFeatureSet(obj,
                         pattern = "^MT-|^Mt-",
                         col.name = "percent.mt")

  transformed_lst[[i]] <- SCTransform(objlist[[i]], 
                                      vars.to.regress = "percent.mt",
                                      verbose = T)
  transformed_lst[[i]]$lib.ident <- lib
}


features <- SelectIntegrationFeatures(object.list = transformed_lst, 
                                      nfeatures = 6000,
                                      assay = rep("SCT", 3)
                                      )
transformed_lst_prepped <- PrepSCTIntegration(object.list = transformed_lst, 
                                              anchor.features = features,
                                              assay = "SCT")
anchors <- FindIntegrationAnchors(object.list = transformed_lst_prepped, 
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  dims = 1:20,
                                  reduction = "cca", 
                                  verbose = T)
# Retained 21300 anchors
#all.genes <- obj@assays$RNA@counts@Dimnames[[1]]
# integrate
integrated <- IntegrateData(anchorset = anchors)

DefaultAssay(integrated) <- "integrated"
unique(integrated$lib.ident)

# downstream analysis
integrated$time.ident <- NA
integrated$time.ident[integrated$lib.ident %in% "SID0"] <- "0hr"
integrated$time.ident[integrated$lib.ident %in% "SID1"] <- "1hr"
integrated$time.ident[integrated$lib.ident %in% "SID6"] <- "6hr"
unique(integrated$time.ident)

unique(RNASeq.integrated$indiv.identity)
CD_05 <- subset(RNASeq.integrated, indiv.identity == "Line 05")
CD_14 <- subset(RNASeq.integrated, indiv.identity == "Line 14")
CD_15 <- subset(RNASeq.integrated, indiv.identity == "Line 15")

integrated$cell.line <- "unmatched"
integrated$cell.line[integrated@assays$RNA@counts@Dimnames[[2]] %in% CD_05@assays$RNA@counts@Dimnames[[2]]] <- "CD_05"
integrated$cell.line[integrated@assays$RNA@counts@Dimnames[[2]] %in% CD_14@assays$RNA@counts@Dimnames[[2]]] <- "CD_14"
integrated$cell.line[integrated@assays$RNA@counts@Dimnames[[2]] %in% CD_15@assays$RNA@counts@Dimnames[[2]]] <- "CD_15"
unique(integrated$cell.line)


# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, 
                     npcs = 30,
                     verbose = T, 
                     seed.use = 10)
integrated <- RunUMAP(integrated, 
                      reduction = "pca", 
                      dims = 1:30,
                      seed.use = 10)
integrated <- FindNeighbors(integrated, 
                            reduction = "pca", 
                            dims = 1:30)
integrated <- FindClusters(integrated, 
                           resolution = 0.25,
                           random.seed = 10)


DimPlot(integrated, 
        label = T) +
  NoLegend() +
  ggtitle("pilot data")

# check time point distribution 
integrated_dimplot <- DimPlot(integrated,
                              #cols = c("steelblue3", "mediumpurple4", "indianred3"),
                              cols = c("cyan", "magenta", "yellow"),
                              # cols = "Set2",
                              # cols = g_2_6_gex@meta.data$cell.line.ident,
                              group.by = "time.ident",
                              label = F) +
  ggtitle("colored by time")
integrated_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
integrated_dimplot
FeaturePlot(integrated, features = c("BDNF", "VGF"), split.by = "lib.ident")
save("integrated", file = "pilot_integrated_obj.RData")

# assign cell types ####
markers <- c("GAD1", "GAD2", "SLC17A6",
                     "EBF1", # striatal
                     "ADCYAP1", "TCERG1L", "SEMA3E", # subcerebral
                     "CUX1", "BCL11B",  # cortical
                     "SST", # inhibitory
                     "CUX2", "NEFM", # excitatory
                     "VIM", "SOX2", "NES", #NPC
                     "MAP2", "DCX",
                     "SLC17A7", "SERTAD4", "FOXG1", # forebrain
                   "POU3F2", "LHX2" # general cortex
)

StackedVlnPlot(obj = integrated, features = markers) +
  coord_flip()


FeaturePlot(integrated, features = c("GAD1", "GAD2"))
FeaturePlot(integrated, features = c("SLC17A6", "SLC17A7"))
FeaturePlot(integrated, features = c("FOXG1", "NES"), max.cutoff = 10)
FeaturePlot(integrated, features = c("VIM", "SOX2"))

new.cluster.ids <- 
  c("GABA", "glut", "glut", "NPC", "forebrain NPC", 
    "glut", "NPC", "GABA", "unknown", "unknown", 
    "NPC", "unknown", "SEMA3E glut", "unknown")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated$seurat_clusters))
names(new.cluster.ids) <- levels(integrated)
integrated <- RenameIdents(integrated, new.cluster.ids)

integrated$cell.type <- integrated@active.ident
unique(integrated$cell.type)

three_celltype <- subset(integrated, cell.type %in% c("glut", "GABA", "NPC"))
three_celltype <- subset(three_celltype, cell.line != "unmatched")

Idents(three_celltype)
three_celltype$time.ident <- NA
three_celltype$time.ident[three_celltype$orig.ident %in% "SID0"] <- "0hr"
three_celltype$time.ident[three_celltype$orig.ident %in% "SID1"] <- "1hr"
three_celltype$time.ident[three_celltype$orig.ident %in% "SID6"] <- "6hr"
unique(three_celltype$time.ident)

# sum(CD_15$orig.ident == "SID6" & CD_15$cell.type == "NPC")
glut <- subset(RNASeq.integrated, cell.type.idents == "Glut")
GABA <- subset(RNASeq.integrated, cell.type.idents == "GABA")
three_celltype$siwei.cell.type <- "others"
three_celltype$siwei.cell.type[three_celltype@assays$RNA@counts@Dimnames[[2]] %in% 
                                 glut@assays$RNA@counts@Dimnames[[2]]] <- "glut"
three_celltype$siwei.cell.type[three_celltype@assays$RNA@counts@Dimnames[[2]] %in% 
                                 GABA@assays$RNA@counts@Dimnames[[2]]] <- "GABA"
unique(three_celltype$siwei.cell.type)
unique(three_celltype$cell.line)

Idents(three_celltype) <- "siwei.cell.type"

CD_05 <- subset(three_celltype, cell.line == "CD_05")
CD_14 <- subset(three_celltype, cell.line == "CD_14")
CD_15 <- subset(three_celltype, cell.line == "CD_15")
DotPlot(CD_15, assay = "SCT", 
        features = rev(c("FOS", "BDNF", "NPAS4", "SP4", "TCF4")), 
        cols = c("blue", "blue", "blue", "blue", "blue"),
        split.by = "orig.ident") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  coord_flip()

# differential expression ####
j <- 0
# calc hour 1 vs 0
for (i in sort(unique(CD_15$siwei.cell.type))) {
  print(i)
  j = j + 1
  if (j == 1) {
    print("first cell type")
    df_1vs0 <- FindMarkers(object = CD_15,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0,
                           features = c("BDNF", "VGF", "FOS", "NPAS4"),
                           test.use = "MAST")
    df_1vs0$gene_symbol <- rownames(df_1vs0)
    df_1vs0$cell_type <- i
  } else {
    print("others")
    df_to_append <- FindMarkers(object = CD_15,
                                ident.1 = "1hr",
                                group.by = 'time.ident',
                                subset.ident = i,
                                ident.2 = "0hr",
                                min.pct = 0.0,
                                logfc.threshold = 0.0,
                                features = c("BDNF", "VGF", "FOS", "NPAS4"),
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
for (i in sort(unique(CD_15$siwei.cell.type))) {
  print(i)
  j = j + 1
  if (j == 1) {
    print("first cell type")
    df_6vs0 <- FindMarkers(object = CD_15,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0,
                           features = c("BDNF", "VGF", "FOS", "NPAS4"),
                           test.use = "MAST")
    df_6vs0$gene_symbol <- rownames(df_6vs0)
    df_6vs0$cell_type <- i
  } else {
    print("others")
    df_to_append <- FindMarkers(object = CD_15,
                                ident.1 = "6hr",
                                group.by = 'time.ident',
                                subset.ident = i,
                                ident.2 = "0hr",
                                min.pct = 0.0,
                                logfc.threshold = 0.0,
                                features = c("BDNF", "VGF", "FOS", "NPAS4"),
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
df_to_plot$source <- c(rep_len("1v0hr", length.out = 8),
                       rep_len("6v0hr", length.out = 8))


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
