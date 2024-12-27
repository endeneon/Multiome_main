# Chuxuan Li 03/17/2023
# After integrating 018-029 RNA data, cluster and assign cell types

# init ####
{
  library(Seurat)
  library(signac)
  library(ggplot2)
  library(pals)
  library(stringr)
  library(future)
  library(edgeR)
  library(harmony)
}

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

load("integrated_obj.RData")

# QC ####
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, assay = "RNA",
                                                   pattern = "^MT-")
DefaultAssay(integrated) <- "RNA"
Idents(integrated) <- "orig.ident"
VlnPlot(integrated, features = "nFeature_RNA", pt.size = 0) + NoLegend()
VlnPlot(integrated, features = "nCount_RNA", pt.size = 0) + NoLegend()
VlnPlot(integrated, features = "percent.mt", pt.size = 0) + NoLegend()

sum(integrated$percent.mt > 15) #18202
sum(integrated$percent.mt > 20) #4314

# clustering ####
unique(integrated$orig.ident)
integrated$time.ident <- paste0(str_sub(integrated$orig.ident, start = -1L), "hr")
unique(integrated$time.ident)

DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, verbose = T, seed.use = 11)

# run UMAP w/o harmony



# run Harmony
integrated <- integrated %>% RunHarmony("orig.ident", plot_convergence = TRUE)
DimPlot(object = integrated, reduction = "harmony", group.by = "orig.ident")

# Run the standard workflow for visualization and clustering
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:20, seed.use = 10)
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.5, random.seed = 10)
DimPlot(integrated, label = T, group.by = "seurat_clusters") +
  NoLegend() +
  ggtitle("Clustering 018-029 RNAseq data after integration") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

DimPlot(integrated, label = F, group.by = "time.ident") +
  ggtitle("by time point") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
DimPlot(integrated, label = F, group.by = "orig.ident") +
  ggtitle("by library") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

# cell type identification ####
DefaultAssay(integrated) <- "SCT"
FeaturePlot(integrated, features = c("SOX2", "VIM"))
FeaturePlot(integrated, features = c("GAD1", "GAD2"))
FeaturePlot(integrated, features = c("SLC17A6", "SLC17A7"))
FeaturePlot(integrated, features = c("POU5F1", "NANOG"))
FeaturePlot(integrated, features = c("MAP2"))
FeaturePlot(integrated, features = c("NEFM", "CUX2"))


trimmed_markers <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
                     "EBF1", # striatal
                     "SEMA3E", # subcerebral
                     "BCL11B",  # cortical
                     "SST", # inhibitory
                     "SATB2",  "NEFM", # excitatory
                     "VIM", "SOX2",  #NPC
                     "SERTAD4", "FOXG1",  # forebrain
                     "POU3F2", "LHX2", # general cortex
                     "ADCYAP1", "CUX1", "CUX2", 
                     "MAP2", "DCX") # pan-neuron

StackedVlnPlot(obj = integrated, features = trimmed_markers) +
  coord_flip()

Idents(integrated) <- "seurat_clusters"

# assign cell types
new.cluster.ids <-
  c("NEFM- glut", "NEFM+ glut",  "GABA", "GABA", "GABA", #0-4
    "?glut", "NEFM+ glut", "?glut", "NPC", "?glut", #5-9
    "unknown", "GABA", "GABA", "GABA" #10-13
    )
unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated$seurat_clusters))

names(new.cluster.ids) <- levels(integrated)
integrated_labeled <- RenameIdents(integrated, new.cluster.ids)

DimPlot(integrated_labeled, reduction = "umap", label = TRUE, repel = T,
        pt.size = 0.3, cols = c("#E6D2B8", #nmglut
                                "#CCAA7A", #npglut
                                "#B33E52", #GABA
                                "#f29116", #?glut
                                "#347545", #NPC
                                "#4293db"#unknown
                              )) +
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 10))


integrated_labeled$fine.cell.type <- integrated_labeled@active.ident
integrated_labeled$cell.type <- as.character(integrated_labeled$fine.cell.type)
# integrated_labeled$cell.type[integrated_labeled$fine.cell.type %in% 
#                                c("SST+ GABA", "GABA")] <- "GABA"
integrated_labeled$cell.type[integrated_labeled$cell.type %in% 
                               c("?glut", "unknown", "stem cell-like",
                                 "NPC", "NPC?")] <- "unidentified"
integrated_labeled$cell.type.forplot <- integrated_labeled$cell.type
integrated_labeled$cell.type[integrated_labeled$fine.cell.type == "NEFM- glut"] <- "nmglut"
integrated_labeled$cell.type[integrated_labeled$fine.cell.type == "NEFM+ glut"] <- "npglut"
unique(integrated_labeled$cell.type)
unique(integrated_labeled$cell.type.forplot)
integrated_labeled$cell.type.counts <- integrated_labeled$cell.type.forplot
types <- unique(integrated_labeled$cell.type.forplot)
for (i in 1:length(types)) {
  count <- sum(integrated_labeled$cell.type.forplot == types[i])
  print(count)
  integrated_labeled$cell.type.counts[integrated_labeled$cell.type.counts == types[i]] <-
    paste0(types[i], "\n", count)
}
unique(integrated_labeled$cell.type.counts)
DimPlot(integrated_labeled, reduction = "umap", group.by = "cell.type.counts", 
        label = TRUE, repel = F, pt.size = 0.3, 
        cols = c("#B33E52", "#E6D2B8", "#CCAA7A", "#54990F")) +
  ggtitle("labeled by cell type + respective cell counts") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

save(integrated_labeled, file = "029_RNA_integrated_labeled_with_harmony.RData")

Idents(integrated_labeled) <- "orig.ident"
integrated_labeled$orig.ident <-
  factor(integrated_labeled$orig.ident,
         levels = sort(unique(integrated_labeled$orig.ident)))


DimPlot(integrated_labeled, 
        reduction = "umap", 
        alpha = 0.8,
        # order = sort(unique(integrated_labeled$orig.ident)),
        pt.size = 0.8,
        label = F, 
        shuffle = T,
        repel = T) +
  # ggtitle("labeled by cell type") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 10)) +
  NoLegend()

sort(unique(integrated_labeled$orig.ident))


DefaultAssay(integrated_labeled) <- "ATAC"

Coverage