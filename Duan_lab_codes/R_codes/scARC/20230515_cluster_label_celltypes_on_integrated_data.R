# Chuxuan Li 05/15/2023
# After integrating all libraries, cluster and assign cell types

# init ####
library(Seurat)
library(ggplot2)
library(pals)
library(stringr)
library(future)
library(edgeR)
library(harmony)

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

load("integrated_018-030_RNAseq_obj_test_QC.RData")

# QC ####
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, assay = "RNA",
                                                   pattern = "^MT-")
DefaultAssay(integrated) <- "RNA"
Idents(integrated) <- "orig.ident"
VlnPlot(integrated, features = "nFeature_RNA", pt.size = 0) + NoLegend()
VlnPlot(integrated, features = "nCount_RNA", pt.size = 0) + NoLegend()
VlnPlot(integrated, features = "percent.mt", pt.size = 0) + NoLegend()


# clustering ####
unique(integrated$orig.ident)
integrated$time.ident <- paste0(str_sub(integrated$orig.ident, start = -1L), "hr")
unique(integrated$time.ident)

DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, verbose = T, seed.use = 11)
# run Harmony
integrated <- integrated %>% RunHarmony("orig.ident", plot_convergence = TRUE)
DimPlot(object = integrated, reduction = "harmony", group.by = "orig.ident")

# Run the standard workflow for visualization and clustering
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:20, seed.use = 10)
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:20)
# cluster resolution to 0.5 produced 22 clusters 
integrated <- FindClusters(integrated, resolution = 0.5, random.seed = 10)

integrated <- FindClusters(integrated, resolution = 0.2, random.seed = 10)

DimPlot(integrated, label = T, repel = T, group.by = "seurat_clusters") +
  NoLegend() +
  ggtitle("Clustering 018-030 RNAseq data after integration") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
DimPlot(integrated, label = F, group.by = "time.ident") +
  ggtitle("by time point") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
DimPlot(integrated, label = F, group.by = "orig.ident") +
  ggtitle("by library") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10)) +
  NoLegend()

# cell type identification ####
DefaultAssay(integrated) <- "SCT"
FeaturePlot(integrated, features = c("SOX2", "VIM"))
FeaturePlot(integrated, features = c("GAD1", "GAD2"))
FeaturePlot(integrated, features = c("SLC17A6", "SLC17A7"))
FeaturePlot(integrated, features = c("POU5F1", "NANOG"))
FeaturePlot(integrated, features = c("NEFM", "MAP2"))


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
  c("GABA", "NEFM- glut",  "NEFM+ glut", "NEFM- glut", "SEMA3E+ GABA", #0-4
    "NEFM+ glut", "GABA", "GABA", "NEFM+ glut", "SST+ GABA", #5-9
    "glut?", "glut?", "glut?", "glut?", "immature neuron", #10-14
    "NEFM- glut", "immature neuron", "VIM+ cells", "CUX2+ GABA", "GABA", #15-19
    "immature neuron", "unknown", "unknown")

unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated$seurat_clusters))

names(new.cluster.ids) <- levels(integrated)
integrated_labeled <- RenameIdents(integrated, new.cluster.ids)

DimPlot(integrated_labeled, reduction = "umap", label = TRUE, repel = T,
        pt.size = 0.3, cols = c("#B33E52", #GABA
                                "#E6D2B8", #nmglut,
                                "#CCAA7A", #npglut
                                "#CC7A88", #SEMA3E GABA
                                "#E6B8BF", #SST GABA
                                "#ad6e15", #glut?
                                "#0075DC", #immature neuron
                                "#11800d", #VIM+ cells
                                "#f5bfc8", #CUX2 GABA
                                "#6200e3" #unknown
        )) +
  ggtitle("018-030 combined data labeled by cell type") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 10))


integrated_labeled$fine.cell.type <- integrated_labeled@active.ident
integrated_labeled$cell.type <- as.character(integrated_labeled$fine.cell.type)
integrated_labeled$cell.type[integrated_labeled$fine.cell.type %in% 
                                c("SST+ GABA", "GABA", "SEMA3E+ GABA", "CUX2+ GABA")] <- "GABA"
integrated_labeled$cell.type[integrated_labeled$cell.type %in% 
                               c("glut", "unknown", "VIM+ cells", 
                                 "immature neuron")] <- "unidentified"
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
        cols = c("#B33E52", "#E6D2B8", "#CCAA7A", "#ad6e15", "#54990F")) +
  ggtitle("labeled by cell type + respective cell counts") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

save(integrated_labeled, file = "018-030_RNA_integrated_labeled_with_harmony.RData")

# examine MERFISH markers ####
MER_genes <- rev(c(#"CTGF", 
  "PDGFRA", "OPALIN", "SELPLG", "AQP4", "SP8", "KLF5",
                   "LGI2", "LAMP5", "GAD1", "TSHZ2", "CPLX3", "C1QL3", "SYT6", 
                   "SMYD1", "TRABD2A", "COL21A1", "RORB", "CUX2", "SLC17A7", "SLC17A6"))
MER_genes <- as.factor(MER_genes)
Idents(integrated_labeled) <- "seurat_clusters"
StackedVlnPlot(obj = integrated_labeled, features = MER_genes) +
  coord_flip() +
  theme(axis.title.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
