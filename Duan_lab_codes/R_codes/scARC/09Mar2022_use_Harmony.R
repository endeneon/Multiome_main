# Chuxuan Li 03/09/2022
# Use Harmony on the integrated object to test normalization

# init ####
library(Seurat)
library(harmony)

library(ggplot2)
library(RCsolorBrewer)
library(cowplot)

library(stringr)
library(future)

# remove lib.ident as variable ####
RNAseq_integrated_labeled <- RNAseq_integrated_labeled %>% 
  RunHarmony("lib.ident", plot_convergence = TRUE)
DimPlot(object = RNAseq_integrated_labeled, reduction = "harmony", group.by = "lib.ident")

RNAseq_integrated_labeled <- RunUMAP(RNAseq_integrated_labeled, 
                        reduction = "harmony", 
                        dims = 1:20,
                        seed.use = 10)
RNAseq_integrated_labeled <- FindNeighbors(RNAseq_integrated_labeled, 
                            reduction = "harmony", 
                            dims = 1:20)
RNAseq_integrated_labeled <- FindClusters(RNAseq_integrated_labeled, 
                           resolution = 0.25,
                           random.seed = 10)
DimPlot(object = RNAseq_integrated_labeled, reduction = "umap", label = T)
by_time_dimplot <- DimPlot(RNAseq_integrated_labeled,
                              cols = c("blue", "tomato", "gold"),
                              group.by = "time.ident",
                              label = F) +
  ggtitle("colored by time")
by_time_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
by_time_dimplot

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

StackedVlnPlot(obj = RNAseq_integrated_labeled, features = markers) +
  coord_flip()


FeaturePlot(integrated, features = c("GAD1", "GAD2"))
FeaturePlot(integrated, features = c("SLC17A6", "SLC17A7"))
FeaturePlot(integrated, features = c("FOXG1", "NES"), max.cutoff = 10)
FeaturePlot(integrated, features = c("VIM", "SOX2"))

new.cluster.ids <- 
  c("NEFM_neg_glut", "NEFM_pos_glut", "GABA", "SEMA3E_pos_glut", "NEFM_pos_glut", 
    "NEFM_neg_glut", "GABA", "GABA", "NEFM_pos_glut", "unknown", 
    "NPC", "FOXG1_neuron", "SST+ glut", "unknown", "NEFM_neg_glut", "unknown", "unknown")

length(new.cluster.ids)
length(unique(RNAseq_integrated_labeled$seurat_clusters))
names(new.cluster.ids) <- levels(RNAseq_integrated_labeled)
labeled <- RenameIdents(RNAseq_integrated_labeled, new.cluster.ids)

labeled$cell.type <- labeled@active.ident
unique(labeled$cell.type)

four_celltype <- subset(labeled, cell.type %in% c("NEFM_neg_glut", "NEFM_pos_glut", "GABA", "NPC"))
DefaultAssay(four_celltype) <- "RNA"
four_celltype <- ScaleData(four_celltype)
DotPlot(four_celltype, 
        features = rev(c("FOS", "BDNF", "NPAS4", "VGF")), 
        cols = c("blue", "blue", "blue", "blue", "blue"),
        split.by = "time.ident") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  coord_flip()



# remove group.ident as variable ####
RNAseq_integrated_labeled <- RNAseq_integrated_labeled %>% 
  RunHarmony("group.ident", plot_convergence = TRUE)
DimPlot(object = RNAseq_integrated_labeled, reduction = "harmony", group.by = "group.ident")
