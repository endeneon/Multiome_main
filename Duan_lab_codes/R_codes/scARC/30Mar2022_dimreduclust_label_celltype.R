# Chuxuan Li 03/30/2022
# After integration of 5-line and 18-line data, perform downstream analyses
#including dimensional reduction, clustering, labeling cell types

# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(SeuratWrappers)
library(harmony)

library(stringr)
library(ggplot2)

library(future)
plan("multisession", workers = 1)

# load data 
load("./5line_18line_combined_after_integration_nfeat_5000.RData")

# dimensional reduction and clustering ####
DefaultAssay(integrated)
# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated,
                     verbose = T,
                     seed.use = 11)
integrated <- RunHarmony(integrated, group.by.vars = "orig.ident")
integrated <- RunUMAP(integrated,
                      reduction = "harmony",
                      dims = 1:30,
                      seed.use = 11)
integrated <- FindNeighbors(integrated,
                            reduction = "harmony",
                            dims = 1:30)
integrated <- FindClusters(integrated,
                           resolution = 0.5,
                           random.seed = 11)

unique(integrated$orig.ident)
integrated$orig.ident <- str_replace_all(integrated$orig.ident, "_", "-")
integrated$time.ident <- paste0(str_extract(integrated$orig.ident, "[0|1|6]$"),
                                "hr")
unique(integrated$time.ident)
integrated$batch.ident <- "18line"
integrated$batch.ident[integrated$orig.ident %in% c("2-0",  "2-1",  "2-6",
                                                    "8-0" , "8-1",  "8-6")] <- "5line"
unique(integrated$batch.ident)

DimPlot(integrated,
        label = T,
        group.by = "seurat_clusters") +
  NoLegend() +
  ggtitle("clustering of 5-line and 18-line combined RNAseq data") +
  theme(text = element_text(size = 12))

p <- DimPlot(integrated,
        label = F,
        group.by = "orig.ident") +
  ggtitle("by library") +
  theme(text = element_text(size = 12))
p[[1]]$layers[[1]]$aes_params$alpha = .2
p

p <- DimPlot(integrated,
        label = F,
        group.by = "time.ident") +
  ggtitle("by time point") +
  theme(text = element_text(size = 12))
p[[1]]$layers[[1]]$aes_params$alpha = .2
p


# cell type identification ####
FeaturePlot(integrated,
            features = c("SOX2", "VIM"), min.cutoff = 0)
FeaturePlot(integrated,
            features = c("GAD1", "SLC17A6"), min.cutoff = 0)
FeaturePlot(integrated,
            features = c("POU5F1", "NANOG"))
FeaturePlot(integrated,
            features = c("MAP2"))
FeaturePlot(integrated,
            features = c("NEFM", "CUX2"), min.cutoff = 0)

DefaultAssay(integrated) <- "SCT"
trimmed_markers <- c("GAD1", "GAD2", "SLC17A6",
                     "EBF1", # striatal
                     "SEMA3E", # subcerebral
                     "BCL11B",  # cortical
                     "SST", # inhibitory
                     "SATB2",  "NEFM", # excitatory
                     "VIM", "SOX2",  #NPC
                     "SLC17A7", "SERTAD4", "FOXG1",  # forebrain
                     "POU3F2", "LHX2", # general cortex
                     "ADCYAP1", "CUX1", "CUX2", 
                     "MAP2", "DCX")

StackedVlnPlot(obj = integrated, features = trimmed_markers) +
  coord_flip()


# assign cell types
new.cluster.ids <-
  c("NEFM_pos_glut", "GABA", "GABA", "NEFM_neg_glut", "NEFM_neg_glut",
    "NEFM_pos_glut", "NEFM_pos_glut", "GABA", "NEFM_neg_glut", "SEMA3E_pos_glut",
    "NEFM_neg_glut", "SST_pos_GABA", "unknown", "unknown", "SEMA3E_pos_glut",
    "SST_pos_GABA", "NEFM_pos_glut", "immature_neuron", "unknown", "GABA", 
    "immature_neuron", "immature_neuron", "GABA", "NEFM_neg_glut", "GABA",
    "NEFM_neg_glut")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated$seurat_clusters))
names(new.cluster.ids) <- levels(integrated)
integrated_labeled <- RenameIdents(integrated, new.cluster.ids)
integrated_labeled$cell.type <- integrated_labeled@active.ident
integrated_labeled$cell.type.to.plot <- as.vector(integrated_labeled@active.ident)

for (i in unique(integrated_labeled$cell.type.to.plot)){
  print(i)
  integrated_labeled$cell.type.to.plot[integrated_labeled$cell.type.to.plot == i] <-
    str_replace_all(str_replace_all(i, "_pos_", "+ "), "_neg_", "- ")
}

unique(integrated_labeled$cell.type.to.plot)

DimPlot(integrated_labeled,
        reduction = "umap",
        label = TRUE,
        repel = F,
        pt.size = 0.3, 
        group.by = "cell.type.to.plot") +
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 11))

save(integrated_labeled, file = "5line_18line_combined_labeled_nfeat5000_obj.RData")


p <- DimPlot(integrated_labeled,
             label = F,
             group.by = "batch.ident") +
  ggtitle("by batch (5 line vs 18 line)") +
  theme(text = element_text(size = 12))
p[[1]]$layers[[1]]$aes_params$alpha = .2
p
