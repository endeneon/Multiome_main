# Chuxuan Li 11/10/2022
# Normalize -> integrate -> cluster -> label 025 data

# init ####
library(Seurat)
library(ggplot2)
library(pals)
library(stringr)
library(future)

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/GRCh38_mapped_after_QC_list.RData")

# normalize ####
transformed_lst <- vector(mode = "list", length(QCed_lst))

for (i in 1:length(QCed_lst)){
  print(paste0("now at ", i))
  # proceed with normalization
  obj <- PercentageFeatureSet(QCed_lst[[i]], pattern = c("^MT-"), col.name = "percent.mt")
  
  obj <- FindVariableFeatures(obj, nfeatures = 8000)
  transformed_lst[[i]] <- SCTransform(obj,
                                      vars.to.regress = "percent.mt",
                                      method = "glmGamPoi", 
                                      return.only.var.genes = F,
                                      variable.features.n = 8000,
                                      seed.use = 42,
                                      verbose = T)
}

# integration using reciprocal PCA ####
# find anchors
features <- SelectIntegrationFeatures(object.list = transformed_lst,
                                      nfeatures = 3000,
                                      fvf.nfeatures = 3000)
transformed_lst <- PrepSCTIntegration(transformed_lst, anchor.features = features)
transformed_lst <- lapply(X = transformed_lst, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(object.list = transformed_lst,
                                  anchor.features = features,
                                  reference = c(1, 2, 3),
                                  reduction = "rpca",
                                  normalization.method = "SCT",
                                  scale = F, 
                                  dims = 1:50)
#Found 9036 anchors
#all.genes <- transformed_lst[[1]]@assays$RNA@counts@Dimnames[[1]]

# integrate
integrated <- IntegrateData(anchorset = anchors,
                            verbose = T)
save(transformed_lst, file = "transformed_and_prepared_for_integration_list.RData")
save(integrated, file = "integrated_obj_nfeature_8000.RData")

# QC ####
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, assay = "RNA",
                                                   pattern = "^MT-")
DefaultAssay(integrated) <- "RNA"
Idents(integrated) <- "orig.ident"
VlnPlot(integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, fill.by = "feature", pt.size = 0)


# clustering ####
unique(integrated$orig.ident)
unique(integrated$time.ident)

DefaultAssay(integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, verbose = T, seed.use = 11)
integrated <- RunUMAP(integrated, reduction = "pca",
                      dims = 1:30, seed.use = 11)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5, random.seed = 11)
DimPlot(integrated, label = T, group.by = "seurat_clusters") +
  NoLegend() +
  ggtitle("Clustering 025 RNAseq data") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

lib_colors <- DiscretePalette(n = length(unique(Idents(integrated))), "alphabet")
DimPlot(integrated, label = F, cols = lib_colors,
        group.by = "orig.ident") +
  ggtitle("by library") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
lib_colors[4:length(lib_colors)] <- "transparent"
lib_colors[1:2] <- "transparent"
DimPlot(integrated, label = F, pt.size = 0.6, 
        group.by = "orig.ident", cols = lib_colors) +
  ggtitle("library 9-6 distribution") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

DimPlot(integrated, label = F, group.by = "time.ident") +
  ggtitle("by time point") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

# cell type identification ####
DefaultAssay(integrated) <- "SCT"
FeaturePlot(integrated, features = c("SOX2", "VIM"))
FeaturePlot(integrated, features = c("GAD1", "GAD2", "SLC17A6", "SLC17A7"), ncol = 2)
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
                     "SLC17A7", "SERTAD4", "FOXG1",  # forebrain
                     "POU3F2", "LHX2", # general cortex
                     "ADCYAP1", "CUX1", "CUX2", 
                     "MAP2", "DCX")

StackedVlnPlot(obj = integrated, features = trimmed_markers) +
  coord_flip()

Idents(integrated) <- "seurat_clusters"

# assign cell types
new.cluster.ids <-
  c("NEFM+ glut", "GABA", "NEFM- glut", "NEFM- glut", "GABA",
    "NEFM+ glut", "SEMA3E+ GABA", "SST+ GABA", "GABA", "NEFM- glut",
    "GABA", "NEFM+ glut", "immature neuron", "unknown", "NEFM+ glut", 
    "NEFM- glut", "GABA", "NEFM- glut", "unknown neuron", "immature neuron", 
    "GABA", "unknown")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated$seurat_clusters))

names(new.cluster.ids) <- levels(integrated)
integrated_labeled <- RenameIdents(integrated, new.cluster.ids)

DimPlot(integrated_labeled, reduction = "umap", label = TRUE, repel = F,
        pt.size = 0.3, cols = c("#CCAA7A", #npglut
                                "#B33E52", #GABA
                                "#E6D2B8", #nmglut
                                "#CC7A88", #SEMA3E GABA
                                "#E6B8BF", #SST GABA
                                "#0075DC", #immature neuron
                                "#993F00", #unknown
                                "#4C005C" #unknown neuron
                                )) +
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

integrated_labeled$fine.cell.type <- integrated_labeled@active.ident
integrated_labeled$cell.type <- as.character(integrated_labeled$fine.cell.type)
integrated_labeled$cell.type[integrated_labeled$fine.cell.type %in% 
                               c("SEMA3E+ GABA", "SST+ GABA", "GABA")] <- "GABA"
integrated_labeled$cell.type[integrated_labeled$cell.type %in% 
                               c("unknown neuron", "immature neuron", "unknown")] <- "unidentified"
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
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

save(integrated_labeled, file = "integrated_labeled.RData")

# check number of cells in group 09 ####
sum(integrated_labeled$orig.ident %in% c("09-0", "09-1", "09-6"))
sum(integrated_labeled$orig.ident %in% c("09-0", "09-1", "09-6") & 
      integrated_labeled$seurat_clusters %in% c("12", "13"))
