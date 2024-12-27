# Chuxuan Li 04/27/2023
# integrate rabbit cells 

# init ####
library(Seurat)
library(ggplot2)
library(pals)
library(stringr)
library(future)

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

load("./oc2_mm10_analysis/rabbit_only_objlist_from_20230421_analysis_before_normalization.Rdata")

# normalize by individual library ####
dimreduclust <- function(obj){
  obj <- ScaleData(obj)
  obj <- RunPCA(obj,
                verbose = T,
                seed.use = 2023)
  return(obj)
}

processed_lst <- vector("list", length(objlist))
for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  print(paste0("now at ", i))
  # proceed with normalization
  #obj <- PercentageFeatureSet(objlist[[i]], pattern = c("^MT-"), col.name = "percent.mt")
  obj <- FindVariableFeatures(obj, nfeatures = 3000)
  obj <- SCTransform(obj, method = "glmGamPoi", 
                     return.only.var.genes = F, variable.features.n = 8000, 
                     seed.use = 2022, verbose = T)
  processed_lst[[i]] <- dimreduclust(obj)
}

# integration ####
features <- SelectIntegrationFeatures(object.list = processed_lst,
                                      nfeatures = 3000,
                                      fvf.nfeatures = 3000)
processed_lst <- PrepSCTIntegration(processed_lst, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = processed_lst,
                                  anchor.features = features,
                                  reference = c(1, 2, 3),
                                  normalization.method = "SCT",
                                  dims = 1:50)
all.genes <- processed_lst[[1]]@assays$RNA@counts@Dimnames[[1]]

# integrate
integrated <- IntegrateData(anchorset = anchors, features.to.integrate = all.genes,
                            verbose = T)
save(integrated, file = "rabbit_cells_only_integrated.RData")

# QC ####
DefaultAssay(integrated) <- "RNA"
Idents(integrated) <- "orig.ident"
VlnPlot(integrated, features = c("nFeature_RNA", "nCount_RNA"), 
        ncol = 3, fill.by = "feature", pt.size = 0)
VlnPlot(integrated, features = c("nFeature_RNA"), 
        ncol = 3, fill.by = "feature", pt.size = 0) +
  geom_hline(yintercept = c(300, 5000))
VlnPlot(integrated, features = c("nCount_RNA"), 
        ncol = 3, fill.by = "feature", pt.size = 0) +
  geom_hline(yintercept = c(400, 15000))



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
  ggtitle("Clustering 029 RNAseq data with rabbit cells only") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

lib_colors <- DiscretePalette(n = length(unique(Idents(integrated))), "alphabet")
DimPlot(integrated, label = F, cols = lib_colors,
        group.by = "orig.ident") +
  ggtitle("by library") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
DimPlot(integrated, label = F, group.by = "time.ident") +
  ggtitle("by time point") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

save(integrated, file = "../Duan_project_029_RNA/oc2_mm10_analysis/integrated_rabbit_only_object_no_filtering.RData")

# filter 1 ####
integrated$discard.cluster <- "keep"
integrated$discard.cluster[integrated$seurat_clusters %in% c(0, 5)] <- "remove"

integrated_filtered <- subset(integrated, 
                              nFeature_RNA < 5000 & nFeature_RNA > 300 &
                                discard.cluster == "keep")
sum(integrated$nFeature_RNA < 5000 & integrated$nFeature_RNA > 300 &
      integrated$discard.cluster == "keep") 
DefaultAssay(integrated_filtered) <- "integrated"
# Run the standard workflow for visualization and clustering
integrated_filtered <- ScaleData(integrated_filtered)
integrated_filtered <- RunPCA(integrated_filtered, verbose = T, seed.use = 11)
integrated_filtered <- RunUMAP(integrated_filtered, reduction = "pca",
                      dims = 1:30, seed.use = 11)
integrated_filtered <- FindNeighbors(integrated_filtered, reduction = "pca", dims = 1:30)
integrated_filtered <- FindClusters(integrated_filtered, resolution = 0.5, random.seed = 11)
DimPlot(integrated_filtered, label = T, group.by = "seurat_clusters") +
  NoLegend() +
  ggtitle("Clustering 029 RNAseq data with rabbit cells only") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

lib_colors <- DiscretePalette(n = length(unique(Idents(integrated_filtered))), "alphabet")
DimPlot(integrated_filtered, label = F, cols = lib_colors,
        group.by = "orig.ident") +
  ggtitle("by library") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
DimPlot(integrated_filtered, label = F, group.by = "time.ident") +
  ggtitle("by time point") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

VlnPlot(integrated_filtered, features = c("nFeature_RNA"), 
        ncol = 2, fill.by = "feature", pt.size = 0) +
  geom_hline(yintercept = 500)

#save(integrated_filtered, file = "../Duan_project_029_RNA/oc2_mm10_analysis/integrated_rabbit_only_object_filtered1.RData")

# filter 2 ####
integrated$discard.cluster <- "keep"
integrated$discard.cluster[integrated$seurat_clusters %in% c(0, 5)] <- "remove"

integrated_filtered <- subset(integrated, 
                              nFeature_RNA < 5000 & nFeature_RNA > 500 &
                                discard.cluster == "keep")
sum(integrated$nFeature_RNA < 5000 & integrated$nFeature_RNA > 500 &
      integrated$discard.cluster == "keep") #9400

DefaultAssay(integrated_filtered) <- "integrated"
# Run the standard workflow for visualization and clustering
integrated_filtered <- ScaleData(integrated_filtered)
integrated_filtered <- RunPCA(integrated_filtered, verbose = T, seed.use = 11)
integrated_filtered <- RunUMAP(integrated_filtered, reduction = "pca",
                               dims = 1:30, seed.use = 11)
integrated_filtered <- FindNeighbors(integrated_filtered, reduction = "pca", dims = 1:30)
integrated_filtered <- FindClusters(integrated_filtered, resolution = 0.5, random.seed = 11)
DimPlot(integrated_filtered, label = T, group.by = "seurat_clusters") +
  NoLegend() +
  ggtitle("Clustering 029 RNAseq data with rabbit cells only") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

lib_colors <- DiscretePalette(n = length(unique(Idents(integrated_filtered))), "alphabet")
DimPlot(integrated_filtered, label = F, cols = lib_colors,
        group.by = "orig.ident") +
  ggtitle("by library") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
DimPlot(integrated_filtered, label = F, group.by = "time.ident") +
  ggtitle("by time point") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

VlnPlot(integrated_filtered, features = c("nFeature_RNA"), 
        ncol = 2, fill.by = "feature", pt.size = 0) +
  geom_hline(yintercept = 500)

save(integrated_filtered, file = "../Duan_project_029_RNA/oc2_mm10_analysis/integrated_rabbit_only_object_filtered2.RData")
