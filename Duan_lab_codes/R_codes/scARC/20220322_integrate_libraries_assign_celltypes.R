# Chuxuan Li 03/22/2022
# From independently normalized 20 libraries, integrate data and assign cell types

# init ####
{
  library(Seurat)
  library(Signac)
  library(ggplot2)
  
  library(future)
  library(RColorBrewer)
}

source("../Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/stacked_vln_copy.R")

plan("multisession", workers = 4)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

# function StackedVlnPlot
source("stacked_vln_copy.R")


load("removed_rat_genes_list.RData")

# integration using reciprocal PCA ####
# find anchors
features <- SelectIntegrationFeatures(object.list = finalobj_lst,
                                      nfeatures = 3000,
                                      fvf.nfeatures = 3000)
finalobj_lst <- PrepSCTIntegration(finalobj_lst, anchor.features = features)
finalobj_lst <- lapply(X = finalobj_lst, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(object.list = finalobj_lst,
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

save("integrated", file = "demuxed_obj_removed_rat_gene_after_integration_nfeat_3000_ftoi_features_updated.RData")

load("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/demuxed_obj_removed_rat_gene_after_integration_nfeat_3000_ftoi_features_updated.RData")
# QC ####
DefaultAssay(integrated) <- "RNA"

integrated[["percent.mt"]] <- 
  PercentageFeatureSet(integrated, 
                       # assay = "RNA",
                       pattern = "^MT-")
DefaultAssay(integrated) <- "RNA"
Idents(integrated) <- "rat.ident"
VlnPlot(integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, fill.by = "feature", pt.size = 0)

sum(integrated$percent.mt > 10) #20731
sum(integrated$percent.mt > 15) #9937
sum(integrated$percent.mt > 20) #750

obj <- 
  integrated[rownames(integrated) %in% 
               rownames(integrated)[str_detect(rownames(integrated), 
                                               "^MT-")],]
rowSums(obj)

sum(integrated$nFeature_RNA > 8000) #1130
sum(integrated$nFeature_RNA < 400) #194
sum(integrated$nCount_RNA < 500) #873
sum(integrated$nCount_RNA > 40000) #765

sum(integrated$percent.mt > 15 | integrated$nFeature_RNA > 8000 | integrated$nCount_RNA > 40000)

# clustering ####
unique(integrated$orig.ident)
unique(integrated$time.ident)

DefaultAssay(integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated)
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
save(integrated,
     file = "Siwei_replot_FigS1_20Aug2024.RData")
# q(save = "no")

load("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/Siwei_replot_FigS1_20Aug2024.RData")

DimPlot(integrated,
        label = T,
        group.by = "seurat_clusters") +
  NoLegend() +
  ggtitle("clustering group 5-31-35-51-63 RNAseq data") +
  theme(text = element_text(size = 12))

DimPlot(integrated,
        label = F,
        group.by = "orig.ident") +
  ggtitle("by library") +
  theme(text = element_text(size = 12))

DimPlot(integrated,
        label = F,
        group.by = "time.ident") +
  ggtitle("by time point") +
  theme(text = element_text(size = 12))

# cell type identification ####
DefaultAssay(integrated) <- "SCT"
FeaturePlot(integrated,
            features = c("SOX2", "VIM"))
FeaturePlot(integrated,
            features = c("GAD1", "SLC17A6"))
FeaturePlot(integrated,
            features = c("POU5F1", "NANOG"))
FeaturePlot(integrated,
            features = c("MAP2"))
FeaturePlot(integrated,
            features = c("NEFM", "CUX2"))


trimmed_markers <- c("GAD1", "GAD2", "SLC17A6",
                     "EBF1", # striatal
                     "SEMA3E", # subcerebral
                     "BCL11B",  # cortical
                     "SST", # inhibitory
                     "SATB2",  "NEFM", # excitatory
                     "VIM",
                     "SOX2",  #NPC
                     "SLC17A7", "SERTAD4", "FOXG1",  # forebrain
                     "POU3F2", "LHX2", # general cortex
                     "ADCYAP1", "CUX1", "CUX2", 
                     "MAP2", "DCX")

DefaultAssay(integrated) <- "SCT"
StackedVlnPlot(obj = integrated, 
               features = trimmed_markers) +
  coord_flip()

Idents(integrated_labeled) <- "seurat_clusters"

# assign cell types
new.cluster.ids <-
  c("NEFM_pos_glut", "GABA", "GABA", "NEFM_neg_glut", "NEFM_neg_glut",
    "GABA", "SEMA3E_pos_GABA", "NEFM_pos_glut", "NEFM_neg_glut", "NEFM_pos_glut",
    "unknown", "NEFM_pos_glut", "NEFM_pos_glut", "unknown", "SST_pos_GABA",
    "unknown", "unknown", "SEMA3E_pos_GABA", "unknown", "unknown")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated$seurat_clusters))


names(new.cluster.ids) <- levels(integrated)
integrated_labeled <- RenameIdents(integrated, new.cluster.ids)
DimPlot(integrated_labeled,
        reduction = "umap",
        label = F,
        repel = T,
        cols = brewer.pal(n = 8,
                          name = "Dark2")) +
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 11))

integrated_labeled$cell.type <- integrated_labeled@active.ident
save(integrated_labeled, file = "labeled_nfeat3000.RData")
