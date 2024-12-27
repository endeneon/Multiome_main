# Siwei 24 Jul 2024
# Extract group 21 from 025 datasets

# init ####
{
  library(Seurat)
  library(Signac)
  
  library(ggplot2)
  library(pals)
  
  library(stringr)
  library(future)
  
  library(dplyr)
  library(harmony)
}

plan("multisession", workers = 8)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

# load data #####
load("~/NVME/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/integrated_labeled.RData")

unique(integrated_labeled$orig.ident)

Idents(integrated_labeled) <- "orig.ident"

subset_group_21 <-
  subset(integrated_labeled,
         subset = (orig.ident %in% c("21-0",
                                     "21-1",
                                     "21-6")))



DimPlot(subset_group_21)
unique(Idents(subset_group_21))

raw_subset_group_21 <-
  CreateSeuratObject(counts = subset_group_21@assays$RNA@counts,
                     meta.data = subset_group_21@meta.data,
                     project = "group_21",
                     min.cells = 0,
                     min.features = 0)

# subset_group_21 <-
#   Seurat::SplitObject(subset_group_21,
#                       split.by = "ident")
raw_subset_group_21 <-
  raw_subset_group_21 %>%
  NormalizeData(verbose = T) #%>%

raw_subset_group_21 <-
  FindVariableFeatures(raw_subset_group_21,
                       verbose = T)



raw_subset_group_21 <-
  raw_subset_group_21 %>%
  ScaleData(verbose = T) %>%
  RunPCA(npcs = 50,
         verbose = T)

raw_subset_group_21 <-
  harmony::RunHarmony(raw_subset_group_21,
                      reduction = "pca",
                      group.by.vars = "orig.ident",
                      ncores = 32,
                      plot_convergence = T,
                      n.seed = 42,
                      reduction.save = "harmony",
                      assay.use = "RNA",
                      verbose = T)
# DefaultAssay(subset_group_21) <- "RNA"

raw_subset_group_21_umap <-
  RunUMAP(raw_subset_group_21,
          dims = 1:50,
          reduction = "harmony")

raw_subset_group_21_umap <-
  raw_subset_group_21_umap %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:50,
                verbose = T) %>%
  FindClusters(resolution = 0.1,
               random.seed = 42,
               verbose = T)

Idents(raw_subset_group_21_umap) <-
  "cell.type.forplot"
Idents(raw_subset_group_21_umap) <-
  "seurat_clusters"
Idents(raw_subset_group_21_umap) <-
  "cell.line.ident"
Idents(raw_subset_group_21_umap) <-
  "time.ident"

DimPlot(raw_subset_group_21_umap,
        label = F,
        repel = T, 
        alpha = 0.5,
        seed = 42) #+
  ggtitle("Resolution = 0.1")
  
FeaturePlot(raw_subset_group_21_umap,
            features = c("GAD1",
                         "SLC17A6",
                         "SLC17A7",
                         "NEFM"),
            alpha = 0.5,
            ncol = 2)
# subset_pca_group21 <-
#   sapply(subset_group_21,
#          FUN = function(x) {
#            x <-
#          })

save(list = c("raw_subset_group_21_umap",
              "subset_group_21"),
     file = "group_21_plot_umap.RData")
