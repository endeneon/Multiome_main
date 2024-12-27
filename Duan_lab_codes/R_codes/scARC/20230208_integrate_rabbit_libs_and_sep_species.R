# Chuxuan Li 02/08/2023
# Integrate all libraries with rabbit together and cluster, identify rabbit and
#mouse cells again

# init ####
library(Seurat)
library(glmGamPoi)
library(sctransform)

library(ggplot2)
library(RColorBrewer)
library(viridis)
library(graphics)
library(ggrepel)

library(dplyr)
library(readr)
library(stringr)
library(readxl)

library(future)
# set threads and parallelization
plan("multisession", workers = 1)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

setwd("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/oc2_mm10_analysis/")
load("labeled_by_species_obj_list.RData")

# integrate objects ####
features <- SelectIntegrationFeatures(object.list = scaled_lst,
                                      nfeatures = 3000,
                                      fvf.nfeatures = 3000)
scaled_lst <- PrepSCTIntegration(scaled_lst, anchor.features = features)
scaled_lst <- lapply(X = scaled_lst, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(object.list = scaled_lst,
                                  anchor.features = features,
                                  reference = c(1, 2, 3),
                                  reduction = "rpca",
                                  normalization.method = "SCT",
                                  scale = F, 
                                  dims = 1:50)

mouse_rabbit_integrated <- IntegrateData(anchorset = anchors,
                            verbose = T)

# clustering ####
unique(mouse_rabbit_integrated$orig.ident)
unique(mouse_rabbit_integrated$time.ident)

DefaultAssay(mouse_rabbit_integrated) <- "integrated"

mouse_rabbit_integrated <- ScaleData(mouse_rabbit_integrated)
mouse_rabbit_integrated <- RunPCA(mouse_rabbit_integrated, verbose = T, seed.use = 11)
mouse_rabbit_integrated <- RunUMAP(mouse_rabbit_integrated, reduction = "pca",
                      dims = 1:30, seed.use = 11)
mouse_rabbit_integrated <- FindNeighbors(mouse_rabbit_integrated, reduction = "pca", dims = 1:30)
mouse_rabbit_integrated <- FindClusters(mouse_rabbit_integrated, resolution = 0.5, random.seed = 11)
DimPlot(mouse_rabbit_integrated, label = T, group.by = "seurat_clusters") +
  NoLegend() +
  ggtitle("Clustering mouse and rabbit cells after integration") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

lib_colors <- DiscretePalette(n = length(unique(Idents(mouse_rabbit_integrated))), "alphabet2")
dimplot <- DimPlot(mouse_rabbit_integrated, label = F, cols = lib_colors,
        group.by = "orig.ident") +
  ggtitle("by library") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
dimplot

dimplot <- DimPlot(mouse_rabbit_integrated, label = F, group.by = "time.ident") +
  ggtitle("by time point") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
dimplot

# expression of markers
DefaultAssay(mouse_rabbit_integrated) <- "SCT"
FeaturePlot(mouse_rabbit_integrated, features = c("Gfap", "S100b",
                                                  "ENSOCUG00000002862", "ENSOCUG00000013326"))

# assign cell types ####
new.cluster.ids <-
  c("rabbit", "mouse", "rabbit", "rabbit", "mouse",
    "rabbit", "mouse", "mouse", "mouse", "unknown",
    "unknown", "unknown", "unknown", "unknown", "unknown",
    "unknown", "unknown", "unknown")
unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(mouse_rabbit_integrated$seurat_clusters))
names(new.cluster.ids) <- levels(mouse_rabbit_integrated)
mouse_rabbit_integrated <- RenameIdents(mouse_rabbit_integrated, new.cluster.ids)
mouse_rabbit_integrated$cell.type <- mouse_rabbit_integrated@active.ident
DimPlot(mouse_rabbit_integrated, group.by = "cell.type") +
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

save(mouse_rabbit_integrated, file = "mouse_rabbit_integrated_obj.RData")


# count number of cells in each species ####
load("./cell_count_by_species_separated_by_lib_analysis.RData")
load("../029_RNA_integrated_labeled.RData")
libs <- sort(unique(mouse_rabbit_integrated$orig.ident))
integrated_labeled <- subset(integrated_labeled, orig.ident %in% libs)
load("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/GRCh38_mapped_raw_list.RData")
j = 0
for (i in 1:length(objlist)) {
  lib <- as.character(unique(objlist[[i]]$orig.ident))
  if (lib %in% libs) {
    j = j + 1
    print(i)
    print(j)
    print(lib)
    subobj_human <- subset(integrated_labeled, orig.ident == lib)
    subobj_nonhuman <- subset(mouse_rabbit_integrated, orig.ident == lib)
    human_bc <- str_extract(subobj_human@assays$RNA@counts@Dimnames[[2]], "[A-Z]+-1")
    rabbit_bc <- str_extract(subobj_nonhuman@assays$RNA@counts@Dimnames[[2]], "[A-Z]+-1")[subobj_nonhuman$cell.type == "rabbit"]
    mouse_bc <- str_extract(subobj_nonhuman@assays$RNA@counts@Dimnames[[2]], "[A-Z]+-1")[subobj_nonhuman$cell.type == "mouse"]
    objlist[[i]]$species.ident <- "unknown"
    objlist[[i]]$species.ident[objlist[[i]]@assays$RNA@counts@Dimnames[[2]] %in% human_bc] <- "human"
    objlist[[i]]$species.ident[objlist[[i]]@assays$RNA@counts@Dimnames[[2]] %in% rabbit_bc] <- "rabbit"
    objlist[[i]]$species.ident[objlist[[i]]@assays$RNA@counts@Dimnames[[2]] %in% mouse_bc] <- "mouse"
    objlist[[i]]$human.ident <- "nonhuman"
    objlist[[i]]$human.ident[objlist[[i]]@assays$RNA@counts@Dimnames[[2]] %in% human_bc] <- "human"
    cell_count_by_species[j, 1] <- length(objlist[[i]]$orig.ident)
    cell_count_by_species[j, 2] <- length(subobj_human$cell.line.ident)
    cell_count_by_species[j, 3] <- sum(objlist[[i]]$human.ident == "nonhuman")
    cell_count_by_species[j, 4] <- sum(subobj_nonhuman$cell.type == "rabbit")
    cell_count_by_species[j, 5] <- sum(subobj_nonhuman$cell.type == "mouse")
    cell_count_by_species[j, 6] <- sum(objlist[[i]]$species.ident == "unknown")
  }
}
write.table(cell_count_by_species, file = "cell_count_by_species_final.csv", 
            quote = F, sep = ",", row.names = T, col.names = T)
