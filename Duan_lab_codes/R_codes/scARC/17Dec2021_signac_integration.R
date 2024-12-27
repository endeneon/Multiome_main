# Chuxuan Li 12/17/2021
# Use Signac to integrate ATACseq data

# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(readr)
library(ggplot2)
library(future)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(GenomicRanges)
library(AnnotationDbi)
library(patchwork)

library(stringr)
library(uwot)

# read in the fragment files in a loop
path_list <- c("/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_0/outs/atac_fragments.tsv.gz",
               "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_1/outs/atac_fragments.tsv.gz")
               # "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_6/outs/atac_fragments.tsv.gz",
               # "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_0/outs/atac_fragments.tsv.gz",
               # "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_1/outs/atac_fragments.tsv.gz",
               # "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_6/outs/atac_fragments.tsv.gz")

# count fragments per cell, then filter
count_lst <- vector(mode = "list", length = length(path_list))
for (i in 1:length(path_list)){
  count_lst[[i]] <- CountFragments(fragments = path_list[i])
  count_lst[[i]] <- count_lst[[i]][count_lst[[i]]$frequency_count > 2000, "CB"]
}

# create the fragment object
frag_lst <- vector(mode = "list", length = length(path_list))
for (i in 1:length(path_list)){
  frag_lst[[i]] <- CreateFragmentObject(path = path_list[i],
                                        cells = count_lst[[i]])
}

# quantify counts in peaks
peak_lst <- vector(mode = "list", length = length(path_list))
for (i in 1:length(path_list)){
  peak_lst[[i]] <- FeatureMatrix(
    fragments = frag_lst[[i]],
    features = peaks_combined,
    cells = count_lst[[i]]
  )
}

# create chromatin assay object
obj_lst <- vector(mode = "list", length = length(path_list))
for (i in 1:length(obj_lst)){
  obj <- CreateChromatinAssay(
    counts = peak_lst[[i]],
    min.features = 500,
    fragments = frag_lst[[i]]
  )
  obj_lst[[i]] <- CreateSeuratObject(counts = obj, assay = "peaks")
}

for (i in 1:length(obj_lst)){
  obj_lst[[i]] <- subset(obj_lst[[i]], nCount_peaks > 2000 & nCount_peaks < 30000)
}

for (i in 1:length(obj_lst)){
  obj_lst[[i]] <- FindTopFeatures(obj_lst[[i]], min.cutoff = 10)
  obj_lst[[i]] <- RunTFIDF(obj_lst[[i]])
  obj_lst[[i]] <- RunSVD(obj_lst[[i]])
}

lib_idents <- c("g_2_0", "g_2_1", "g_2_6", "g_8_0", "g_8_1", "g_8_6")
for (i in 1:length(obj_lst)){
  obj_lst[[i]]$library <- lib_idents[i]
}


# merge first to get an uncorrected LSI embedding
# failed attempt to do recursion
# merging <- function(l, i){
#   if (i == 1){
#     return(merge(l[[i]], l[[i + 1]]))
#   } else {
#     i = i - 1
#     return(merge(merging(l, i)), l[[i + 1]])
#   }
# }
# merged_signac <- merging(obj_lst, length(obj_lst))
merged_signac <- merge(obj_lst[[1]], obj_lst[[2]])
# for (i in 3:length(obj_lst)){
#   merged_signac <- merge(merged_signac, obj_lst[[i]])
# }

# process the combined dataset
merged_signac <- FindTopFeatures(merged_signac, min.cutoff = 10)
merged_signac <- RunTFIDF(merged_signac)
merged_signac <- RunSVD(merged_signac)
merged_signac <- RunUMAP(merged_signac, reduction = "lsi", dims = 2:30)

DimPlot(merged_signac, group.by = "library")

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = c(obj_lst[[1]], obj_lst[[2]]),
  anchor.features = rownames(obj_lst[[1]]),
  reduction = "rlsi",
  normalization.method = "SCT",
  dims = 2:30
)
merged_signac@assays$peaks@counts@Dimnames[[2]] <-
  str_sub(merged_signac@assays$peaks@counts@Dimnames[[2]], end = -3L)
merged_signac@assays$peaks@data@Dimnames[[2]] <-
  str_sub(merged_signac@assays$peaks@data@Dimnames[[2]], end = -3L)
rownames(merged_signac@reductions$lsi@cell.embeddings) <-
  str_sub(rownames(merged_signac@reductions$lsi@cell.embeddings), end = -3L)

for (i in 1:length(integration.anchors@object.list)){
  rownames(integration.anchors@object.list[[i]]@reductions$lsi@cell.embeddings) <-
    str_sub(rownames(integration.anchors@object.list[[i]]@reductions$lsi@cell.embeddings), end = -3L)
  integration.anchors@object.list[[i]]@assays$peaks@counts@Dimnames[[2]] <-
    str_sub(integration.anchors@object.list[[i]]@assays$peaks@counts@Dimnames[[2]], end = -3L)
  integration.anchors@object.list[[i]]@assays$peaks@data@Dimnames[[2]] <-
    str_sub(integration.anchors@object.list[[i]]@assays$peaks@data@Dimnames[[2]], end = -3L)
}

load("22Dec2021_integration_partial.RData")
# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = merged_signac[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(integrated, group.by = "library")


save.image("17Dec2021_use_signac_integration_after_featurematrix.RData")
