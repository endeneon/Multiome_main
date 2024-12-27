# Chuxuan Li 01/21/2022
# Start with individual libraries (g20, g21, ...), first normalize them individually,
#then integrate the ones from the same time point using Signac/Seurat IntegrateData()

# init ####
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)

library(stringr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)

library(future)

set.seed(1021)
# set threads and parallelization
plan("multisession", workers = 1)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)



# load data #### 
setwd("~/NVME/scARC_Duan_018/GRCh38_mapped_only")

# Read files, separate the ATACseq data from the .h5 matrix
h5list <- list.files(path = ".", 
                     pattern = "filtered_feature_bc_matrix.h5",
                     recursive = T,
                     include.dirs = T)
gzlist <- list.files(path = ".", 
                     pattern = "atac_fragments.tsv.gz$",
                     recursive = T,
                     include.dirs = T)
objlist <- vector(mode = "list", length = length(h5list))



for (i in 1:6){
  h5file <- Read10X_h5(filename = h5list[i])
  chrom_assay <- CreateChromatinAssay(
    counts = h5file$Peaks,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = gzlist[i],
    min.cells = 10,
    min.features = 200
  )
  
  objlist[[i]] <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    project = str_extract(string = h5list[i], 
                          pattern = "libraries_[0-9]_[0-6]"))
  
  
}

# apply the 10x aggr raw peak set to the libraries ####
newpeakslist <- vector(mode = "list", length = length(objlist))
setwd("~/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs")
bedfile <- read.table("atac_peaks_cleaned.bed", sep = "\t")
peaks_aggr <- GRanges(seqnames = bedfile$V1,
                      ranges = IRanges(start = bedfile$V2,
                                       end = bedfile$V3))

setwd("~/NVME/scARC_Duan_018/GRCh38_mapped_only")
for (i in 1:length(objlist)){
  # apply another set of peaks
  macs2_counts_complete <- FeatureMatrix(
    fragments = Fragments(objlist[[i]]),
    features = peaks_aggr,
    cells = colnames(objlist[[i]])
  )
  newpeakslist[[i]] <- objlist[[i]]
  newpeakslist[[i]][["macs2"]] <- CreateChromatinAssay(
    counts = macs2_counts_complete,
    fragments = gzlist[[i]],
    annotation = annot_aggr_signac_ucsc
  )
}

save(newpeakslist, file = "after_assigning_10xaggr_raw_peakset_to_6lib.RData")
# process 0hr libraries first to test ####
hist(newpeakslist[[1]]$nCount_peaks, breaks = 1000)
hist(newpeakslist[[4]]$nCount_peaks, breaks = 1000)

# look at the objects first, before subsetting
test <- subset(newpeakslist[[1]], nCount_peaks > 2000 & nCount_peaks < 30000)
test2 <- subset(newpeakslist[[4]], nCount_peaks > 2000 & nCount_peaks < 30000)

DefaultAssay(test) <- "macs2"
DefaultAssay(test2) <- "macs2"
hist(test$nCount_peaks, breaks = 1000)

# load annotation file because annotations cannot be directly accessed using GetGRangesFromEnsDb
# test with time=0 first
Annotation(test) <- annot_aggr_signac_ucsc
Annotation(test2) <- annot_aggr_signac_ucsc

DefaultAssay(test)
DefaultAssay(test2)

test <- TSSEnrichment(object = test, fast = FALSE)
test2 <- TSSEnrichment(object = test2, fast = FALSE)
test <- NucleosomeSignal(object = test)
test <- subset(test, nCount_peaks > 2000 & 
                            nCount_peaks < 30000 &
                            nucleosome_signal < 4 &
                            TSS.enrichment > 2)
test2 <- NucleosomeSignal(object = test2)
test2 <- subset(test2, nCount_peaks > 2000 & 
                 nCount_peaks < 30000 &
                 nucleosome_signal < 4 &
                 TSS.enrichment > 2)

test <- FindTopFeatures(test, min.cutoff = 10)
test <- RunTFIDF(test)
test <- RunSVD(test)

test2 <- FindTopFeatures(test2, min.cutoff = 10)
test2 <- RunTFIDF(test2)
test2 <- RunSVD(test2)

combined <- merge(test, test2)
combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunTFIDF(combined)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:30)


# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(test, test2),
  anchor.features = rownames(test),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
DimPlot(integrated, group.by = "orig.ident")



# process all objects ####
subobjlist <- vector(mode = "list", length = length(newpeakslist))

for (i in 1:length(newpeakslist)){
  DefaultAssay(newpeakslist[[i]]) <- "macs2"
  Annotation(newpeakslist[[i]]) <- annot_aggr_signac_ucsc
  newpeakslist[[i]] <- TSSEnrichment(object = newpeakslist[[i]], fast = FALSE)
  newpeakslist[[i]] <- NucleosomeSignal(object = newpeakslist[[i]])
  subobjlist[[i]] <- subset(newpeakslist[[i]], nCount_peaks > 2000 & 
                              nCount_peaks < 30000 &
                              nucleosome_signal < 4 &
                              TSS.enrichment > 2)
}


# compute LSI
for (i in 1:length(subobjlist)){
  subobjlist[[i]] <- FindTopFeatures(subobjlist[[i]], min.cutoff = 10)
  subobjlist[[i]] <- RunTFIDF(subobjlist[[i]])
  subobjlist[[i]] <- RunSVD(subobjlist[[i]])
}

save(subobjlist, file = "after_QC_lsi_on_6libs.RData")

# rename the cells so that when merging, weird subscripts are not added
for (i in 1:length(subobjlist)){
  subobjlist[[i]] <- RenameCells(subobjlist[[i]], add.cell.id = paste0(i, "_"))
}

# first merge the six objects
combined <- merge(subobjlist[[1]] , subobjlist[[2]])
DefaultAssay(combined)
combined <- merge(combined, subobjlist[[3]])
combined <- merge(combined, subobjlist[[4]])
combined <- merge(combined, subobjlist[[5]])
combined <- merge(combined, subobjlist[[6]])
DefaultAssay(combined)

combined <- FindTopFeatures(combined, min.cutoff = 10)
combined <- RunTFIDF(combined)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:30)
DefaultAssay(combined)

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = subobjlist,
  anchor.features = rownames(combined),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
integrated <- FindNeighbors(object = integrated, 
                              reduction = 'integrated_lsi', 
                              dims = 2:30) # remove the first component

integrated <- FindClusters(object = integrated, 
                             verbose = FALSE, 
                             algorithm = 3, 
                             resolution = 0.5,
                             random.seed = 99)


DimPlot(object = integrated, label = TRUE) + 
  NoLegend() +
  ggtitle("integrated data colored by cluster")


dimplot_by_lib <- DimPlot(integrated, 
                          group.by = "orig.ident",
                          cols = "Set1") + 
  ggtitle("integrated data colored by libraries")
dimplot_by_lib[[1]]$layers[[1]]$aes_params$alpha = .3
dimplot_by_lib

unique(integrated$orig.ident) 
#"libraries_2_0" "libraries_2_1" "libraries_2_6" "libraries_8_0" "libraries_8_1" "libraries_8_6"

integrated$time.ident <- NA
integrated$time.ident[integrated$orig.ident %in% c("libraries_2_0", "libraries_8_0")] <- "0hr"
integrated$time.ident[integrated$orig.ident %in% c("libraries_2_1", "libraries_8_1")] <- "1hr"
integrated$time.ident[integrated$orig.ident %in% c("libraries_2_6", "libraries_8_6")] <- "6hr"
unique(integrated$time.ident)
# "0hr" "1hr" "6hr"

dimplot_by_time <- DimPlot(integrated, 
                           group.by = "time.ident", 
                           cols = c("cyan", "magenta", "gold")) + 
  ggtitle("integrated data colored by time points")
dimplot_by_time[[1]]$layers[[1]]$aes_params$alpha = .3
dimplot_by_time
DimPlot(integrated, 
        group.by = "time.ident", 
        cols = c("cyan", "transparent", "transparent")) + 
  ggtitle("integrated data - 0hr")
DimPlot(integrated, 
        group.by = "time.ident", 
        cols = c("transparent", "magenta", "transparent")) + 
  ggtitle("integrated data - 1hr")
DimPlot(integrated, 
        group.by = "time.ident", 
        cols = c("transparent", "transparent", "gold")) + 
  ggtitle("integrated data - 6hr")

save("integrated", file = "object_after_integration_UMAP_clustering.RData")



# assign cell type with RNAseq data ####

# load RNAseq object from Analysis_RNAseq_v3
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/nfeature5000_labeled_obj.RData")
head(four_celltype_obj@assays$RNA@data@Dimnames[[2]])
head(integrated@assays$macs2@data@Dimnames[[2]])

# assign cell types and clusters in RNAseq onto ATACseq
celltypes <- unique(RNAseq_integrated_labeled$broad.cell.type)
RNA_clusters <- unique(RNAseq_integrated_labeled$seurat_clusters)

integrated$broad.cell.type <- NA
for (i in 1:length(celltypes)){
  id <- celltypes[i]
  sub <- subset(RNAseq_integrated_labeled, broad.cell.type == id)
  integrated$broad.cell.type[str_sub(integrated@assays$macs2@data@Dimnames[[2]], start = 4L) %in%
                               str_sub(sub@assays$RNA@data@Dimnames[[2]], end = -3L)] <- id
}
unique(integrated$broad.cell.type)

integrated$RNA.clusters <- NA
for (i in 1:length(RNA_clusters)){
  c <- RNA_clusters[i]
  sub <- subset(RNAseq_integrated_labeled, seurat_clusters == c)
  integrated$RNA.clusters[str_sub(integrated@assays$macs2@data@Dimnames[[2]], start = 4L) %in%
                               str_sub(sub@assays$RNA@data@Dimnames[[2]], end = -3L)] <- c
}
unique(integrated$RNA.clusters)

# ATACseq clusters onto RNAseq obj
ATAC_clusters <- unique(integrated$seurat_clusters)

RNAseq_integrated_labeled$ATAC.clusters <- NA
for (i in 1:length(ATAC_clusters)){
  c <- ATAC_clusters[i]
  sub <- subset(integrated, seurat_clusters == c)
  RNAseq_integrated_labeled$ATAC.clusters[str_sub(RNAseq_integrated_labeled@assays$RNA@data@Dimnames[[2]], end = -3L) %in% 
                               str_sub(sub@assays$macs2@data@Dimnames[[2]], start = 4L)] <- c
}
unique(RNAseq_integrated_labeled$ATAC.clusters)

# check projection
DimPlot(object = RNAseq_integrated_labeled, 
        label = TRUE, 
        repel = T,
        group.by = "ATAC.clusters") + 
  ggtitle("RNAseq data \n projected by ATACseq clusters")

DimPlot(object = integrated, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.clusters") + 
  ggtitle("ATACseq data \n projected by RNAseq clusters")

DimPlot(object = integrated, 
        label = TRUE, 
        repel = T,
        group.by = "broad.cell.type") + 
  ggtitle("ATACseq data \n projected by RNAseq cell type ident")

Idents(integrated)

# check number of cells in each cell type
sum(integrated$broad.cell.type == "SEMA3E_pos_glut", na.rm = T) #SEMA3Ep glut: 5000
sum(integrated$broad.cell.type == "NPC", na.rm = T) #NPC: 4845
sum(integrated$broad.cell.type == "GABA", na.rm = T) #GABA: 16113
sum(integrated$broad.cell.type == "NEFM_pos_glut", na.rm = T) #NEFM_pos_glut: 16270
sum(integrated$broad.cell.type == "NEFM_neg_glut", na.rm = T) #NEFM_neg_glut: 13581
sum(integrated$broad.cell.type == "forebrain_NPC", na.rm = T) #forebrain_NPC: 2294
sum(integrated$broad.cell.type == "unknown", na.rm = T) #unknown: 10506



# differentially accessible peaks ####

# first add idents for cell type and time
celltypes <- unique(integrated$broad.cell.type)
integrated$celltype.time.ident <- NA
for (ct in celltypes){
  integrated$celltype.time.ident[integrated$time.ident == "0hr" &
                                   integrated$broad.cell.type == ct] <- paste0(ct, "0hr")
  integrated$celltype.time.ident[integrated$time.ident == "1hr" &
                                   integrated$broad.cell.type == ct] <- paste0(ct, "1hr")
  integrated$celltype.time.ident[integrated$time.ident == "6hr" &
                                   integrated$broad.cell.type == ct] <- paste0(ct, "6hr")
}
unique(integrated$celltype.time.ident)

integrated$spec.cell.type <- integrated$broad.cell.type
integrated$spec.cell.type[integrated$spec.cell.type %in% c("SEMA3E_pos_glut",
                                                            "forebrain_NPC",
                                                            "unknown")] <- "others"
celltypes <- sort(unique(integrated$spec.cell.type))
celltypes <- celltypes[celltypes != "others"]

# call differentially accessible peaks
Idents(integrated) <- "spec.cell.type"
j <- 0
# calc hour 1 vs 0
for (i in celltypes) {
  print(i)
  j = j + 1
  if (j == 1) {
    print("first cell type")
    df_1vs0 <- FindMarkers(object = integrated,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0,
                           test.use = "MAST")
    df_1vs0$peak <- rownames(df_1vs0)
    df_1vs0$cell_type <- i
  } else {
    print("others")
    df_to_append <- FindMarkers(object = integrated,
                                ident.1 = "1hr",
                                group.by = 'time.ident',
                                subset.ident = i,
                                ident.2 = "0hr",
                                min.pct = 0.0,
                                logfc.threshold = 0.0,
                                test.use = "MAST")
    df_to_append$peak <- rownames(df_to_append)
    df_to_append$cell_type <- i
    df_1vs0 <- rbind(df_1vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
  }
}

j <- 0
# calculate 6hr vs 0hr
for (i in celltypes) {
  print(i)
  j = j + 1
  if (j == 1) {
    print("first cell type")
    df_6vs0 <- FindMarkers(object = integrated,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_6vs0$peak <- rownames(df_6vs0)
    df_6vs0$cell_type <- i
  } else {
    print("others")
    df_to_append <- FindMarkers(object = integrated,
                                ident.1 = "6hr",
                                group.by = 'time.ident',
                                subset.ident = i,
                                ident.2 = "0hr",
                                min.pct = 0.0,
                                logfc.threshold = 0.0)
    df_to_append$peak <- rownames(df_to_append)
    df_to_append$cell_type <- i
    df_6vs0 <- rbind(df_6vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
  }
}
unique(df_6vs0$cell_type)
write.table(df_6vs0[df_6vs0$cell_type == "GABA", ], 
            file = "GABA_DE_markers_6v0hr.csv",
            col.names = T,
            row.names = T,
            sep = ",",
            quote = F)
