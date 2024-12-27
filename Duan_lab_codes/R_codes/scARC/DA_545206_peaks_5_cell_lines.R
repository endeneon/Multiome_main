# Siwei 23 Feb 2022
# run DA peaks and fragment counts using the 545206
# non-overlapping, unique peak set across 4 cell types
# for the 5-cell line data set

# init
library(Seurat)
library(Signac)

library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize = 200 * 1024 ^ 3)

library(readr)
library(GenomicRanges)
# library(En)


# load object

setwd("/nvmefs/scARC_Duan_018/R_peak_permutation")

load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/ATACseq_5line_objonly_labeled.RData")
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/annot_aggr_signac_EnsDb_UCSC.RData")

fragpath <- "~/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/atac_fragments.tsv.gz"

raw_peaks_all <-
  read_delim("~/NVME/scARC_Duan_018/R_peak_permutation/scARC_5_lines_summit/unique_non_overlap_peaks/sum_5_lines_all_unique_non_overlap_peaks.bed", 
             delim = "\t", escape_double = FALSE, 
             col_names = FALSE, trim_ws = TRUE)
# assign colnames to raw_peaks_all
colnames(raw_peaks_all) <- c("CHR", "START", "END", "NAME", "logQ", "STRAND")

# convert raw peaks to GRanges
all_peaks_GRanges <-
  makeGRangesFromDataFrame(df = raw_peaks_all,
                           keep.extra.columns = T,
                           ignore.strand = T,
                           seqnames.field = "CHR",
                           start.field = "START",
                           end.field = "END",
                           strand.field = "STRAND",
                           starts.in.df.are.0based = T)

# try recount fragments from the obj_complete
DefaultAssay(obj_complete)

non_overlapping_peaks_recounted <- 
  FeatureMatrix(fragments = Fragments(obj_complete),
                features = all_peaks_GRanges,
                cells = colnames(obj_complete))

save(file = "count_from_501_bp_peaks.RData",
     list = "non_overlapping_peaks_recounted",
     compression_level = 9)

# create a new object
Signac_501bp_recounted <-
  CreateChromatinAssay(counts = non_overlapping_peaks_recounted,
                       ranges = all_peaks_GRanges,
                       min.cells = 0,
                       min.features = 0)

Signac_501bp_recounted <-
  CreateSeuratObject(counts = Signac_501bp_recounted,
                     assay = "peaks")

Signac_501bp_recounted@assays$peaks@counts@Dimnames[[2]][1:10]
DefaultAssay(Signac_501bp_recounted)

# assign cell properties from Lexi's obj_complete
Signac_501bp_recounted$time.group.ident <- obj_complete$time.group.ident
Signac_501bp_recounted$time.ident <- obj_complete$time.ident
Signac_501bp_recounted$RNA.cluster.ident <- obj_complete$RNA.cluster.ident
Signac_501bp_recounted$trimmed.barcodes <- obj_complete$trimmed.barcodes
Signac_501bp_recounted$RNA.cell.type.ident <- obj_complete$RNA.cell.type.ident
Signac_501bp_recounted$broad.cell.type <- obj_complete$broad.cell.type
Signac_501bp_recounted$fine.cell.type <- obj_complete$fine.cell.type
Signac_501bp_recounted$celltype.time.ident <- obj_complete$celltype.time.ident
Signac_501bp_recounted$broad.celltype.time.ident <- obj_complete$broad.celltype.time.ident
###

rm(annot_aggr_signac_ucsc)
rm(obj_complete)
rm(non_overlapping_peaks_recounted)
rm(raw_peaks_all)

save.image(file = "Signac_501bp_raw_object.RData")

# run normalisation
Signac_501bp_recounted <- RunTFIDF(Signac_501bp_recounted)
Signac_501bp_recounted <- FindTopFeatures(Signac_501bp_recounted)
Signac_501bp_recounted <- RunSVD(object = Signac_501bp_recounted)

DepthCor(Signac_501bp_recounted)

Signac_501bp_recounted <- RunUMAP(object = Signac_501bp_recounted,
                                  reduction = 'lsi',
                                  dims = 2:30)
Signac_501bp_recounted <- FindNeighbors(object = Signac_501bp_recounted,
                                        reduction = 'lsi',
                                        dims = 2:30)
Signac_501bp_recounted <- FindClusters(object = Signac_501bp_recounted,
                                       algorithm = 3,
                                       resolution = 1.2,
                                       verbose = T)
DimPlot(Signac_501bp_recounted, label = F)
Idents(Signac_501bp_recounted) <- "time.ident"

write.table(colnames(Signac_501bp_recounted),
            file = "scATAC_cell_list_from_Lexi.txt",
            quote = F, row.names = F, col.names = F, sep = "\t")

colnames(Signac_501bp_recounted)[1:10]
Signac_501bp_recounted$time.group.ident[1:10]
Signac_501bp_recounted$time.ident[1:10]
Signac_501bp_recounted$trimmed.barcodes[1:10]

##### re-build the integrated object from 6 ATAC-seq libraries
library(stringr)

# init ####
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
# library()
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)

library(stringr)
library(ggplot2)
library(readr)
library(RColorBrewer)
# library(dplyr)
library(viridis)
library(graphics)

library(future)

set.seed(42)
# set threads and parallelization
plan("multisession", workers = 4)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 100 * 1024 ^ 3)
options(ucscChromosomeNames = T)

# load data #### 
setwd("~/NVME/scARC_Duan_018/GRCh38_mapped_only")

load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/annot_aggr_signac_EnsDb_UCSC.RData")
rm(annot_aggr_signac)

# Read files, separate the ATACseq data from the .h5 matrix
# h5list <- list.dirs(path = ".", 
#                      # pattern = "outs$",
#                      recursive = F,
#                      include.dirs = T)
h5list <- list.dirs(path = ".", 
                    # pattern = "outs$",
                    recursive = F)
h5list <- h5list[str_detect(string = h5list,
                            pattern = "libraries")]
h5list <- paste(h5list, 'outs/filtered_feature_bc_matrix', 
                sep = '/')

gzlist <- list.files(path = ".", 
                     pattern = "atac_fragments.tsv.gz$",
                     recursive = T,
                     include.dirs = T)
# objlist <- vector(mode = "list", length = length(h5list))

barcode_library_list <- vector(mode = "list", length = length(h5list))
names(barcode_library_list) <- c("g_2_0", "g_2_1", "g_2_6",
                                 "g_8_0", "g_8_1", "g_8_6")



pseudo_bulk_unique_non_overlap_dataset <- 
  vector(mode = "list", length = length(h5list))


i <- 1
for (i in 1:6) {
  # direct count the fragment file from here

  # h5file <- Read10X_h5(filename = h5list[i])
  h5file <- Read10X(data.dir = h5list[i])
  
  chrom_assay_get_cell_list <-
    as.character(h5file$Peaks@Dimnames[[2]])
  
  # chrom_assay_get_cell_list <- CreateChromatinAssay(
  #   counts = h5file$Peaks,
  #   sep = c(":", "-"),
  #   genome = 'hg38',
  #   fragments = gzlist[i],
  #   min.cells = 10,
  #   min.features = 200
  # )
  
  # record cell barcodes per group
  barcode_library_list[[i]] <- chrom_assay_get_cell_list
    # colnames(chrom_assay_get_cell_list)
  # create fragment file
  chrom_frag <- 
    CreateFragmentObject(path = gzlist[i],
                         cells = chrom_assay_get_cell_list)
  
  # create featureMatrix
  chrom_fmatrix <-
    FeatureMatrix(fragments = chrom_frag,
                  features = all_peaks_GRanges,
                  cells = chrom_assay_get_cell_list)
  
  print(i)
  # create chromAssay
  pseudo_bulk_chromAssay <- 
    CreateChromatinAssay(counts = chrom_fmatrix,
                         # fragments = chrom_frag, # need a Fragment object here
                         annotation = annot_aggr_signac_ucsc,
                         sep = c("-", "-"),
                         genome = "hg38",
                         min.cells = 10,
                         min.features = 200)
  # colnames(pseudo_bulk_chromAssay) <- chrom_assay_get_cell_list
  print(i)
  # assign chromAssay to objList as a list
  pseudo_bulk_unique_non_overlap_dataset[[i]] <- 
    CreateSeuratObject(
    counts = pseudo_bulk_chromAssay,
    names.field = 1,
    names.delim = "-",
    assay = "pseudo_bulk_ATAC",
    project = str_extract(string = h5list[i], 
                          pattern = "libraries_[0-9]_[0-6]"))
  
  # this object needs to have its fragment object assigned separately
  pseudo_bulk_unique_non_overlap_dataset[[i]]@assays$pseudo_bulk_ATAC@fragments[[1]] <-
    CreateFragmentObject(path = gzlist[i],
                         cells = chrom_assay_get_cell_list)
  
}

setwd("~/NVME/scARC_Duan_018/R_peak_permutation")
rm(chrom_fmatrix)
rm(chrom_frag)
rm(pseudo_bulk_chromAssay)
rm(h5file)
rm(chrom_assay_get_cell_list)

save.image(file = "raw_reads_directly_540k_peaks_6_objects_28Feb2022.RData")

#### normalise 6 libraries independently

load(file = "raw_reads_directly_540k_peaks_6_objects_28Feb2022.RData")

# change working directory to count fragment files and filter
setwd("~/NVME/scARC_Duan_018/GRCh38_mapped_only")

i <- 1
for (i in 1:6) {
  print(i)
  DefaultAssay(pseudo_bulk_unique_non_overlap_dataset[[i]]) <- "pseudo_bulk_ATAC"
  
  Annotation(pseudo_bulk_unique_non_overlap_dataset[[i]]) <- 
    annot_aggr_signac_ucsc
  
  pseudo_bulk_unique_non_overlap_dataset[[i]] <- 
    TSSEnrichment(object = pseudo_bulk_unique_non_overlap_dataset[[i]], 
                  fast = F)
  pseudo_bulk_unique_non_overlap_dataset[[i]] <- 
    NucleosomeSignal(object = pseudo_bulk_unique_non_overlap_dataset[[i]])

  temp_frag <- CountFragments(fragments = gzlist[i],
                              cells = barcode_library_list[[i]])
  temp_count_peaks <- 
    CountsInRegion(object = pseudo_bulk_unique_non_overlap_dataset[[i]],
                   assay = "pseudo_bulk_ATAC",
                   regions = pseudo_bulk_unique_non_overlap_dataset[[i]]@assays$pseudo_bulk_ATAC@ranges)
  
  pseudo_bulk_unique_non_overlap_dataset[[i]]$peak_region_fragment <- 
    temp_count_peaks
  temp_cells <-
    data.frame(CB = names(pseudo_bulk_unique_non_overlap_dataset[[i]]$nCount_pseudo_bulk_ATAC))
  temp_cells <-
    merge(temp_cells,
          temp_frag,
          by = "CB")
  pseudo_bulk_unique_non_overlap_dataset[[i]]$reads_count <- 
    temp_cells$reads_count

  pseudo_bulk_unique_non_overlap_dataset[[i]]$pct_reads_in_peaks <-
    pseudo_bulk_unique_non_overlap_dataset[[i]]$peak_region_fragment / 
    pseudo_bulk_unique_non_overlap_dataset[[i]]$reads_count * 100
  
  pseudo_bulk_unique_non_overlap_dataset[[i]] <- 
    subset(pseudo_bulk_unique_non_overlap_dataset[[i]], 
           nCount_pseudo_bulk_ATAC > 1500 & 
                   nCount_pseudo_bulk_ATAC < 10000 &
                   nucleosome_signal < 1 &
                   TSS.enrichment > 2 &
             pct_reads_in_peaks > 18 &
             pct_reads_in_peaks < 42)
}

setwd("/nvmefs/scARC_Duan_018/R_peak_permutation")
save.image(file = "raw_reads_540k_peaks_6_objects_filtered_01Mar2022.RData")

#### calculate the LSI of each element
load(file = "raw_reads_540k_peaks_6_objects_filtered_01Mar2022.RData")

i <- 1L
# compute LSI
for (i in 1:length(pseudo_bulk_unique_non_overlap_dataset)) {
  print(i)
  pseudo_bulk_unique_non_overlap_dataset[[i]] <- 
    FindTopFeatures(pseudo_bulk_unique_non_overlap_dataset[[i]], 
                    min.cutoff = "q90", 
                    verbose = T)
  pseudo_bulk_unique_non_overlap_dataset[[i]] <- 
    RunTFIDF(pseudo_bulk_unique_non_overlap_dataset[[i]], 
             method = 1, 
             verbose = T)
  pseudo_bulk_unique_non_overlap_dataset[[i]] <- 
    RunSVD(pseudo_bulk_unique_non_overlap_dataset[[i]], 
           features = pseudo_bulk_unique_non_overlap_dataset[[i]]@assays$pseudo_bulk_ATAC@var.features,
           # n = 5000,
           verbose = T)
}

# rename the cells so that when merging, weird subscripts are not added
for (i in 1:length(pseudo_bulk_unique_non_overlap_dataset)) {
  print(i)
  pseudo_bulk_unique_non_overlap_dataset[[i]] <- 
    RenameCells(pseudo_bulk_unique_non_overlap_dataset[[i]], 
                add.cell.id = paste0(i, "_"))
}



# make the combined object
i <- 1L
for (i in 1:length(pseudo_bulk_unique_non_overlap_dataset)) {
  print(i)
  if (i == 1) {
    pseudo_bulk_merged <- pseudo_bulk_unique_non_overlap_dataset[[i]]
  } else {
    pseudo_bulk_merged <-
      merge(pseudo_bulk_merged,
            pseudo_bulk_unique_non_overlap_dataset[[i]])
  }
}

save(list = "pseudo_bulk_merged",
     file = "merged_540k_peaks_6_objects_filtered_01Mar2022.RData")

# process the merged dataset
pseudo_bulk_merged <- 
  FindTopFeatures(pseudo_bulk_merged,
                  min.cutoff = "q90",
                  verbose = T)
pseudo_bulk_merged <-
  RunTFIDF(pseudo_bulk_merged,
           method = 1,
           verbose = T)
pseudo_bulk_merged <-
  RunSVD(pseudo_bulk_merged,
         features = pseudo_bulk_merged@assays$pseudo_bulk_ATAC@var.features,
         verbose = T)
pseudo_bulk_merged <-
  RunUMAP(pseudo_bulk_merged,
          reduction = 'lsi',
          dims = 2:30)
save(list = "pseudo_bulk_merged",
     file = "merged_540k_peaks_6_objects_filtered_01Mar2022.RData")

integration_anchors <-
  FindIntegrationAnchors(object.list = pseudo_bulk_unique_non_overlap_dataset,
                         anchor.features = rownames(pseudo_bulk_merged),
                         reduction = 'rlsi',
                         dims = 2:30)

pseudo_bulk_integrated <-
  IntegrateEmbeddings(anchorset = integration_anchors,
                      reductions = pseudo_bulk_merged[["lsi"]],
                      new.reduction.name = "integrated_lsi",
                      dims.to.integrate = 1:30)
pseudo_bulk_integrated <-
  RunUMAP(object = pseudo_bulk_integrated,
          reduction = 'integrated_lsi',
          dims = 2:30, 
          verbose = T)

DimPlot(pseudo_bulk_integrated,
        group.by = "orig.ident",
        cells = colnames(pseudo_bulk_integrated)[pseudo_bulk_integrated$orig.ident %in%
                                                            "libraries_2_0"])
DimPlot(pseudo_bulk_integrated,
        split.by = "orig.ident",
        ncol = 3) +
  NoLegend()

colnames(pseudo_bulk_integrated)

# assign cell type with RNAseq data ####

# load RNAseq object from Analysis_RNAseq_v3 and project the identities from RNAseq

# change cell barcodes
library(dplyr)
library(stringr)

lib_orig_ident <-
  unique(pseudo_bulk_integrated$orig.ident)

unified_colnames <-
  colnames(pseudo_bulk_integrated) %>%
  str_sub(start = 4L, end = -2L) %>%
  paste0(., i)
unified_colnames
length(unified_colnames)
length(colnames(pseudo_bulk_integrated))

pseudo_bulk_integrated$orig.colnames <- 
  colnames(pseudo_bulk_integrated)

colnames(pseudo_bulk_integrated) <- unified_colnames
lib_name <- ""
i <- 1L

pseudo_bulk_integrated$corrected.barcodes <- NA
for (i in 1:length(lib_orig_ident)) {
  pseudo_bulk_integrated$corrected.barcodes[pseudo_bulk_integrated$orig.ident %in%
                                     lib_orig_ident[i]] <-
    colnames(pseudo_bulk_integrated)[pseudo_bulk_integrated$orig.ident %in%
                                       lib_orig_ident[i]] %>%
    str_sub(start = 4L, end = -2L) %>%
    paste0(., i)
}



# load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/nfeature5000_labeled_obj.RData")
# load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/10Feb2022_increased_nfeat_after_transformation.RData")
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/human_only_obj_with_cell_type_labeling.RData")
four_celltype_obj <- human_only
head(four_celltype_obj@assays$RNA@data@Dimnames[[2]])
head(pseudo_bulk_integrated$orig.ident)
head(pseudo_bulk_integrated$corrected.barcodes)

unique(str_split(colnames(four_celltype_obj), 
                 pattern = "-",
                 simplify = T)[, 2])
ncol(four_celltype_obj)

four_celltype_obj$corrected.barcodes <-
  paste0(str_sub(colnames(four_celltype_obj),
                 end = -4),
         str_sub(colnames(four_celltype_obj),
                 start = -1))
four_celltype_obj$corrected.barcodes

# assign cell types and clusters in RNAseq onto ATACseq
celltypes <- unique(four_celltype_obj$cell.type.for.plot)
RNA_clusters <- unique(four_celltype_obj$seurat_clusters)

# pseudo_bulk_integrated$broad.cell.type <- NA
pseudo_bulk_integrated$cell.types <- NA
for (i in 1:length(celltypes)) {
  # id <- celltypes[i]
  pseudo_bulk_integrated$cell.types[pseudo_bulk_integrated$corrected.barcodes %in%
                                           four_celltype_obj$corrected.barcodes[four_celltype_obj$cell.type.for.plot %in%
                                                                                  celltypes[i]]] <- 
                                           celltypes[i]
}

celltypes <- unique(pseudo_bulk_integrated$broad.cell.type)

pseudo_bulk_integrated$RNA.clusters <- NA

for (i in 1:length(RNA_clusters)) {
  pseudo_bulk_integrated$RNA.clusters[pseudo_bulk_integrated$corrected.barcodes %in%
                            four_celltype_obj$corrected.barcodes[four_celltype_obj$seurat_clusters %in% as.character(i)]] <- 
    as.character(i)
}

pseudo_bulk_integrated$human.origin <- T
pseudo_bulk_integrated$human.origin[is.na(pseudo_bulk_integrated$cell.types)] <- F


DimPlot(object = pseudo_bulk_integrated, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.clusters") + 
  NoLegend()

DimPlot(object = pseudo_bulk_integrated, 
        label = TRUE, 
        repel = T,
        group.by = "cell.types") + 
  NoLegend()

DimPlot(object = Integrated_540k_peaks_human_only, 
        label = TRUE, 
        repel = T,
        group.by = "cell.types") + 
  NoLegend()

# ggtitle("ATACseq data \n projected by RNAseq clusters")


Integrated_540k_peaks_human_only <-
  subset(pseudo_bulk_integrated,
         subset = human.origin == T)

save(list = "Integrated_540k_peaks_human_only",
     file = "Integrated_540k_peaks_human_only.RData",
     compression_level = 9)

save.image(file = "Integrated_540k_peaks_human_only.RData")

load("Integrated_540k_peaks_human_only.RData")  
### work with Integrate_540k

DimPlot(Integrated_540k_peaks_human_only,
  reduction = "umap",
  group.by = "seurat_clusters")
DimPlot(Integrated_540k_peaks_human_only,
        reduction = "umap",
        group.by = "cell.types")


Integrated_540k_peaks_human_only$broad.cell.type.ATAC <- NA
Integrated_540k_peaks_human_only$broad.cell.type.ATAC[Integrated_540k_peaks_human_only$seurat_clusters %in%
                                                        c("1", "6", "7")] <- "NEFMp_glut"
Integrated_540k_peaks_human_only$broad.cell.type.ATAC[Integrated_540k_peaks_human_only$seurat_clusters %in%
                                                        c("2", "3", "10")] <- "NEFMm_glut"
Integrated_540k_peaks_human_only$broad.cell.type.ATAC[Integrated_540k_peaks_human_only$seurat_clusters %in%
                                                        c("0", "4")] <- "GABA_ATAC"
Integrated_540k_peaks_human_only$broad.cell.type.ATAC[Integrated_540k_peaks_human_only$cell.types %in%
                                                        "NPC"] <- "NPC"
unique(Integrated_540k_peaks_human_only$broad.cell.type.ATAC)


Integrated_540k_peaks_human_only <-
  FindNeighbors(object = Integrated_540k_peaks_human_only,
                reduction = 'integrated_lsi',
                dims = 2:30)

Integrated_540k_peaks_human_only <-
  FindClusters(object = Integrated_540k_peaks_human_only,
               algorithm = 2,
               resolution = 1.0,
               verbose = T,
               random.seed = 42)
DimPlot(object = Integrated_540k_peaks_human_only,
        label = T) +
  NoLegend()

# add timing information
Integrated_540k_peaks_human_only$sample.time <-
  paste0(unlist(str_split(Integrated_540k_peaks_human_only$orig.ident,
                          pattern = "_",
                          simplify = T)[, 3]),
         "hr")

# Find peaks
Idents(Integrated_540k_peaks_human_only) <- "broad.cell.type.ATAC"

DA_540k_1vs0_GABA <- 
  FindMarkers(object = Integrated_540k_peaks_human_only,
              # slot = "pseudo_bulk_ATAC",
              ident.1 = "1hr",
              ident.2 = "0hr",
              group.by = "sample.time",
              subset.ident = "GABA_ATAC",
              # cells.1 = colnames(Integrated_540k_peaks_human_only)[(Integrated_540k_peaks_human_only$sample.time %in% "1hr") &
              #                                                        (Integrated_540k_peaks_human_only$seurat_clusters %in% c(1, 6, 7))],
              # cells.2 = colnames(Integrated_540k_peaks_human_only)[(Integrated_540k_peaks_human_only$sample.time %in% "0hr") &
              #                                                        (Integrated_540k_peaks_human_only$seurat_clusters %in% c(1, 6, 7))],
              logfc.threshold = 0.01,
              min.pct = 0.01,
              random.seed = 42,
              min.cells.feature = 10,
              test.use = "MAST")
# DA_540k_1vs0_GABA <- DA_540k_0vs1_GABA
DA_540k_6vs0_GABA <- 
  FindMarkers(object = Integrated_540k_peaks_human_only,
              ident.1 = "6hr",
              ident.2 = "0hr",
              group.by = "sample.time",
              subset.ident = "GABA_ATAC",
              logfc.threshold = 0.01,
              min.pct = 0.01,
              random.seed = 42,
              min.cells.feature = 10,
              test.use = "MAST")
# save.image()

DA_540k_1vs0_NEFMp_Glut <- 
  FindMarkers(object = Integrated_540k_peaks_human_only,
              ident.1 = "1hr",
              ident.2 = "0hr",
              group.by = "sample.time",
              subset.ident = "NEFMp_glut",
              logfc.threshold = 0.01,
              min.pct = 0.01,
              random.seed = 42,
              min.cells.feature = 10,
              test.use = "MAST")

DA_540k_6vs0_NEFMp_Glut <- 
  FindMarkers(object = Integrated_540k_peaks_human_only,
              ident.1 = "6hr",
              ident.2 = "0hr",
              group.by = "sample.time",
              subset.ident = "NEFMp_glut",
              logfc.threshold = 0.01,
              min.pct = 0.01,
              random.seed = 42,
              min.cells.feature = 10,
              test.use = "MAST")
save.image()

DA_540k_1vs0_NEFMm_Glut <- 
  FindMarkers(object = Integrated_540k_peaks_human_only,
              ident.1 = "1hr",
              ident.2 = "0hr",
              group.by = "sample.time",
              subset.ident = "NEFMm_glut",
              logfc.threshold = 0.01,
              min.pct = 0.01,
              random.seed = 42,
              min.cells.feature = 10,
              test.use = "MAST")

DA_540k_6vs0_NEFMm_Glut <- 
  FindMarkers(object = Integrated_540k_peaks_human_only,
              ident.1 = "6hr",
              ident.2 = "0hr",
              group.by = "sample.time",
              subset.ident = "NEFMm_glut",
              logfc.threshold = 0.01,
              min.pct = 0.01,
              random.seed = 42,
              min.cells.feature = 10,
              test.use = "MAST")  
save.image()

DA_540k_1vs0_NPC <- 
  FindMarkers(object = Integrated_540k_peaks_human_only,
              ident.1 = "1hr",
              ident.2 = "0hr",
              group.by = "sample.time",
              subset.ident = "NPC",
              logfc.threshold = 0.01,
              min.pct = 0.01,
              random.seed = 42,
              min.cells.feature = 10,
              test.use = "MAST")

DA_540k_6vs0_NPC <- 
  FindMarkers(object = Integrated_540k_peaks_human_only,
              ident.1 = "6hr",
              ident.2 = "0hr",
              group.by = "sample.time",
              subset.ident = "NPC",
              logfc.threshold = 0.01,
              min.pct = 0.01,
              random.seed = 42,
              min.cells.feature = 10,
              test.use = "MAST")

write.table(DA_540k_1vs0_NPC,
            file = "output_540k/DA_540k_1vs0_NPC.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)
write.table(DA_540k_6vs0_NPC,
            file = "output_540k/DA_540k_6vs0_NPC.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)




write.table(DA_540k_1vs0_GABA,
            file = "output_540k/DA_540k_1vs0_GABA.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)
write.table(DA_540k_6vs0_GABA,
            file = "output_540k/DA_540k_6vs0_GABA.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)
write.table(DA_540k_6vs0_NEFMm_Glut,
            file = "output_540k/DA_540k_6vs0_NEFMm.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)
write.table(DA_540k_1vs0_NEFMm_Glut,
            file = "output_540k/DA_540k_1vs0_NEFMm.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)
write.table(DA_540k_1vs0_NEFMp_Glut,
            file = "output_540k/DA_540k_1vs0_NEFMp.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)
write.table(DA_540k_6vs0_NEFMp_Glut,
            file = "output_540k/DA_540k_6vs0_NEFMp.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)

library(stringr)
library(ggplot2)
sum(str_detect(string = raw_peaks_all$NAME,
               pattern = "GABA")) # 108013
sum(str_detect(string = raw_peaks_all$NAME,
               pattern = "NEFMm")) # 96116
sum(str_detect(string = raw_peaks_all$NAME,
               pattern = "NEFMp")) # 187782
sum(str_detect(string = raw_peaks_all$NAME,
               pattern = "NPC")) # 153295

DA_540k_1vs0_GABA$PeakID <- rownames(DA_540k_6vs0_NPC)
temp_table <- 
  read_delim("output_540k/output/DA_540k_6vs0_NPC.annot.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
DA_540k_1vs0_GABA <- merge(x = DA_540k_1vs0_GABA,
                           y = temp_table,
                           by.x = "PeakID",
                           by.y = 'PeakID (cmd=annotatePeaks.pl temp.bed hg38)')
DA_540k_1vs0_GABA$respond <- DA_540k_1vs0_GABA$`Gene Name`

DA_540k_1vs0_GABA$respond[!(DA_540k_1vs0_GABA$`Gene Name` %in%
                                  c("FOS", "NPAS4", "IGF1", "NR4A1",
                                    "FOS2", "EGR1", "BDNF", "FOSB") &
                              (str_detect(string = DA_540k_1vs0_GABA$Annotation,
                                          pattern = c("promoter-TSS"),
                                          negate = F)))] <- ""

# View((DA_540k_1vs0_GABA)[(DA_540k_1vs0_GABA$`Gene Name` %in%
#                             c("FOS", "NPAS4", "IGF1", "PNOC")), ]) 
View(DA_540k_1vs0_GABA[DA_540k_1vs0_GABA$respond != "", ])
# make volcano plots
library(ggplot2)
library(ggrepel)

df <- DA_540k_6vs0_NPC
# temp_table$`PeakID (cmd=annotatePeaks.pl temp.bed hg38)`

df$PeakID <- rownames(df)
temp_table <- 
  read_delim("output_540k/output/DA_540k_6vs0_NPC.annot.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)
df <- merge(x = df,
            y = temp_table,
            by.x = "PeakID",
            by.y = 'PeakID (cmd=annotatePeaks.pl temp.bed hg38)')
df$respond <- df$`Gene Name`

df$respond[!(df$`Gene Name` %in%
               c("FOS", "NPAS4", "IGF1", "NR4A1",
                 "FOSB", "BDNF") &
               (str_detect(string = df$Annotation,
                           pattern = c("Intergenic|non-coding"),
                           negate = T)))] <- ""

# View((DA_540k_1vs0_GABA)[(DA_540k_1vs0_GABA$`Gene Name` %in%
#                             c("FOS", "NPAS4", "IGF1", "PNOC")), ]) 
# View(df[df$respond != "", ])
ggplot(df, 
       aes(x = avg_log2FC,
           y = (0 - log10(p_val)),
           label = respond,
           color = ifelse(p_val_adj > 0.05,
                          yes = "FDR > 0.05",
                          no = ifelse(avg_log2FC > 0,
                                      yes = "log2FC > 0,\n FDR < 0.05",
                                      no = "log2FC < 0,\n FDR < 0.05")))) +
  scale_colour_manual(values = c("grey70", "steelblue", "red")) +
  labs(colour = "") +
  geom_point(shape = 20,
             size = 0.5) +
  ylab('-log10P') +
  ylim(0, 100) +
  theme_bw() +
  theme(legend.text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  geom_text_repel(max.overlaps = Inf,
                  # box.padding = 1,
                  min.segment.length = 0.5,
                  show.legend = F,
                  force = 500) +
  ggtitle(label = paste("500 bp peaks from",
                        "NPC",
                        "6hr vs 0hr,\n",
                        sum(df$p_val_adj < 0.05),
                        "DA peaks in total"))


###
# write.table(df,
#             file = "output_540k/DA_540k_1vs0_NPC_annotated.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)
write.table(df,
            file = "output_540k/DA_540k_6vs0_NPC_annotated.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)
# write.table(df,
#             file = "output_540k/DA_540k_1vs0_GABA_annotated.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)
# write.table(df,
#             file = "output_540k/DA_540k_6vs0_GABA_annotated.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)
# write.table(df,
#             file = "output_540k/DA_540k_6vs0_NEFMm_annotated.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)
# write.table(df,
#             file = "output_540k/DA_540k_1vs0_NEFMm_annotated.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)
# write.table(df,
#             file = "output_540k/DA_540k_1vs0_NEFMp_annotated.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)
# write.table(df,
#             file = "output_540k/DA_540k_6vs0_NEFMp_annotated.tsv",
#             quote = F, sep = "\t",
#             row.names = T, col.names = T)

# check barcode consistency
i <- 1L
for (i in 1:6) {
  print(i)
  print(sum((pseudo_bulk_integrated$corrected.barcodes[pseudo_bulk_integrated$orig.ident %in%
                                                         lib_orig_ident[i]]) %in% 
              (colnames(four_celltype_obj)[four_celltype_obj$time.group.ident %in%
                                             as.character(i)])))
}






save.image()
# find integration anchors


# sample(x = names(test$orig.ident),
#        size = 200)

### test FindTopFeatures, etc.
test <- pseudo_bulk_unique_non_overlap_dataset[[1]]
# Idents(test)
# test <- subset(test, 
#                cells = sample(x = names(test$orig.ident),
#                               size = 20000))

test <-
  FindTopFeatures(test, 
                  min.cutoff = "q90", 
                  verbose = T)
length(test@assays$pseudo_bulk_ATAC@var.features) # 1000:275648; 2000:293584

#######
# look at the objects first, before subsetting
test <- pseudo_bulk_unique_non_overlap_dataset[[1]]
test2 <- pseudo_bulk_unique_non_overlap_dataset[[4]]

DefaultAssay(test) <- "pseudo_bulk_ATAC"

test@assays$pseudo_bulk_ATAC@fragments[[1]] <-
  CreateFragmentObject(path = gzlist[i],
                       cells = barcode_library_list[[1]])
DefaultAssay(test2) <- "pseudo_bulk_ATAC"
hist(test$nCount_pseudo_bulk_ATAC, breaks = 1000)

# load annotation file because annotations cannot be directly accessed using GetGRangesFromEnsDb
# test with time=0 first
test@assays$pseudo_bulk_ATAC@seqinfo@genome <- "GRCh38"
Annotation(test) <- 
  GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
EnsDb_Hs_v86 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# Annotation(test2) <- annot_aggr_signac_ucsc
Annotation(test) <- annot_aggr_signac_ucsc

DefaultAssay(test)
DefaultAssay(test2)

test <- TSSEnrichment(object = test, fast = F)
test <- NucleosomeSignal(object = test)
test$perc_reads_in_peak
test_frag <- CountFragments(fragments = gzlist[1],
                       cells = barcode_library_list[[1]])
test_count_peaks <- CountsInRegion(object = test,
                                   assay = "pseudo_bulk_ATAC",
                                   regions = test@assays$pseudo_bulk_ATAC@ranges)
test$peak_region_fragment <- test_count_peaks
test_cells <-
  data.frame(CB = names(test$nCount_pseudo_bulk_ATAC))
test_cells <-
  merge(test_cells,
        test_frag,
        by = "CB")
test$reads_count <- test_cells$reads_count
test$nCount_pseudo_bulk_ATAC
test$pct_reads_in_peaks <-
  test$peak_region_fragment / test$reads_count * 100
# test_FRIP <- FRiP(object = test,
#              assay = 'pseudo_bulk_ATAC',
#              total.fragments = "fragments",
#              col.name = "FRIP_perc")
annot_aggr_signac_ucsc

test$pct_reads_in_peaks <-
  test$
test$high.tss <- ifelse(test$TSS.enrichment > 2, 
                        'High', 'Low')
TSSPlot(test, group.by = 'high.tss') +
  NoLegend()
TSSPlot(test)

VlnPlot(object = test,
        features = c('pct_reads_in_peaks', 'nCount_pseudo_bulk_ATAC',
                     'TSS.enrichment', 'nucleosome_signal'),
        pt.size = 0.1,
        ncol = 4)
hist(test$pct_reads_in_peaks)
hist(test$nCount_pseudo_bulk_ATAC, 
     breaks = 1000, 
     xlim = c(0, 10000))
hist(test$nucleosome_signal, 
     breaks = 100)

test2 <- TSSEnrichment(object = test2, fast = T)
test <- NucleosomeSignal(object = test)
test <- subset(test, nCount_pseudo_bulk_ATAC > 2000 & 
                 nCount_pseudo_bulk_ATAC < 30000 &
                 nucleosome_signal < 4 &
                 TSS.enrichment > 2)
test2 <- NucleosomeSignal(object = test2)
test2 <- subset(test2, nCount_peaks > 2000 & 
                  nCount_peaks < 30000 &
                  nucleosome_signal < 4 &
                  TSS.enrichment > 2)




setwd("/nvmefs/scARC_Duan_018/R_peak_permutation")
