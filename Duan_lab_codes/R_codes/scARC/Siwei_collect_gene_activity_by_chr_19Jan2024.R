# Siwei 19 Jan 2024
# Use Lexi's multiomic data of 28 Apr 2023 (new peakset)
# Collect Gene activity scores and plot graphs accordingly
# Count gene activity matrix from FeatureMatrix
# Make FeatureMatrix in blocks, segmented by chromosomes
# (memory index exhausted if count as one)


# init ####
library(Seurat)
library(Signac)
# library(sctransform)
# library(glmGamPoi)
library(EnsDb.Hsapiens.v86)
# library(GenomeInfoDb)
# library(GenomicFeatures)
# library(AnnotationDbi)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

# library(patchwork)
# library(readr)
# library(ggplot2)
library(RColorBrewer)
# library(dplyr)
# library(viridis)
# library(graphics)
library(stringr)
library(future)
library(Matrix)

plan("multisession", workers = 6)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# load data #####
# get peaksets
hg38_gene_activity_10k_1k <-
  readRDS("gRanges_hg38_gene_activity_10k_1k_19Jan2024.RDs")
hg38_gene_activity_10k_1k <-
  sortSeqlevels(hg38_gene_activity_10k_1k)

## split into a GRangesList
List_hg38_gene_activity_10k_1k <-
  split(hg38_gene_activity_10k_1k,
        f = seqnames(hg38_gene_activity_10k_1k))
names(List_hg38_gene_activity_10k_1k)

test_gRangesList <-
  GRangesList(gr1 = hg38_gene_activity_10k_1k[1],
              gr2 = hg38_gene_activity_10k_1k[2])

# load cell list
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/df_Lexi_cells_470K.RData")

# load fragment file
# fragments_raw <-
#   readRDS(file = "Lexi_fragments_raw_obj.RDs")
cell_list_downsampled_10K <-
  sample(x = df_Lexi_cells_470K,
         size = 10000,
         replace = F)
fragments_raw <-
  CreateFragmentObject(path = "./018_to_029_combined/outs/atac_fragments.tsv.gz",
                       cells = cell_list_downsampled_10K,
                       validate.fragments = F,
                       verbose = T)

system.time({
  raw_featureMatrix <-
    FeatureMatrix(fragments = fragments_raw,
                  features = hg38_gene_activity_10k_1k,
                  process_n = 4000,
                  verbose = T)
})
save(raw_featureMatrix,
     file = "raw_featureMatrix_10K_cell.RData")

## 
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/raw_featureMatrix_10K_cell.RData")

annotations_hg38 <-
  GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) 

seqlevels(annotations_hg38) <-
  paste0("chr", seqlevels(annotations_hg38))
annotations_hg38 <-
  sortSeqlevels(annotations_hg38)
genome(annotations_hg38) <-
  "hg38"


# make gene activity assay
subsetted_gene_activity_assay <-
  CreateAssayObject(counts = raw_featureMatrix,
                    annotation = annotations_hg38)

subsetted_gene_activity_assay <-
  CreateSeuratObject(counts = subsetted_gene_activity_assay,
                     assay = "gact")
subsetted_gene_activity_assay <-
  NormalizeData(object = subsetted_gene_activity_assay,
                normalization.method = "LogNormalize",
                scale.factor = median(subsetted_gene_activity_assay$nCount_gact),
                verbose = T)
subsetted_gene_activity_assay$barcode_wo_suffix <-
  str_split(string = unique(colnames(subsetted_gene_activity_assay)),
            pattern = "-",
            simplify = T)[, 1]


DefaultAssay(subsetted_gene_activity_assay) <- "gact"

# load the original RNA object (550K) to get dimensions ####
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/018-029_RNA_integrated_labeled_with_harmony.RData")
integrated_labeled$barcode_wo_suffix <-
  str_split(string = unique(colnames(integrated_labeled)),
            pattern = "-",
            simplify = T)[, 1]

## count cells in embeddings (UMAP)
DefaultAssay(integrated_labeled) <-
  "integrated"
sum(subsetted_gene_activity_assay$barcode_wo_suffix %in%
      integrated_labeled$barcode_wo_suffix)

integrated_labeled_subsetted <-
  integrated_labeled_subsetted[, integrated_labeled$barcode_wo_suffix %in%
                                 subsetted_gene_activity_assay$barcode_wo_suffix
                       ]

head(colnames(subsetted_gene_activity_assay))
head(colnames(integrated_labeled_subsetted))

umap_original <-
  integrated_labeled_subsetted@reductions$umap

integrated_labeled_subsetted <-
  integrated_labeled_subsetted[,
                               match(subsetted_gene_activity_assay$barcode_wo_suffix,
                                     integrated_labeled_subsetted$barcode_wo_suffix)]
umap_sorted <-
  integrated_labeled_subsetted@reductions$umap

unique(colnames(subsetted_gene_activity_assay))
length(unique(colnames(subsetted_gene_activity_assay))) # 10000
length(str_split(string = unique(colnames(subsetted_gene_activity_assay)),
                 pattern = "-",
                 simplify = T)[, 1]) # 10000

length(unique(colnames(integrated_labeled))) # 533333
length(unique(str_split(string = colnames(integrated_labeled),
                        pattern = "-",
                        simplify = T)[, 1])) # 533333
# project 

# # test_cell_list <-
# #   df_Lexi_cells_470K[1:2]
# 
# test_matrices_list <-
#   vector(mode = "list",
#          length = length(test_gRangesList))
# names(test_matrices_list) <-
#   c("peak1",
#     "peak2")



# for (i in 1:length(test_matrices_list)) {
#   print(paste("i =",
#               i))
#   system.time({
#     test_matrices_list[[i]] <-
#       FeatureMatrix(fragments = fragments_raw,
#                     features = test_gRangesList[[i]],
#                     process_n = 4000,
#                     verbose = T)
#   })
# 
# }
Signac::GRangesToString(peaks_uncombined[[1]])

test_range <-
  peaks_uncombined[[1]]

peaks_500_width <-
  vector(mode = "list",
         length = length(peaks_uncombined))

convert_GRanges_to_500_peaks <-
  function(x) {
    return_df <-
      data.frame(chr = x@seqnames,
                 start = (x@ranges@start +
                   x@elementMetadata$relative_summit_position -
                   250),
                 end = (x@ranges@start +
                          x@elementMetadata$relative_summit_position +
                          250),
                 peak_ID = x@elementMetadata@listData$name,
                 score = x@elementMetadata@listData$score,
                 strand = x@strand@values,
                 stringsAsFactors = F)
    return(return_df)
  }

for (i in 1:length(peaks_500_width)) {
  print(i)
  peaks_500_width[[i]] <-
    convert_GRanges_to_500_peaks(peaks_uncombined[[i]])
}

save(peaks_500_width,
     file = "df_Lexi_new_peak_set_500bp.RData")

peaks_500_width[[1]]$peak_ID[1]

for (i in 1:length(peaks_500_width)) {
  names(peaks_500_width)[i] <-
    unlist(str_split(string = peaks_500_width[[i]]$peak_ID[1],
                     pattern = "_peak",
                     simplify = T)[, 1])
}

dir.create("Lexi_500bp_new_peaks")

for (i in 1:length(peaks_500_width)) {
  print(i)
  write.table(peaks_500_width[[i]],
              file = paste0("Lexi_500bp_new_peaks",
                            "/",
                            names(peaks_500_width)[i],
                            "_",
                            "500bp.bed"),
              quote = F, sep = "\t",
              row.names = F, col.names = F)
}
