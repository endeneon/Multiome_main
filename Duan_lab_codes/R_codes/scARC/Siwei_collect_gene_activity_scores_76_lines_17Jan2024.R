# Siwei 17 Jan 2024
# Use Lexi's multiomic data of 28 Apr 2023 (new peakset)
# Collect Gene activity scores and plot graphs accordingly


# init ####
library(Seurat)
library(Signac)
# library(sctransform)
# library(glmGamPoi)
# library(EnsDb.Hsapiens.v86)
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

plan("multisession", workers = 8)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# options(Use)

# load data ####
load("multiomic_obj_with_new_peaks_labeled.RData")

# CoveragePlot(object = multiomic_obj_new,
#              region = c("chr1:150065620-150069621"),
#              assay = "ATAC",
#              annotation = T,
#              group.by = "RNA.fine.cell.type",
#              show.bulk = T,
#              peaks = T,
#              # annotation = T,
#              # group.by = "time.ident",
#              # split.by = "RNA.fine.cell.type",
#              sep = c(":", "-"),
#              assay.scale = "separate",
#              region = c("chr1:150067620-150067621"),
#              region.highlight = ranges.show)
# 
# multiomic_obj_new$time.ident

## Collect Gene activity scores #####
## GeneActivity will return a sparse matrix

# counts <-
#   Read10X_h5(filename = "./018_to_029_combined/outs/filtered_feature_bc_matrix.h5")
# 
# 
# 
# multiomic_obj_new_backup <-
#   multiomic_obj_new

# multiomic_obj_new <-
#   multiomic_obj_new_backup

multiomic_obj_new <-
  UpdateSeuratObject(multiomic_obj_new)
# 
DefaultAssay(multiomic_obj_new) <- "ATAC"

fragments_raw <-
  CreateFragmentObject(path = "./018_to_029_combined/outs/atac_fragments.tsv.gz",
                       cells = colnames(multiomic_obj_new)[1:10],
                       validate.fragments = F,
                       verbose = T)

Fragments(multiomic_obj_new) <- NULL
Fragments(multiomic_obj_new) <-
  fragments_raw
# saveRDS(fragments_raw,
#         file = "Lexi_fragments_raw_obj.RDs")


###



fragments_raw <-
  readRDS(file = "Lexi_fragments_raw_obj.RDs")
Fragments(multiomic_obj_new) <-
  fragments_raw

hg38_gene_activity_10k_1k <-
  readRDS("gRanges_hg38_gene_activity_10k_1k_19Jan2024.RDs")
hg38_gene_activity_10k_1k[1]

test_featureMatrix <-
  FeatureMatrix(fragments = fragments_raw,
                features = hg38_gene_activity_10k_1k[1],
                process_n = 2000,
                verbose = T)
df_Lexi_cells_470K <-
  colnames(multiomic_obj_new)
save(df_Lexi_cells_470K,
     file = "df_Lexi_cells_470K.RData")
# multiomic_obj_new <-
#   multiomic_obj_new_backup

# DefaultAssay(multiomic_obj_new) <- "RNA"
# DefaultAssay(multiomic_obj_new) <- "ATAC"
# DefaultAssay(multiomic_obj_new) <- "peaks"

# multiomic_obj_new@assays$peaks@seqinfo <-
#   multiomic_obj_new@assays$peaks@ranges@seqnames


# Error in validObject(result) : invalid class “GRanges” object: 1: 
#   'seqnames(x)' must be a 'factor' Rle
# invalid class “GRanges” object: 2: 
#   'seqlevels(seqinfo(x))' and 'levels(seqnames(x))' are not identical
# 
# multiomic_obj_new@assays[["ATAC"]]@seqinfo@seqnames
# multiomic_obj_new@assays[["ATAC"]]@annotation@seqinfo
# Signac::seqinfo(multiomic_obj_new)
# Signac::seqnames(multiomic_obj_new)
# Signac::seqlengths(multiomic_obj_new)
# 
# Signac::genome(multiomic_obj_new) <-
#   'hg38'
# multiomic_obj_new@assays$ATAC@seqinfo@seqnames <-
#   seqnames(multiomic_obj_new@assays$ATAC@seqinfo@seqnames)
# 
# 
# Signac::seqlevels(Signac::seqinfo(multiomic_obj_new))
# levels(seqnames(multiomic_obj_new))
# seqnames(multiomic_obj_new)
# 
# 
# GenomeInfoDb::seqlevels(GenomeInfoDb::seqinfo(multiomic_obj_new@assays[["ATAC"]]@annotation@seqinfo))
# 
# factor(multiomic_obj_new@assays[["ATAC"]]@annotation@seqinfo)
# 
# temp_seqnames_df <-
#   multiomic_obj_new@assays[["ATAC"]]@annotation@seqinfo
# temp_seqnames_df <-
#   as.data.frame(temp_seqnames_df)
# temp_seqnames_df$chr <-
#   rownames(temp_seqnames_df)
# temp_seqnames_df$chr <-
#   as.numeric(temp_seqnames_df$chr)
# temp_seqnames_df <-
#   temp_seqnames_df[order(temp_seqnames_df$chr), ]
# 
# multiomic_obj_new@assays[["ATAC"]]@seqinfo <-
#   multiomic_obj_new@assays[["ATAC"]]@ranges@seqnames
# 
# seqinfo(multiomic_obj_new)
# multiomic_obj_new@assays[["ATAC"]]@seqinfo@lengths <-
#   temp_seqnames_df$seqlengths
# multiomic_obj_new@assays[["ATAC"]]@seqinfo@is_circular <-
#   temp_seqnames_df$isCircular
# 
# 
# 
# test_type <-
#   test_atac@assays$peaks@annotation@seqinfo@seqnames
# test_atac@assays[["peaks"]]@seqinfo@seqnames
# seqnames(test_atac)
# 
# seqlevels(seqinfo(test_atac))
# levels(seqnames(test_atac))
# length(seqlevels(seqinfo(test_atac)))
# seqinfo(multiomic_obj_new)
# 
# seqnames(multiomic_obj_new) 
# seqlevels(multiomic_obj_new)

# df_gene_activity_raw <-
#   GeneActivity(object = multiomic_obj_new,
#                assay = "ATAC",
#                extend.upstream = 10000,
#                extend.downstream = 1000,,
#                process_n = 20000,
#                verbose = T)

# Lexi_new_peak_set <-
#   multiomic_obj_new@assays$ATAC@data@Dimnames[[1]]
# Lexi_new_peak_set <-
#   matrix(Lexi_new_peak_set,
#          ncol = 1)
# 
# df_Lexi_new_peak_set <-
#   data.frame(CHR = str_split(Lexi_new_peak_set,
#                              pattern = "-",
#                              simplify = T)[, 1],
#              START = str_split(Lexi_new_peak_set,
#                                pattern = "-",
#                                simplify = T)[, 2],
#              END = str_split(Lexi_new_peak_set,
#                                pattern = "-",
#                                simplify = T)[, 3])
# 
# saveRDS(df_Lexi_new_peak_set,
#         file = "df_Lexi_new_peak_set.RData")
# 
# test_chr_seqinfo <- 
#   Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chrM"),
#           seqlengths=c(100, 200, NA, 15),
#           isCircular=c(NA, FALSE, FALSE, TRUE),
#           genome="sasquatch")
# 
# 
# 
# test_atac <-
#   atac_small
# 
# fragments_raw <-
#   CreateFragmentObject(path = "./018_to_029_combined/outs/atac_fragments.tsv.gz",
#                        cells = colnames(multiomic_obj_new),
#                        validate.fragments = F,
#                        verbose = T)
# Fragments(test_atac) <-
#   fragments_raw
# genome(test_atac) <-
#   


# 
test_annotation <-
  Annotation(multiomic_obj_new[["ATAC"]])
seqlevelsStyle(test_annotation) <- 'UCSC'
genome(test_annotation) <- "hg38"
# test_annotation$gene_id
# test_transcripts <-
#   makeGr

## custom functions #####


SetIfNull <-
  function(x, y) {
    if (is.null(x = x)) {
      return(y)
    }
    else {
      return(x)
    }
  }


CollapseToLongestTranscript_2 <-
  function(ranges) {
    range.df <- data.table::as.data.table(x = ranges)
    range.df$strand <- as.character(x = range.df$strand)
    range.df$strand <- ifelse(test = range.df$strand == "*", 
                              yes = "+", no = range.df$strand)
    collapsed <- range.df[, .(unique(seqnames), min(start), 
                              max(end), strand[[1]], gene_biotype[[1]], gene_name[[1]]), 
                          "gene_id"]
    colnames(x = collapsed) <- c("gene_id", "seqnames", "start", 
                                 "end", "strand", "gene_biotype", "gene_name")
    collapsed$gene_name <- .Internal(make.unique(names = collapsed$gene_name,
                                                 sep = "."))
    gene.ranges <- GenomicRanges::makeGRangesFromDataFrame(df = collapsed, 
                                            keep.extra.columns = TRUE)
    return(gene.ranges)
  }


GeneActivity_2 <-
  function(object, assay = NULL, features = NULL, extend.upstream = 2000, 
            extend.downstream = 0, biotypes = "protein_coding", max.width = 5e+05, 
            process_n = 2000, gene.id = FALSE, verbose = TRUE) {
    print("==0==")
    if (!is.null(x = features)) {
      if (length(x = features) == 0) {
        stop("Empty list of features provided")
      }
    }
    print("==1==")
    assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
    if (!inherits(x = object[[assay]], what = "ChromatinAssay")) {
      stop("The requested assay is not a ChromatinAssay.")
    }
    print(object)
    print(assay)
    print("==2==")
    # annotation <- Signac::Annotation(object = object[[assay]])
    annotation <- test_annotation
    print("==2.5==")
    # annotation$gene_name <- ifelse(test = is.na(x = annotation$gene_name) | 
    #                                  (annotation$gene_name == ""), yes = annotation$gene_id, 
    #                                no = annotation$gene_name)
    if (length(x = annotation) == 0) {
      stop("No gene annotations present in object")
    }
    print("==3==")
    if (verbose) {
      message("Extracting gene coordinates")
    }
    transcripts <- CollapseToLongestTranscript_2(ranges = annotation)
    if (gene.id) {
      transcripts$gene_name <- transcripts$gene_id
    }
    if (!is.null(x = biotypes)) {
      transcripts <- transcripts[transcripts$gene_biotype %in% 
                                   biotypes]
      if (length(x = transcripts) == 0) {
        stop("No genes remaining after filtering for requested biotypes")
      }
    }
    if (!is.null(x = features)) {
      transcripts <- transcripts[transcripts$gene_name %in% 
                                   features]
      if (length(x = transcripts) == 0) {
        stop("None of the requested genes were found in the gene annotation")
      }
    }
    if (!is.null(x = max.width)) {
      transcript.keep <- which(x = width(x = transcripts) < 
                                 max.width)
      transcripts <- transcripts[transcript.keep]
      if (length(x = transcripts) == 0) {
        stop("No genes remaining after filtering for max.width")
      }
    }
    transcripts <- Signac::Extend(x = transcripts, upstream = extend.upstream, 
                          downstream = extend.downstream)
    frags <- Signac::Fragments(object = object[[assay]])
    if (length(x = frags) == 0) {
      stop("No fragment information found for requested assay")
    }
    # saveRDS(transcripts,
    #         file = "gRanges_hg38_gene_activity_10k_1k_19Jan2024.RDs")
    cells <- colnames(x = object[[assay]])
    counts <- FeatureMatrix(fragments = frags, features = transcripts,
                            process_n = process_n, cells = cells, verbose = verbose)
    gene.key <- transcripts$gene_name
    names(x = gene.key) <- GenomicRanges::GRangesToString(grange = transcripts)
    rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
    counts <- counts[rownames(x = counts) != "", ]
    return(1)
  }

## Count ####

df_gene_activity_raw <-
  GeneActivity_2(object = multiomic_obj_new,
                 assay = "ATAC",
                 extend.upstream = 10000,
                 extend.downstream = 1000,
                 process_n = 5000,
                 verbose = T)
saveRDS(df_gene_activity_raw,
        file = "df_gene_activity_raw_multiomic_obj_new_18Jan2024.RDs")
