# Siwei 22 Jan 2024
# Use 18 lines, add GeneActivity assay

# init ####
library(Seurat)
library(Signac)

# library(harmony)

library(future)

library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)

library(stringr)
# library(plyranges)

plan("multisession", workers = 16)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# load data #####
# load("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/multiomic_obj.RData")
load("/nvmefs/scARC_Duan_018/018-029_combined_analysis/018-029_multiomioc_obj_added_TSS_nucsig.RData")

# DefaultAssay(integrated_labeled) <- "integrated"
# FeaturePlot(integrated_labeled,
#             features = "SLC17A6")
# multiomic_obj <-
#   UpdateSeuratObject(multiomic_obj)
DefaultAssay(multiomic_obj) <- "ATAC"

# levels(seqnames(multiomic_obj))
# seqlevels(seqinfo(multiomic_obj))

# multiomic_obj_backup <-
#   multiomic_obj

# seqlevels(ens_use) <-
#   Rle(factor(paste0('chr',
#                     seqlevels(ens_use))))



# DefaultAssay(multiomic_obj) <- "peaks"
# levels(seqnames(multiomic_obj))
# seqlevels(seqinfo(multiomic_obj))

# 
# seqnames(multiomic_obj) <-
#   Rle(factor(multiomic_obj@assays[["ATAC"]]@ranges@seqnames@values))
# 
# seqlevels(seqinfo(multiomic_obj)) <-
#   factor()
# # calculate gene activity  ####
## custom functions #####

test_annotation <-
  Signac::Annotation(multiomic_obj[["ATAC"]])
# seqlevels(test_annotation)
# seqlevels(test_annotation) <-
#   paste0("chr", seqlevels(test_annotation))
# test_annotation <-
#   sortSeqlevels(test_annotation)
genome(test_annotation) <- "hg38"

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
    # annotation$gene_name <-
    #   test_annotation$gene_name
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
    save(counts,
         file = "df_gene_activity_counts_76_lines_26Jan2024.RData")
    gene.key <- transcripts$gene_name
    names(x = gene.key) <- Signac::GRangesToString(grange = transcripts)
    rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
    counts <- counts[rownames(x = counts) != "", ]
    return(1)
  }

Get_gene_names <-
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
    # annotation$gene_name <-
    #   test_annotation$gene_name
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
    # counts <- FeatureMatrix(fragments = frags, features = transcripts,
    #                         process_n = process_n, cells = cells, verbose = verbose)
    # save(counts,
    #      file = "df_gene_activity_counts_26Jan2024.RData")
    gene.key <- transcripts$gene_name
    names(x = gene.key) <- Signac::GRangesToString(grange = transcripts)
    key_2_return <- gene.key
      # as.vector(x = gene.key[rownames(x = counts)])
    # rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
    # counts <- counts[rownames(x = counts) != "", ]
    return(key_2_return)
  }
## Count ####
# 

df_gene_activity_raw <-
  GeneActivity_2(object = multiomic_obj,
                 assay = "ATAC",
                 extend.upstream = 10000,
                 extend.downstream = 1000,
                 process_n = 2000,
                 verbose = T)
save(df_gene_activity_raw,
        file = "df_gene_activity_raw_multiomic_obj__76_lines_26Jan2024.RData")

rm(df_gene_activity_raw)

gene_activity_names <-
  Get_gene_names(object = multiomic_obj,
                 assay = "ATAC",
                 extend.upstream = 10000,
                 extend.downstream = 1000,
                 process_n = 2000,
                 verbose = T)
save(gene_activity_names,
     file = "gene_activity_names_76_lines_26Jan2024.RData")

# load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/df_gene_activity_counts_26Jan2024.RData")
# 
# # rownames(x = counts) <- as.vector(x = gene_activity_names[rownames(x = counts)])
# # counts <- counts[rownames(x = counts) != "", ]
# names(gene_activity_names)
# rownames(counts)
# counts_filtered <-
#   counts[rownames(counts) %in% names(gene_activity_names), ]
# 
# head(rownames(counts_filtered))
# head((names(gene_activity_names)))
# rownames(counts_filtered) <-
#   names(gene_activity_names)
# counts_filtered@Dimnames[[1]] <-
#   names(gene_activity_names)
# 
# gene_activity_18lines <-
#   CreateAssayObject(counts = counts_filtered)
# seurat_gene_activity_18lines <-
#   CreateSeuratObject(counts = gene_activity_18lines,
#                      assay = "gact")
# seurat_gene_activity_18lines <-
#   NormalizeData(object = seurat_gene_activity_18lines,
#                 normalization.method = "LogNormalize",
#                 scale.factor = median(seurat_gene_activity_18lines$nCount_gact))
# save(seurat_gene_activity_18lines,
#      file = "seurat_gene_activity_18_lines_26Jan2024.RData")
# # GeneActivity()
# # gene.activities <- GeneActivity(multiomic_obj)
# # multiomic_obj[['GACT']] <- CreateAssayObject(counts = gene.activities)
# # multiomic_obj <- NormalizeData(
# #   object = multiomic_obj,
# #   assay = 'GACT',
# #   normalization.method = 'LogNormalize'
# # )
# # 
# 
# 
# # FeaturePlot()
# 
# # assign annotations
