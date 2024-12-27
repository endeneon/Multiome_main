# Siwei 17 Jan 2024
# Use Lexi's multiomic data of 28 Apr 2023 (new peakset)
# and replot ATACseq data w new peak set projected by RNAseq cell types


# init ####
{
  library(Seurat)
  library(Signac)
  library(sctransform)
  library(glmGamPoi)
  # library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
  # library(GenomicFeatures)
  library(AnnotationDbi)
  # library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  
  library(patchwork)
  library(readr)
  library(ggplot2)
  library(RColorBrewer)
  library(dplyr)
  library(viridis)
  library(graphics)
  library(stringr)
  library(future)
  
}

plan("multisession", workers = 4)
set.seed(2001)
options(future.globals.maxSize = 229496729600)

# load data ####
load("multiomic_obj_with_new_peaks_labeled.RData")

# RNA.cell.type
unique(multiomic_obj_new$RNA.fine.cell.type)

multiomic_obj_new$RNA.cell.type.as.100 <-
  multiomic_obj_new$RNA.fine.cell.type

multiomic_obj_new$RNA.cell.type.as.100[multiomic_obj_new$RNA.cell.type.as.100 %in% "unknown"] <-
  'unidentified'
multiomic_obj_new$RNA.cell.type.as.100[multiomic_obj_new$RNA.cell.type.as.100 %in% '?glut'] <-
  'glut? (ambiguous)'
multiomic_obj_new$RNA.cell.type.as.100[multiomic_obj_new$RNA.cell.type.as.100 %in% 'NPC'] <-
  'unidentified'
multiomic_obj_new$RNA.cell.type.as.100[multiomic_obj_new$RNA.cell.type.as.100 %in% 'NEFM+ glut'] <-
  'NEFM+ glut (npglut)'
multiomic_obj_new$RNA.cell.type.as.100[multiomic_obj_new$RNA.cell.type.as.100 %in% 'NEFM- glut'] <-
  'NEFM- glut (nmglut)'

unique(multiomic_obj_new$RNA.cell.type.as.100)
multiomic_obj_new$RNA.cell.type.as.100 <-
  factor(multiomic_obj_new$RNA.cell.type.as.100,
         levels = c("GABA",'NEFM+ glut (npglut)',
                    'NEFM- glut (nmglut)',
                    'glut? (ambiguous)',
                    "unidentified"))

# output resolution 1000x800
DefaultAssay(multiomic_obj_new) <- "ATAC"
DimPlot(object = multiomic_obj_new, 
        label = TRUE, 
        repel = T,
        alpha = 0.5,
        cols = brewer.pal(n = 5,
                          name = "Dark2"),
        # cols = c("#B33E52",#GABA
        #          "#E6D2B8", #?glut
        #          "#CCAA7A", #nmglut
        #          "#f29116", #npglut
        #          "#347545"), #unknown
        group.by = "RNA.cell.type.as.100") + 
  ggtitle("ATACseq data projected\n by RNAseq cell types") +
  theme(text = element_text(size = 10)) 

# remove cluster 7, 8, 10, 11, 16, 18
multiomic_obj_new_plot <-
  subset(multiomic_obj_new,
         subset = (seurat_clusters %in% c(7, 8, 10, 11, 16, 18)),
         invert = T)

# multiomic_obj_new_plot_umap <-
#   RunPCA(multiomic_obj_new_plot,
#          npcs = 50,
#          seed.use = 42)

multiomic_obj_new_plot_umap <-
  RunUMAP(multiomic_obj_new_plot,
          reduction = "lsi",
          reduction.name = "umap_v2",
          seed.use = 42,
          dims = 1:50,
          # graph = "ATAC_snn",
          verbose = T)

DimPlot(object = multiomic_obj_new_plot_umap, 
        reduction = "umap_v2",
        label = TRUE, 
        repel = T,
        alpha = 0.5,
        
        # cols = c("#B33E52",#GABA
        #          "#E6D2B8", #?glut
        #          "#CCAA7A", #nmglut
        #          "#f29116", #npglut
        #          "#347545"), #unknown
        group.by = "seurat_clusters") + 
  ggtitle("ATACseq data projected\n by RNAseq cell types") +
  theme(text = element_text(size = 10)) 

multiomic_obj_new_plot_umap$RNA.cell.type.as.100[multiomic_obj_new_plot_umap$seurat_clusters %in% c(15)] <- "unidentified"
multiomic_obj_new_plot_umap$RNA.cell.type.as.100[multiomic_obj_new_plot_umap$seurat_clusters %in% c(14, 17)] <- 'glut? (ambiguous)'



DimPlot(object = multiomic_obj_new_plot_umap, 
        reduction = "umap_v2",
        label = F, 
        # repel = T,
        alpha = 0.5,
        cols = brewer.pal(n = 5,
                          name = "Dark2"),
        # cols = c("#B33E52",#GABA
        #          "#E6D2B8", #?glut
        #          "#CCAA7A", #nmglut
        #          "#f29116", #npglut
        #          "#347545"), #unknown
        group.by = "RNA.cell.type.as.100") + 
  ggtitle("ATACseq data projected\n by RNAseq cell types") +
  theme(text = element_text(size = 10)) 

# ! use this one for ATAC projection !
saveRDS(object = multiomic_obj_new_plot_umap,
        file = "multiomic_obj_78_ATAC_new_plot_umap_20Aug2024.RDs")

multiomic_obj_new_plot_umap <-
  readRDS("multiomic_obj_78_ATAC_new_plot_umap_20Aug2024.RDs")

DimPlot(object = multiomic_obj_new_plot_umap, 
        reduction = "umap_v2",
        label = F, 
        # repel = T,
        alpha = 0.5,
        # cols = c("#B33E52",#GABA
        #          "#E6D2B8", #?glut
        #          "#CCAA7A", #nmglut
        #          "#f29116", #npglut
        #          "#347545"), #unknown
        group.by = "seurat_clusters") + 
  ggtitle("ATACseq data projected\n by RNAseq cell types") +
  theme(text = element_text(size = 10)) 

FeaturePlot(multiomic_obj_new_plot_umap,
            assay)


DimPlot(object = multiomic_obj_new, 
        label = F, 
        repel = T,
        alpha = 0.5,
        cols = c("#B33E52",#GABA
                 "#E6D2B8", #?glut
                 "#CCAA7A", #nmglut
                 "#f29116", #npglut
                 "#347545"), #unknown
        group.by = "RNA.cell.type.as.100") + 
  ggtitle("ATACseq data projected\n by RNAseq cell types") +
  theme(text = element_text(size = 10)) 

multiomic_obj_new_plot_umap <-
  TSSEnrichment(multiomic_obj_new_plot_umap,
                fast = T,
                assay = "ATAC",
                process_n = 2000)
unique(multiomic_obj_new_plot_umap$orig.ident)
TSSPlot(multiomic_obj_new_plot_umap,
        group.by = "orig.ident")

multiomic_obj_new_plot_umap[['GACT']] <-
  CreateAssayObject(counts = GeneActivity(multiomic_obj_new_plot_umap))




# save(df_gene_activity_raw,
#      file = "df_gene_activity_raw_multiomic_obj__76_lines_26Jan2024.RData")

# rm(df_gene_activity_raw)

# gene_activity_names <-
#   Get_gene_names(object = multiomic_obj,
#                  assay = "ATAC",
#                  extend.upstream = 10000,
#                  extend.downstream = 1000,
#                  process_n = 2000,
#                  verbose = T)


## RNA #####
load("../Integrated_all_RNAseq_data_analysis/018-030_RNA_integrated_labeled_with_harmony.RData")
DefaultAssay(integrated_labeled) <- "integrated"

unique(integrated_labeled$cell.type.forplot)
integrated_labeled$RNA.cell.type.as.100 <-
  factor(integrated_labeled$cell.type.forplot,
         levels = c("GABA",
                    'glut?',
                    'NEFM- glut',
                    'NEFM+ glut',
                    "unidentified"))

DimPlot(object = integrated_labeled, 
        label = TRUE, 
        repel = T,
        alpha = 0.5,
        cols = c("#B33E52",#GABA
                 "#E6D2B8", #?glut
                 "#CCAA7A", #nmglut
                 "#f29116", #npglut
                 "#347545"), #unknown
        group.by = "RNA.cell.type.as.100") + 
  ggtitle("Cells labelled by type") +
  theme(text = element_text(size = 10)) 
DimPlot(object = integrated_labeled, 
        label = F, 
        repel = T,
        alpha = 0.5,
        cols = c("#B33E52",#GABA
                 "#E6D2B8", #?glut
                 "#CCAA7A", #nmglut
                 "#f29116", #npglut
                 "#347545"), #unknown
        group.by = "RNA.cell.type.as.100") + 
  ggtitle("Cells labelled by type") +
  theme(text = element_text(size = 10)) 





test_annotation <-
  Signac::Annotation(multiomic_obj_new_plot_umap[["ATAC"]])
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

## Collect Gene activity scores #####
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
    counts <- FeatureMatrix(fragments = frags, features = transcripts,
                            process_n = process_n, cells = cells, verbose = verbose)
    # save(counts,
    #      file = "df_gene_activity_counts_26Jan2024.RData")
    gene.key <- transcripts$gene_name
    names(x = gene.key) <- Signac::GRangesToString(grange = transcripts)
    key_2_return <- 
      list(gene.key,
           counts)
    # as.vector(x = gene.key[rownames(x = counts)])
    # rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
    # counts <- counts[rownames(x = counts) != "", ]
    return(key_2_return)
  }

test_annotation <-
  Signac::Annotation(multiomic_obj_new_plot_umap[["ATAC"]])
gene_activity_names <-
  Get_gene_names(object = multiomic_obj_new_plot_umap,
                 assay = "ATAC",
                 extend.upstream = 10000,
                 extend.downstream = 1000,
                 process_n = 2000,
                 verbose = T)
gene_activity_names[1]

df_raw_counts <-
  FeatureMatrix(fragments = multiomic_obj_new_plot_umap@assays$ATAC@fragments,
                features = names(gene_activity_names),
                verbose = T)
saveRDS(df_raw_counts,
        file = "df_gact_76_lines_raw_counts.RDs")
rownames(df_raw_counts)
counts_filtered <-
  df_raw_counts[rownames(df_raw_counts) %in% names(gene_activity_names), ]

gene_activity_76lines <-
  CreateAssayObject(counts = counts_filtered)
# gene_activity_76lines <-
#   CreateSeuratObject(counts = gene_activity_76lines,
#                      assay = 'gact')
# gene_activity_76lines <-
#   NormalizeData(object = gene_activity_76lines,
#                 normalization.method = "LogNormalize",
#                 scale.factor = median(gene_activity_76lines$nCount_gact))
# 
# DefaultAssay(multiomic_obj_new_plot_umap) <- "RNA"
# multiomic_obj_new_plot_umap[["gact"]] <-
#   gene_activity_76lines
save.image(file = "multiomic_obj_new_plot_umap_w_gact.RData")
multiomic_obj_new_plot_umap <-
  NormalizeData(object = multiomic_obj_new_plot_umap,
                assay = "gact",
                normalization.method = "LogNormalize",
                scale.factor = median(multiomic_obj_new_plot_umap$nCount_gact))

DefaultAssay(multiomic_obj_new_plot_umap) <- "gact"
FeaturePlot(multiomic_obj_new_plot_umap,
            reduction = "umap_v2",
            features = c(names(gene_activity_names)[gene_activity_names == "SLC17A6"],
                         names(gene_activity_names)[gene_activity_names == "GAD1"],
                         names(gene_activity_names)[gene_activity_names == "FOS"]),
            max.cutoff = 2,
            ncol = 3,
            alpha = 0.8)
# DefaultAssay(multiomic_obj_new_plot_umap) <- "RNA"

DefaultAssay(multiomic_obj_new_plot_umap) <- "ATAC"
multiomic_obj_new_plot_umap <-
  TSSEnrichment(multiomic_obj_new_plot_umap)

multiomic_obj_new_plot_umap <-
  readRDS("multiomic_obj_w_gact_TSS_NS.RDs")
DefaultAssay(multiomic_obj_new_plot_umap) <- "ATAC"
TSSPlot(multiomic_obj_new_plot_umap)
# NucleosomeSignal(multiomic_obj_new_plot_umap)
Idents(multiomic_obj_new_plot_umap)
unique(multiomic_obj_new_plot_umap$group.ident)
unique(multiomic_obj_new_plot_umap$orig.ident)

Idents(multiomic_obj_new_plot_umap) <- "group.ident"
FragmentHistogram(multiomic_obj_new_plot_umap,
                  group.by = "group.ident")
multiomic_obj_new_plot_umap$plot_groups <-
  str_c('Library',
        multiomic_obj_new_plot_umap$group.ident,
        "hr",
        sep = ' ')
multiomic_obj_new_plot_umap$subset_groups <-
  str_split(multiomic_obj_new_plot_umap$group.ident,
            pattern = '-',
            simplify = T)[, 1]
unique(multiomic_obj_new_plot_umap$subset_groups)

FragmentHistogram(multiomic_obj_new_plot_umap,
                  group.by = "plot_groups")

Idents(multiomic_obj_new_plot_umap) <- "subset_groups"
multiomic_obj_new_plot_umap_good <-
  subset(multiomic_obj_new_plot_umap,
         idents = c("22",
                    "36",
                    "39",
                    "44",
                    "49"),
         invert = T)

Idents(multiomic_obj_new_plot_umap_good) <- "plot_groups"
FragmentHistogram(multiomic_obj_new_plot_umap_good,
                  group.by = "plot_groups")

# TSSPlot ####
multiomic_obj_new_plot_umap <-
  readRDS("multiomic_obj_w_gact_TSS.RData")
DefaultAssay(multiomic_obj_new_plot_umap) <- "ATAC"
TSSPlot(multiomic_obj_new_plot_umap)
# NucleosomeSignal(multiomic_obj_new_plot_umap)
Idents(multiomic_obj_new_plot_umap)
unique(multiomic_obj_new_plot_umap$group.ident)
unique(multiomic_obj_new_plot_umap$orig.ident)

Idents(multiomic_obj_new_plot_umap) <- "group.ident"
FragmentHistogram(multiomic_obj_new_plot_umap,
                  group.by = "group.ident")
multiomic_obj_new_plot_umap$plot_groups <-
  str_c('Library',
        multiomic_obj_new_plot_umap$group.ident,
        "hr",
        sep = ' ')
multiomic_obj_new_plot_umap$subset_groups <-
  str_split(multiomic_obj_new_plot_umap$group.ident,
            pattern = '-',
            simplify = T)[, 1]
unique(multiomic_obj_new_plot_umap$subset_groups)

FragmentHistogram(multiomic_obj_new_plot_umap,
                  group.by = "plot_groups")

Idents(multiomic_obj_new_plot_umap) <- "subset_groups"
multiomic_obj_new_plot_umap_good <-
  subset(multiomic_obj_new_plot_umap,
         idents = c("22",
                    "36",
                    "39",
                    "44",
                    "49"),
         invert = T)

Idents(multiomic_obj_new_plot_umap_good) <- "plot_groups"
FragmentHistogram(multiomic_obj_new_plot_umap_good,
                  group.by = "plot_groups")
