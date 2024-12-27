# Siwei 02 Sept 2024
# Use Lexi's multiomic data of 28 Apr 2023 (new peakset)
# and replot ATACseq data w new peak set projected by RNAseq cell types


# init ####
{
  library(Seurat)
  library(Signac)
  library(sctransform)
  library(glmGamPoi)
  library(EnsDb.Hsapiens.v86)
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

setwd("/nvmefs/scARC_Duan_018/018-029_combined_analysis")

# load
load("multiomic_obj_new_plot_umap_w_gact.RData")

load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")

ens_use <-
  ens_use[ens_use$tx_biotype == "protein_coding"]
ens_use <-
  ens_use[ens_use$gene_biotype == "protein_coding"]
ens_use <-
  ens_use[!duplicated(ens_use$gene_name)]

# gene.ranges <-
#   keepStandardChromosomes(ens_use,
# #                           pruning.mode = "coarse")
# tss.ranges <-
#   resize(ens_use,
#          width = 1,
#          fix = "start")
ens_df <-
  as.data.frame(ens_use)

ens_remake <-
  makeGRangesFromDataFrame(df = ens_df,
                           keep.extra.columns = T)
# ens_remake <-
tss.ranges <-
  resize(ens_remake,
         width = 1,
         fix = "start")

# tss.rages <-
  
# gene.ranges <-
#   genes(EnsDb.Hsapiens.v86)
# annot_aggr_signac <-
#   GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86,
#                       standard.chromosomes = T)
# seqlevelsStyle(gene.ranges) <- 'UCSC'




DefaultAssay(multiomic_obj_new_plot_umap) <- "ATAC"
multiomic_obj_new_plot_umap <-
  TSSEnrichment(multiomic_obj_new_plot_umap, 
                tss.positions = tss.ranges,
                fast = F,
                verbose = T)
saveRDS(multiomic_obj_new_plot_umap,
        file = "multiomic_obj_w_gact_TSS.RData")
q(save = "no")

multiomic_obj_new_plot_umap <-
  readRDS("multiomic_obj_w_gact_TSS.RData")

multiomic_obj_new_plot_umap$plot_groups <-
  str_c('Library',
        multiomic_obj_new_plot_umap$group.ident,
        "hr",
        sep = ' ')
multiomic_obj_new_plot_umap$subset_groups <-
  str_split(multiomic_obj_new_plot_umap$group.ident,
            pattern = '-',
            simplify = T)[, 1]

Idents(multiomic_obj_new_plot_umap) <- "subset_groups"
unique(multiomic_obj_new_plot_umap$subset_groups)

multiomic_obj_new_plot_umap_good <-
  subset(multiomic_obj_new_plot_umap,
         idents = c("22",
                    "36",
                    "39",
                    "44",
                    "49"),
         invert = T)
Idents(multiomic_obj_new_plot_umap_good) <- "plot_groups"
TSSEnrichment(multiomic_obj_new_plot_umap_good,
              fast = T,
              process_n = 20)
TSSPlot(multiomic_obj_new_plot_umap_good, 
        group.by = "cell.line.ident")

# load the RData file with nuclear signal #####
multiomic_obj_new_plot_umap_NS <-
  readRDS("multiomic_obj_w_gact_TSS_NS.RDs")
unique(multiomic_obj_new_plot_umap_NS$group.ident)


multiomic_obj_new_plot_umap_NS$plot_groups <-
  str_c('Library',
        multiomic_obj_new_plot_umap_NS$group.ident,
        "hr",
        sep = ' ')
multiomic_obj_new_plot_umap_NS$subset_groups <-
  str_split(multiomic_obj_new_plot_umap_NS$group.ident,
            pattern = '-',
            simplify = T)[, 1]

Idents(multiomic_obj_new_plot_umap_NS) <- "subset_groups"
unique(multiomic_obj_new_plot_umap_NS$subset_groups)

multiomic_obj_new_plot_umap_good <-
  subset(multiomic_obj_new_plot_umap_NS,
         idents = c("22",
                    "36",
                    "39",
                    "44",
                    "49"),
         invert = T)
Idents(multiomic_obj_new_plot_umap_good) <- "plot_groups"
unique(multiomic_obj_new_plot_umap_good$plot_groups)
unique(multiomic_obj_new_plot_umap_good$subset_groups)

unique(multiomic_obj_new_plot_umap_NS)

Idents(multiomic_obj_new_plot_umap_good) <- "subset_groups"

# saveRDS(multiomic_obj_new_plot_umap_good,
#         file = "multiomic_obj_w_gact_25Nov2024.RData")

# load
multiomic_obj_new_plot_umap_good <-
  readRDS("multiomic_obj_w_gact_25Nov2024.RData")

multiomic_obj_new_plot_umap_good <-
  TSSEnrichment(multiomic_obj_new_plot_umap_good,
                fast = T,
                process_n = 40)
TSSPlot(multiomic_obj_new_plot_umap_good)
DefaultAssay(multiomic_obj_new_plot_umap_good)
multiomic_obj_new_plot_umap_good@assays$ATAC@annotation
seqnames(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation)
typeof(seqlevels(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation))
seqlevels(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation) <-
  factor(unique(seqnames(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation)))

typeof(seqlevels(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation))
typeof(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation)

seqlevels(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation)
levels(seqlevels(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation@seqnames))
seqnames(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation)
levels(seqnames(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation))
seqinfo(seqnames(multiomic_obj_new_plot_umap_good@assays$ATAC))
seqinfo(multiomic_obj_new_plot_umap_good@assays$ATAC)
seqlevels(multiomic_obj_new_plot_umap_good@assays$ATAC)
seqnames(multiomic_obj_new_plot_umap_good@assays$ATAC)
typeof(multiomic_obj_new_plot_umap_good@assays$ATAC)
Annotation(multiomic_obj_new_plot_umap_good@assays$ATAC)
Annotation(multiomic_obj_new_plot_umap_good@assays$gact)

multiomic_obj_new_plot_umap_good@assays$ATAC@annotation@metadata
multiomic_obj_new_plot_umap_good@assays$ATAC@annotation@seqinfo@seqnames <-
  unique(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation@seqnames)
typeof(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation@seqinfo@seqnames)

multiomic_obj_new_plot_umap_good@assays$ATAC@annotation@seqinfo@seqnames <-
  unique(seqnames(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation))

unique(str_remove_all(string = seqnames(multiomic_obj_new_plot_umap_good@assays$ATAC@annotation),
                      pattern = "^chr"))

grange_tss_position



multiomic_obj_new_plot_umap_recount <-
  TSSEnrichment(multiomic_obj_new_plot_umap_good,
                tss.positions = tss.ranges,
                assay = "ATAC",
                fast = F,
                process_n = 2000)

Idents(multiomic_obj_new_plot_umap_recount) <- "time.ident"
unique(multiomic_obj_new_plot_umap_recount$subset_groups)



multiomic_obj_new_plot_umap_recount$library.ident <-
  str_c('Group-',
        multiomic_obj_new_plot_umap_recount$orig.ident,
        '-',
        multiomic_obj_new_plot_umap_recount$time.ident)
unique(multiomic_obj_new_plot_umap_recount$library.ident)
unique(multiomic_obj_new_plot_umap_recount$group.ident)

multiomic_obj_new_plot_umap_recount$group.ident.4plot <-
  str_c("Library",
        multiomic_obj_new_plot_umap_recount$group.ident,
        "hr",
        sep = ' ')

TSSPlot(multiomic_obj_new_plot_umap_recount, 
        # idents = "time.ident", 
        group.by = "group.ident.4plot")

saveRDS(multiomic_obj_new_plot_umap_recount,
        file = "multiomic_obj_w_gact_TSS_27Nov2024.RDs")

multiomic_obj_new_plot_umap_recount <-
  readRDS("multiomic_obj_w_gact_TSS_27Nov2024.RDs")
# TSS_value_2_plot <-
#   multiomic_obj_new_plot_umap_recount@meta.data %>%
#   group_by(group.ident.4plot) %>%
#   summarise(TSS_enrich_value = )

# > unique(multiomic_obj_new_plot_umap_recount$subset_groups)
# [1] "60060" "50040" "5"     "13"    "51"    "35"    "20087" 
# "53"    "09"    "70179" "20088" "21"    "33"    "63"    "23"   

multiomic_obj_new_plot_umap_recount_subset <-
  multiomic_obj_new_plot_umap_recount[, multiomic_obj_new_plot_umap_recount$subset_groups %in% c("5", 
                                                                                                 "13",
                                                                                                 "20088",
                                                                                                 "20087",
                                                                                                 "60060",
                                                                                                 "09",
                                                                                                 "51", 
                                                                                                 "33", 
                                                                                                 "35", 
                                                                                                 "63")]

TSSPlot(multiomic_obj_new_plot_umap_recount_subset, 
        # idents = "time.ident", 
        group.by = "group.ident.4plot") +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))

FragmentHistogram(multiomic_obj_new_plot_umap_recount_subset,
                  group.by = "group.ident.4plot") +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))


ATAC_embed <-
  as.data.frame(multiomic_obj_new_plot_umap_recount@reductions$umap_v2@cell.embeddings)


DefaultAssay(multiomic_obj_new_plot_umap_recount) <- "gact"
# DefaultAssay(multiomic_obj_new_plot_umap_recount) <- "ATAC"
FeaturePlot(multiomic_obj_new_plot_umap_recount,
            reduction = "umap_v2",
            features = c("SLC17A6",
                         "GAD1",
                         "FOS"),
            max.cutoff = 'q75',
            keep.scale = "feature",
            # blend.threshold = c(0, 0.75),
            dims = c(1, 2),
            ncol = 3)

colnames(multiomic_obj_new_plot_umap_recount@reductions$umap_v2@cell.embeddings)
length(multiomic_obj_new_plot_umap_recount@reductions$umap_v2@cell.embeddings)


FeaturePlot(multiomic_obj_new_plot_umap_recount,
            reduction = "umap_v2",
            features = c(names(gene_activity_names)[gene_activity_names == "SLC17A6"],
                         names(gene_activity_names)[gene_activity_names == "GAD1"],
                         names(gene_activity_names)[gene_activity_names == "FOS"]),
            max.cutoff = 2,
            ncol = 3,
            alpha = 0.8)

DefaultAssay(multiomic_obj_new_plot_umap_recount) <- "RNA"
Idents(multiomic_obj_new_plot_umap_recount) <- "RNA.cluster.ident"
DimPlot(multiomic_obj_new_plot_umap_recount,
        reduction = "lsi")

# 03Dec2024 #####

multiomic_obj_new_plot_umap_recount <-
  readRDS("multiomic_obj_w_gact_TSS_27Nov2024.RDs")
Idents(multiomic_obj_new_plot_umap_recount) <- "timextype.ident"
DefaultAssay(multiomic_obj_new_plot_umap_recount) <- "ATAC"
DefaultAssay(multiomic_obj_new_plot_umap_recount) <- "RNA"
Features(multiomic_obj_new_plot_umap_recount)


Annotation(multiomic_obj_new_plot_umap_recount)

CoveragePlot(multiomic_obj_new_plot_umap_recount,
             region = makeGRangesFromDataFrame(df = data.frame(chr = "chr7",
                                                               start = 18086825,
                                                               end = 19002416),
                                               ignore.strand = T),
             # region = "HDAC9",
             # features = "HDAC9",
             expression.assay = "RNA",
             idents = unique(multiomic_obj_new_plot_umap_recount$timextype.ident)[1])


CoveragePlot(multiomic_obj_new_plot_umap_recount,
             # region = c("chr7-18086825-19002416"),
             assay = "ATAC",
             # annotation = F,
             region = "HDAC9",
             features = "HDAC9",
             expression.assay = "RNA"
             # idents = unique(multiomic_obj_new_plot_umap_recount$timextype.ident)[1]
             )


Idents(multiomic_obj_new_plot_umap_recount) <-
  "RNA.cell.type.as.100"
DefaultAssay(multiomic_obj_new_plot_umap_recount) <-
  "RNA"

DimPlot(multiomic_obj_new_plot_umap_recount,
        reduction = "umap_v2",
        cols = c("#B33E52", 
                 "darkorange",

                 "#CCAA7A",
                 "#E6D2B8", 
                 "darkgreen"),
        # cols = brewer.pal(n = 8,
        #                   name = "Dark2")[c(2,1,4,3,5)],
        alpha = 0.5,
        pt.size = 1) +
  NoLegend()


Idents(multiomic_obj_new_plot_umap_recount) <-
  "RNA.cell.type.as.100"
DefaultAssay(multiomic_obj_new_plot_umap_recount) <-
  "RNA"

DimPlot(multiomic_obj_new_plot_umap_recount,
        reduction = "umap_v2",
        cols = c("#B33E52", 
                 "darkorange",
                 
                 "#CCAA7A",
                 "#E6D2B8", 
                 "darkgreen"),
        # cols = brewer.pal(n = 8,
        #                   name = "Dark2")[c(2,1,4,3,5)],
        alpha = 0.2,
        pt.size = 1) +
  NoLegend()
