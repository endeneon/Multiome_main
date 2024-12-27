# Siwei 18 Jan 2024
# try to run DA analysis and link Lexi's 76 line RNASeq data to ATACSeq data

# init ####
library(Seurat)
library(Signac)

library(harmony)

library(future)

# library(EnsDb.Hsapiens.v86)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)

library(stringr)
# library(plyranges)

plan("multisession", workers = 8)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# load data ####
# npglut 0hr = 8
# load("../018-029_combined_analysis/multiomic_obj_with_new_peaks_labeled.RData")

load("~/NVME/scARC_Duan_018/018-029_combined_analysis/018-029_macs2_called_new_peaks_sep_by_timextype.RData")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/018-029_multiomic_obj_clustered_added_metadata.RData")

# need to subset cells: npglut-0hr
DefaultAssay(multiomic_obj) <- "ATAC"
unique(multiomic_obj$timextype.ident)
table(multiomic_obj$timextype.ident)

npglut_0hr_cells <-
  colnames(multiomic_obj)[multiomic_obj$timextype.ident %in% "npglut_0hr"]

npglut_0hr_peaks <-
  peaks_uncombined[[8]]
# confirm peak identity
unique(npglut_0hr_peaks$ident)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
npglut_0hr_peaks <-
  keepStandardChromosomes(npglut_0hr_peaks,
                          pruning.mode = "coarse")
npglut_0hr_peaks <-
  subsetByOverlaps(x = npglut_0hr_peaks,
                   ranges = blacklist_hg38_unified,
                   invert = TRUE)
unique(npglut_0hr_peaks@seqnames)
npglut_0hr_peaks_autosomes <-
  npglut_0hr_peaks %>%
  plyranges::filter(!(seqnames %in% c("chrX",
                                    "chrY")),
                    .preserve = T)
unique(npglut_0hr_peaks_autosomes@seqnames)


fragments_raw <-
  readRDS(file = "../018-029_combined_analysis/Lexi_fragments_raw_obj.RDs")


# collect raw counts ####
counts_raw <-
  FeatureMatrix(fragments = fragments_raw,
                features = npglut_0hr_peaks_autosomes,
                cells = npglut_0hr_cells,
                process_n = 5000,
                verbose = T)
# saveRDS(counts_raw,
#         file = "npglut_0hr_specific_peaks_autosome_featurematrix.RDs")

load("../R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")

assay_npglut_0hr <-
  CreateChromatinAssay(counts = counts_raw,
                       fragments = fragments_raw,
                       min.cells = 0,
                       min.features = 0,
                       sep = c(":", "-"),
                       # fragments = frag,
                       annotation = ens_use)
seurat_npglut_0hr <-
  CreateSeuratObject(counts = assay_npglut_0hr,
                     assay = "peaks")

master_metadata_table <-
  multiomic_obj@meta.data
npglut_0hr_cells <-
  colnames(seurat_npglut_0hr)

seurat_metadata_table <-
  master_metadata_table[rownames(master_metadata_table) %in% npglut_0hr_cells, ]

seurat_npglut_0hr@meta.data <-
  seurat_metadata_table

## load covar table ####
load("~/NVME/scARC_Duan_018/case_control_DE_analysis/018-029_covar_table_final.RData")
### cut 0hr
covar_table_0hr <-
  covar_table_final[covar_table_final$time %in% "0hr", ]
covar_table_0hr <-
  covar_table_0hr[covar_table_0hr$cell_line %in% unique(seurat_npglut_0hr$cell.line.ident), ]

seurat_npglut_0hr <-
  seurat_npglut_0hr[, seurat_npglut_0hr$cell.line.ident %in% unique(covar_table_0hr$cell_line)]

## Run TFIDF ####
seurat_npglut_0hr <-
  RunTFIDF(seurat_npglut_0hr,
           verbose = T)
seurat_npglut_0hr <-
  FindTopFeatures(seurat_npglut_0hr,
                  min.cutoff = "q5",
                  verbose = T)
seurat_npglut_0hr <-
  RunSVD(seurat_npglut_0hr,
         verbose = T)
seurat_npglut_0hr <-
  RunUMAP(seurat_npglut_0hr,
          reduction = "lsi",
          dims = 2:30)

## DA Analysis ####
head(seurat_npglut_0hr)
## assemble identity table
covar_table_0hr_short <-
  covar_table_0hr[ , c("cell_line",
                       "aff")]
seurat_indiv_table <-
  data.frame(barcode = colnames(seurat_npglut_0hr),
             cell_line = seurat_npglut_0hr$cell.line.ident)
seurat_indiv_table <-
  merge(x = seurat_indiv_table,
        y = covar_table_0hr_short,
        by = "cell_line",
        all.x = T)
sum(is.na(seurat_indiv_table$aff))
## assign aff identity to all cells
head(colnames(seurat_npglut_0hr))
head(seurat_indiv_table$barcode)

seurat_indiv_table <-
  seurat_indiv_table[match(x = colnames(seurat_npglut_0hr),
                           table = seurat_indiv_table$barcode), ]
head(seurat_indiv_table$barcode)
seurat_npglut_0hr$aff <-
  seurat_indiv_table$aff
unique(seurat_npglut_0hr$aff)

Idents(seurat_npglut_0hr) <- "aff"

da_seurat_peaks <-
  FindMarkers(object = seurat_npglut_0hr,
              ident.1 = "case",
              ident.2 = "control", # if avg_logFC > 0 then ident.1 > ident.2
              test.use = 'LR',
              latent.vars = 'nCount_peaks',
              base = 2,
              random.seed = 42,
              verbose = T)
save.image("DA_seurat_peaks_19Jan2023_w_results.RData")

df_DA_peaks_4_gRange <-
  rownames(da_seurat_peaks)
df_DA_peaks_4_gRange <-
  data.frame(seqnames = str_split(df_DA_peaks_4_gRange,
                                  pattern = "-",
                                  simplify = T)[, 1],
             start = str_split(df_DA_peaks_4_gRange,
                               pattern = "-",
                               simplify = T)[, 2],
             end = str_split(df_DA_peaks_4_gRange,
                             pattern = "-",
                             simplify = T)[, 3])
gRanges_DA_peaks <-
  makeGRangesFromDataFrame(df_DA_peaks_4_gRange,
                           keep.extra.columns = F,
                           ignore.strand = T,
                           seqnames.field = "seqnames",
                           start.field = "start",
                           end.field = "end",
                           na.rm = T)

ens_annotation <-
  GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(ens_annotation) <- 'UCSC' # requires ftp access
seqlevels(ens_annotation) <-
  paste0('chr',
         seqlevels(ens_annotation))
genome(ens_annotation) <-
  "hg38"
ens_annotation <-
  sortSeqlevels(ens_annotation)
# txdb_hg38 <-
#   TxDb.Hsapiens.UCSC.hg38.knownGene
# columns(TxDb.Hsapiens.UCSC.hg38.knownGene)
# TxDb.Hsapiens.UCSC.hg38.knownGene
# keys(TxDb.Hsapiens.UCSC.hg38.knownGene)
# keytypes(TxDb.Hsapiens.UCSC.hg38.knownGene)
# head(keys(TxDb.Hsapiens.UCSC.hg38.knownGene,
#           "GENEID"))
#
# AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg38.knownGene,
#                       keys = head(keys(TxDb.Hsapiens.UCSC.hg38.knownGene)),
#                       keytype = "GENEID",
#                       columns = c("TXID",
#                                   "TXNAME"))
#
# AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg38.knownGene,
#                       keys = keys(TxDb.Hsapiens.UCSC.hg38.knownGene),
#                       keytype = "GENEID",
#                       columns = c("CDSCHROM",
#                                   "TXNAME"))
#
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
keys(org.Hs.eg.db)

# head(keys(org.Hs.eg.db))
# AnnotationDbi::select(org.Hs.eg.db,
#        keys = head(keys(org.Hs.eg.db)),
#        columns = c("ALIAS",
#                    "SYMBOL"))

seurat_npglut_0hr_backup <-
  seurat_npglut_0hr

genome(seurat_npglut_0hr) <-
  'hg38'
Annotation(seurat_npglut_0hr)
seqlevels(seurat_npglut_0hr)
seqlevelsStyle(seurat_npglut_0hr)

genome(seurat_npglut_0hr)

genome(gRanges_DA_peaks) <- 'hg38'
genome(gRanges_DA_peaks)

da_seurat_peaks$closest_genes <-
  ClosestFeature(seurat_npglut_0hr,
                 regions = gRanges_DA_peaks,
                 annotation = ens_annotation)
da_seurat_peaks$peak_coord <-
  da_seurat_peaks$closest_genes$query_region

peak_ID_meta <-
  data.frame(peak_coord = GRangesToString(npglut_0hr_peaks),
             peak_id = npglut_0hr_peaks@elementMetadata@listData$name)
da_seurat_peaks_id <-
  merge(x = da_seurat_peaks,
        y = peak_ID_meta,
        by = "peak_coord")
da_seurat_peaks_id <-
  da_seurat_peaks_id[order(da_seurat_peaks_id$p_val), ]

write.table(da_seurat_peaks,
            file = "npglut_0hr_aff_vs_ctrl_DAG_annotated.tsv",
            quote = F, sep = "\t",
            row.names = F, col.names = T)

write.table(da_seurat_peaks_id,
            file = "npglut_0hr_aff_vs_ctrl_DAG_annotated_w_peakid.tsv",
            quote = F, sep = "\t",
            row.names = F, col.names = T)

save.image("DA_seurat_peaks_19Jan2023_w_results.RData")
# run RNA-seq data normalisation with Harmony
# DefaultAssay(multiomic_obj_new) <- "RNA"
# Idents(multiomic_obj_new) <- "orig.ident"
#
# multiomic_obj_new <-
#   UpdateSeuratObject(multiomic_obj_new)
#
# multiomic_obj_new[["percent.mt"]] <-
#   PercentageFeatureSet(object = multiomic_obj_new,
#                        assay = "RNA",
#                        pattern = "^MT-")
#
# multiomic_obj_new <-
#   SCTransform(multiomic_obj_new,
#               assay = "RNA",
#               new.assay.name = "SCT",
#               vars.to.regress = "percent.mt",
#               verbose = T,
#               seed.use = 42)
#
# multiomic_obj_new <-
#   IntegrateLayers(object = multiomic_obj_new,
#                   assay = "SCT",
#                   method = HarmonyIntegration,
#                   orig.reduction = "pca",
#                   verbose = T)

## Run Harmony #####
Idents(seurat_npglut_0hr) <- "group.ident"

data_matrix_raw <-
  as.matrix(seurat_npglut_0hr@assays$peaks@counts)

df_data_matrix_raw_4_gRange <-
  rownames(data_matrix_raw)
df_data_matrix_raw_4_gRange <-
  data.frame(seqnames = str_split(df_data_matrix_raw_4_gRange,
                                  pattern = "-",
                                  simplify = T)[, 1],
             start = str_split(df_data_matrix_raw_4_gRange,
                               pattern = "-",
                               simplify = T)[, 2],
             end = str_split(df_data_matrix_raw_4_gRange,
                             pattern = "-",
                             simplify = T)[, 3])
gRanges_data_matrix_raw <-
  makeGRangesFromDataFrame(df_data_matrix_raw_4_gRange,
                           keep.extra.columns = F,
                           ignore.strand = T,
                           seqnames.field = "seqnames",
                           start.field = "start",
                           end.field = "end",
                           na.rm = T)



seurat_npglut_0hr_Harmony <-
  CreateChromatinAssay(counts = data_matrix_raw,
                       ranges = gRanges_data_matrix_raw,
                       # sep = c("-", "-"),
                       genome = "hg38",
                       annotation = ens_annotation)
seurat_npglut_0hr_Harmony <-
  CreateSeuratObject(counts = seurat_npglut_0hr_Harmony,
                     assay = "peaks",
                     meta.data = seurat_npglut_0hr@meta.data)
DefaultAssay(seurat_npglut_0hr_Harmony)
Idents(seurat_npglut_0hr_Harmony) <- "group.ident"

# seurat_npglut_0hr_Harmony <-
#   seurat_npglut_0hr
# seurat_npglut_0hr_Harmony <-
#   UpdateSeuratObject(seurat_npglut_0hr_Harmony)
# # # DefaultAssay(seurat_npglut_0hr_Harmony)
# Idents(seurat_npglut_0hr_Harmony) <- "group.ident"
# seurat_npglut_0hr_Harmony@assays$peaks@scale.data <-
#   as.matrix(as.data.frame(seurat_npglut_0hr_Harmony@assays$peaks@data))
seurat_npglut_0hr_Harmony <-
  NormalizeData(seurat_npglut_0hr_Harmony,
                normalization.method = "LogNormalize")
seurat_npglut_0hr_Harmony <-
  FindVariableFeatures(seurat_npglut_0hr_Harmony,
                       selection.method = "vst")
# save.image()
save.image("DA_seurat_peaks_19Jan2023_w_results_harmony.RData")

## restart from here
seurat_npglut_0hr_Harmony <-
  ScaleData(seurat_npglut_0hr)

# rm(data_matrix_raw)

seurat_npglut_0hr_Harmony <-
  RunPCA(seurat_npglut_0hr_Harmony,
         seed.use = 42,
         verbose = T)
# > colnames(seurat_npglut_0hr@meta.data)
# [1] "orig.ident"            "nCount_RNA"            "nFeature_RNA"
# [4] "nCount_ATAC"           "nFeature_ATAC"         "group.ident"
# [7] "cell.line.ident"       "nucleosome_signal"     "nucleosome_percentile"
# [10] "TSS.enrichment"        "TSS.percentile"        "ATAC_snn_res.0.5"
# [13] "seurat_clusters"       "RNA.cell.type"         "RNA.fine.cell.type"
# [16] "time.ident"            "timextype.ident"       "nCount_peaks"
# [19] "nFeature_peaks"        "aff"

seurat_npglut_0hr_Harmony <-
  RunHarmony(seurat_npglut_0hr_Harmony,
             group.by.vars = c("group.ident",
                               "TSS.percentile",
                               "nucleosome_percentile"))

save.image("DA_seurat_peaks_19Jan2023_w_results_harmony.RData")
