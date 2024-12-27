# Siwei 18 Jan 2024
# try to run DA analysis and link Lexi's 76 line RNASeq data to ATACSeq data

# init ####
library(Seurat)
library(Signac)

library(harmony)

library(future)

library(EnsDb.Hsapiens.v86)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)

library(stringr)
library(plyranges)

plan("multisession", workers = 16)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# load data ####
# npglut 0hr = 8
# load("../018-029_combined_analysis/multiomic_obj_with_new_peaks_labeled.RData")

load("~/NVME/scARC_Duan_018/018-029_combined_analysis/018-029_macs2_called_new_peaks_sep_by_timextype.RData")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/018-029_multiomic_obj_clustered_added_metadata.RData")

# need to subset cells: npglut-6hr
DefaultAssay(multiomic_obj) <- "ATAC"
unique(multiomic_obj$timextype.ident)
table(multiomic_obj$timextype.ident)

npglut_6hr_cells <-
  colnames(multiomic_obj)[multiomic_obj$timextype.ident %in% "npglut_6hr"]

npglut_6hr_peaks <-
  peaks_uncombined[[9]]
# confirm peak identity
unique(npglut_6hr_peaks$ident)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
npglut_6hr_peaks <-
  keepStandardChromosomes(npglut_6hr_peaks,
                          pruning.mode = "coarse")
npglut_6hr_peaks <-
  subsetByOverlaps(x = npglut_6hr_peaks,
                   ranges = blacklist_hg38_unified,
                   invert = TRUE)
unique(npglut_6hr_peaks@seqnames)
npglut_6hr_peaks_autosomes <-
  npglut_6hr_peaks %>%
  plyranges::filter(!(seqnames %in% c("chrX",
                                      "chrY")),
                    .preserve = T)
unique(npglut_6hr_peaks_autosomes@seqnames)


fragments_raw <-
  readRDS(file = "../018-029_combined_analysis/Lexi_fragments_raw_obj.RDs")


# collect raw counts ####
counts_raw <-
  FeatureMatrix(fragments = fragments_raw,
                features = npglut_6hr_peaks_autosomes,
                cells = npglut_6hr_cells,
                process_n = 5000,
                verbose = T)
# saveRDS(counts_raw,
#         file = "npglut_0hr_specific_peaks_autosome_featurematrix.RDs")

load("../R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")

assay_npglut_6hr <-
  CreateChromatinAssay(counts = counts_raw,
                       fragments = fragments_raw,
                       min.cells = 0,
                       min.features = 0,
                       sep = c(":", "-"),
                       # fragments = frag,
                       annotation = ens_use)
seurat_npglut_6hr <-
  CreateSeuratObject(counts = assay_npglut_6hr,
                     assay = "peaks")

master_metadata_table <-
  multiomic_obj@meta.data
npglut_6hr_cells <-
  colnames(seurat_npglut_6hr)

seurat_metadata_table <-
  master_metadata_table[rownames(master_metadata_table) %in% npglut_6hr_cells, ]

seurat_npglut_6hr@meta.data <-
  seurat_metadata_table

## load covar table ####
load("~/NVME/scARC_Duan_018/case_control_DE_analysis/018-029_covar_table_final.RData")
### cut 6hr
covar_table_6hr <-
  covar_table_final[covar_table_final$time %in% "6hr", ]
covar_table_6hr <-
  covar_table_6hr[covar_table_6hr$cell_line %in% unique(seurat_npglut_6hr$cell.line.ident), ]

# seurat_npglut_6hr <-
#   seurat_npglut_6hr[, seurat_npglut_6hr$cell.line.ident %in% unique(seurat_npglut_6hr$cell_line)]

## Run TFIDF ####
seurat_npglut_6hr <-
  RunTFIDF(seurat_npglut_6hr,
           verbose = T)
seurat_npglut_6hr <-
  FindTopFeatures(seurat_npglut_6hr,
                  min.cutoff = "q5",
                  verbose = T)
seurat_npglut_6hr <-
  RunSVD(seurat_npglut_6hr,
         verbose = T)
seurat_npglut_6hr <-
  RunUMAP(seurat_npglut_6hr,
          reduction = "lsi",
          dims = 2:30)

## DA Analysis ####
head(seurat_npglut_6hr)
## assemble identity table
covar_table_6hr_short <-
  covar_table_6hr[ , c("cell_line",
                       "aff")]
seurat_indiv_table <-
  data.frame(barcode = colnames(seurat_npglut_6hr),
             cell_line = seurat_npglut_6hr$cell.line.ident)
seurat_indiv_table <-
  merge(x = seurat_indiv_table,
        y = covar_table_6hr_short,
        by = "cell_line",
        all.x = T)
sum(is.na(seurat_indiv_table$aff))
## assign aff identity to all cells
head(colnames(seurat_npglut_6hr))
head(seurat_indiv_table$barcode)

seurat_indiv_table <-
  seurat_indiv_table[match(x = colnames(seurat_npglut_6hr),
                           table = seurat_indiv_table$barcode), ]
head(seurat_indiv_table$barcode)
seurat_npglut_6hr$aff <-
  seurat_indiv_table$aff
unique(seurat_npglut_6hr$aff)

Idents(seurat_npglut_6hr) <- "aff"

da_seurat_peaks <-
  FindMarkers(object = seurat_npglut_6hr,
              ident.1 = "case",
              ident.2 = "control", # if avg_logFC > 0 then ident.1 > ident.2
              test.use = 'LR',
              latent.vars = 'nCount_peaks',
              base = 2,
              random.seed = 42,
              verbose = T)
save(da_seurat_peaks,
     file = "DA_seurat_peaks_npglut_6hr_25Jan2023_w_results.RData")


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

seurat_npglut_6hr_backup <-
  seurat_npglut_6hr

genome(seurat_npglut_6hr) <-
  'hg38'
genome(gRanges_DA_peaks) <- 'hg38'


da_seurat_peaks$closest_genes <-
  ClosestFeature(seurat_npglut_6hr,
                 regions = gRanges_DA_peaks,
                 annotation = ens_annotation)
da_seurat_peaks$peak_coord <-
  da_seurat_peaks$closest_genes$query_region

peak_ID_meta <-
  data.frame(peak_coord = GRangesToString(npglut_6hr_peaks_autosomes),
             peak_id = npglut_6hr_peaks_autosomes@elementMetadata@listData$name)
da_seurat_peaks_id <-
  merge(x = da_seurat_peaks,
        y = peak_ID_meta,
        by = "peak_coord")
da_seurat_peaks_id <-
  da_seurat_peaks_id[order(da_seurat_peaks_id$p_val), ]

# write.table(da_seurat_peaks,
#             file = "npglut_6hr_aff_vs_ctrl_DAG_annotated.tsv",
#             quote = F, sep = "\t",
#             row.names = F, col.names = T)

write.table(da_seurat_peaks_id,
            file = "npglut_6hr_aff_vs_ctrl_DAG_annotated_w_peakid.tsv",
            quote = F, sep = "\t",
            row.names = F, col.names = T)

rm(multiomic_obj)
rm(counts_raw)
rm(seurat_npglut_6hr_backup)

save.image("DA_seurat_peaks_npglut_6hr_26Jan2023_w_results.RData")
