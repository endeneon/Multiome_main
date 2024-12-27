# Siwei 12 Feb 2024
# Calc chromatin accessibility correlation using Yifan's code
# Use peaks flanking 1 hr npglut SNPs as reference

# init ####
{
  library(Seurat)
  library(Signac)
  
  library(EnsDb.Hsapiens.v86)
  library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  
  library(edgeR)
  library(sva)
  
  library(RColorBrewer)
  library(stringr)
  library(future)

  library(ggplot2)

}

# param #####
plan("multisession", workers = 8)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# load data #####
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/peaks_inside_npglut_1hr_SNPs_500bp.RData")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/npglut_hr1_count_FeatureMatrix_76_lines_08Feb2024.RData")
load("multiomic_obj_new_470K_cells_frag_file.RData")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/metadata_multiomic_obj_470K_cells.RData")

# load("multiomic_obj_with_new_peaks_labeled.RData")
# metadata_multiomic_obj <-
#   multiomic_obj_new@meta.data
# save(metadata_multiomic_obj,
#      file = "metadata_multiomic_obj_470K_cells.RData")

# dedup
npglut_hr1_peaks_count <-
  npglut_hr1_peaks_count[!duplicated(rownames(npglut_hr1_peaks_count)), ]
signac_npglut_hr1_peaks <-
  CreateChromatinAssay(counts = npglut_hr1_peaks_count,
                       sep = c("-", "-"),
                       genome = "hg38",
                       fragments = ATAC_fragment)
Seurat_npglut_hr1_peaks <-
  CreateSeuratObject(counts = signac_npglut_hr1_peaks,
                     assay = "ATAC",
                     meta.data = metadata_multiomic_obj)

DefaultAssay(Seurat_npglut_hr1_peaks)
hg38_annot <-
  GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(hg38_annot) <-
  paste0('chr',
         seqlevels(hg38_annot))
genome(hg38_annot) <-
  'hg38'
Annotation(Seurat_npglut_hr1_peaks[['ATAC']]) <- hg38_annot

Seurat_npglut_hr1_peaks <-
  FindTopFeatures(Seurat_npglut_hr1_peaks,
                  min.cutoff = 10,
                  verbose = T)
Seurat_npglut_hr1_peaks <-
  RunTFIDF(Seurat_npglut_hr1_peaks,
           verbose = T)
Seurat_npglut_hr1_peaks <-
  RunSVD(Seurat_npglut_hr1_peaks,
         verbose = T)
Seurat_npglut_hr1_peaks <-
  RunUMAP(Seurat_npglut_hr1_peaks,
          reduction = 'lsi',
          dims = 2:30,
          reduction.name = 'umap.atac',
          seed.use = 42)

DimPlot(Seurat_npglut_hr1_peaks,
        reduction = 'umap.atac') +
  NoLegend()

unique(Seurat_npglut_hr1_peaks$RNA.cell.type)
subsetted_seurat_npglut_cellss_only <-
  Seurat_npglut_hr1_peaks[, Seurat_npglut_hr1_peaks$RNA.cell.type == "npglut"]

unique(subsetted_seurat_npglut_cellss_only$time.ident)
subsetted_seurat_npglut_cellss_only$time.ident <-
  factor(subsetted_seurat_npglut_cellss_only$time.ident,
         levels = sort(unique(subsetted_seurat_npglut_cellss_only$time.ident)))
unique(subsetted_seurat_npglut_cellss_only$cell.line.ident)
subsetted_seurat_npglut_cellss_only$cell.line.ident <-
  factor(subsetted_seurat_npglut_cellss_only$cell.line.ident,
         levels = sort(unique(subsetted_seurat_npglut_cellss_only$cell.line.ident)))
levels(subsetted_seurat_npglut_cellss_only$cell.line.ident)[1]

# head(rowSums(subsetted_seurat_npglut_cellss_only@assays$ATAC@scale.data))

SNP_flanking_peak_count_split_by_time <-
  vector(mode = "list",
         length = length(levels(subsetted_seurat_npglut_cellss_only$time.ident)))
names(SNP_flanking_peak_count_split_by_time) <-
  levels(subsetted_seurat_npglut_cellss_only$time.ident)

for (i in 1:length(levels(subsetted_seurat_npglut_cellss_only$time.ident))) {
  seurat_2_count <-
    subsetted_seurat_npglut_cellss_only[, 
                                        subsetted_seurat_npglut_cellss_only$time.ident %in%
                                          levels(subsetted_seurat_npglut_cellss_only$time.ident)[i]]
  raw_rowcounts_matrix <-
    matrix(data = 0,
           nrow = nrow(seurat_2_count),
           ncol = length(levels(seurat_2_count$cell.line.ident)),
           dimnames = list(rownames(seurat_2_count),
                           levels(seurat_2_count$cell.line.ident)))
  for (j in 1:length(levels(seurat_2_count$cell.line.ident))) {
    print(paste(i, j))
    subsetted_count_bcmatrix <-
      seurat_2_count[, seurat_2_count$cell.line.ident == levels(seurat_2_count$cell.line.ident)[j]]
    raw_rowcounts_matrix[, j] <-
      unlist(rowSums(subsetted_count_bcmatrix@assays$ATAC@counts))
  }
  SNP_flanking_peak_count_split_by_time[[i]] <-
    raw_rowcounts_matrix
}

SNP_flanking_peak_cpm <-
  vector(mode = "list",
         length = length(levels(subsetted_seurat_npglut_cellss_only$time.ident)))
names(SNP_flanking_peak_cpm) <-
  levels(subsetted_seurat_npglut_cellss_only$time.ident)

meta_line_match_group <-
  metadata_multiomic_obj[!duplicated(metadata_multiomic_obj$cell.line.ident), ]
meta_line_match_group <-
  meta_line_match_group[order(meta_line_match_group$cell.line.ident), ]
meta_line_match_group$batch <-
  str_split(meta_line_match_group$group.ident,
            pattern = "-",
            simplify = T)[, 1]

for (i in 1:length(levels(subsetted_seurat_npglut_cellss_only$time.ident))) {
  print(i)
  raw_count_post_combat <-
    ComBat_seq(counts = SNP_flanking_peak_count_split_by_time[[i]],
               batch = meta_line_match_group$batch,
               full_mod = T)
  cpm_count <-
    DGEList(counts = raw_count_post_combat,
            samples = colnames(SNP_flanking_peak_count_split_by_time[[i]]),
            genes = rownames(SNP_flanking_peak_count_split_by_time[[i]]),
            remove.zeros = F)
  cpm_count <-
    calcNormFactors(cpm_count)
  cpm_count <-
    estimateDisp(cpm_count)
  cpm_count <-
    cpm(cpm_count, 
        log = T,
        prior.count = 0) / 2
  SNP_flanking_peak_cpm[[i]] <-
    as.data.frame(cpm_count)
}

## 1vs0 hr ####
df_2_plot <-
  data.frame(x = rowMeans(SNP_flanking_peak_cpm[[2]],
                          na.rm = T),
             y = rowMeans(SNP_flanking_peak_cpm[[1]],
                          na.rm = T),
             stringsAsFactors = F)
perc_comparable <-
  sum(abs(df_2_plot$y - df_2_plot$x) < 0.5) / nrow(df_2_plot)


# df_polygon_mapping <-
#   data.frame(x = c(0, 0, 0.5,
#                    max(df_2_plot$x) + 0.5,
#                    max(df_2_plot$x) + 0.5,
#                    max(df_2_plot$x)), 
#              y = c(0.5, 0, 0, 
#                    max(df_2_plot$y),
#                    max(df_2_plot$y) + 0.5,
#                    max(df_2_plot$y) + 0.5))
df_polygon_mapping <-
  data.frame(x = c(0, 0, 0.5,
                   4.7 + 0.5,
                   4.7 + 0.5,
                   4.7), 
             y = c(0.5, 0, 0, 
                   4.7,
                   4.7 + 0.5,
                   4.7 + 0.5))
max(df_2_plot$x)
# max(df_2_plot$y)

ggplot(df_2_plot,
       aes(x = x,
           y = y)) +
  geom_point(size = 0.25,
             alpha = 0.3) +
  geom_abline(slope = 1,
              intercept = 0,
              colour = "darkred") +
  xlim(c(0, max(df_2_plot$x) + 0.5)) +
  ylim(c(0, max(df_2_plot$y) + 0.5)) +
  xlab("Accessibility in 1 hr npglut (log2 ATAC-seq peak CPM)") +
  ylab("Accessibility in 0 hr npglut (log2 ATAC-seq peak CPM)") +
  geom_polygon(data = df_polygon_mapping,
               aes(x = x, 
                   y = y),
               fill = "blue", 
               alpha = 0.2, inherit.aes = F) +
  geom_text(x = 1,
            y = 1,
            label = paste0(signif(perc_comparable * 100,
                                  digits = 3),
                           '%')) +
  theme_minimal() +
  ggtitle(paste("npglut 1hr-specific ASoC SNP-flanking peaks",
                "n = 10489",
                sep = "\n"))
