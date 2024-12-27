# Siwei 08 Feb 2024
# Calc chromatin accessibility correlation using Yifan's code

# init ####
{
  library(Seurat)
  library(Signac)
  # library(sctransform)
  # library(glmGamPoi)
  library(EnsDb.Hsapiens.v86)
  # library(GenomeInfoDb)
  library(GenomicFeatures)
  # library(AnnotationDbi)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  
  library(edgeR)
  # library(patchwork)
  # library(readr)
  # library(ggplot2)
  library(RColorBrewer)
  # library(dplyr)
  # library(viridis)
  # library(graphics)
  library(stringr)
  library(future)
  # library(Seurat)
  # library(Signac)
  # library(GenomicRanges)
  # library(BSgenome.Hsapiens.UCSC.hg38)
  # 
  # # library(readr)
  # # library(vcfR)
  # library(stringr)
  
  library(ggplot2)
  # 
  # library(parallel)
  # library(future)
  # 
  # library(MASS)
  # 
  # library(RColorBrewer)
  # library(grDevices)
  # 
  # library(data.table)
  # library(tidyverse)
  # 
  # library(ggvenn)
  
  # library(dplyr)
}

# param #####
plan("multisession", workers = 8)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# load lexi's pre-made 500 bp peak data set ####
load("/nvmefs/scARC_Duan_018/018-029_combined_analysis/df_Lexi_new_peak_set_500bp.RData")
## remove unidentified
peaks_500_width <-
  peaks_500_width[str_detect(string = names(peaks_500_width),
                             pattern = "unidentified",
                             negate = T)]

## merge as one df
df_master_peak_set <-
  do.call(rbind,
          args = peaks_500_width)
sum(is.na(df_master_peak_set$chr))
unique(df_master_peak_set$chr)
df_master_peak_set <-
  df_master_peak_set[str_detect(string = df_master_peak_set$chr,
                                pattern = "^chr[1-9].*"), ]
length(sort(unique(df_master_peak_set$chr)))
## re-sort by chr names
df_master_peak_set$chr <-
  factor(df_master_peak_set$chr,
         levels = sort(unique(df_master_peak_set$chr)))

## split by chr name
all_peaks_by_chr <-
  split(x = df_master_peak_set,
        f = df_master_peak_set$chr)
names(all_peaks_by_chr)

discontinuous_peaks_list_by_chr <-
  vector(mode = "list",
         length = length(all_peaks_by_chr))
names(discontinuous_peaks_list_by_chr) <-
  names(all_peaks_by_chr)

for (i in 1:length(all_peaks_by_chr)) {
  peaks_by_chr <-
    all_peaks_by_chr[[i]]
  peaks_by_chr <-
    peaks_by_chr[order(peaks_by_chr$start), ]
  for (j in 1:(nrow(peaks_by_chr))) {
    print(paste0("chr",
                 i,
                 ", peak ",
                 j))
    current_peak <-
      peaks_by_chr[j, ]
    if (j == 1) {
      previous_peak <- current_peak
      df_combined_peak <- current_peak
    } else {
      if (previous_peak$end <= current_peak$start) {
        df_combined_peak <-
          as.data.frame(rbind(df_combined_peak,
                              current_peak))
        previous_peak <- current_peak
      } else {
        if (previous_peak$score < current_peak$score) {
          df_combined_peak[nrow(df_combined_peak), ] <-
            unlist(current_peak)
          previous_peak <- current_peak
        }
      }
    }
  }
  discontinuous_peaks_list_by_chr[[i]] <-
    df_combined_peak
}

full_discontinuous_peaks_list <-
  do.call(rbind,
          args = discontinuous_peaks_list_by_chr)
sum(is.na(full_discontinuous_peaks_list$chr))

df_2_write <-
  full_discontinuous_peaks_list
df_2_write$chr <-
  str_split(string = rownames(df_2_write),
            pattern = "\\.",
            simplify = T)[, 1]
unique(df_2_write$chr)
# df_2_write <-
#   df_2_write[, c(7, 2:6)]
df_2_write$start <-
  as.integer(df_2_write$start)
df_2_write$end <-
  as.integer(df_2_write$end)
df_2_write$strand <- '*'

GRanges_all_peaks <-
  makeGRangesFromDataFrame(df = df_2_write,
                           keep.extra.columns = T, 
                           na.rm = T)

# load all SNPs ####
load("~/NVME/scARC_Duan_018/R_ASoC/Sum_100_lines_21Dec2023.RData")

## extract the lists w ASoC SNPs only ([[3]])
vcf_lst_ASoC_only <-
  vector(mode = "list",
         length = length(master_vcf_list))

for (i in 1:length(vcf_lst_ASoC_only)) {
  print(i)
  vcf_lst_ASoC_only[[i]] <-
    master_vcf_list[[i]][[3]]
}
names(vcf_lst_ASoC_only) <-
  names(master_vcf_list)

full_ASoC_SNP_list <-
  do.call(rbind,
          args = vcf_lst_ASoC_only)
full_ASoC_SNP_list <-
  full_ASoC_SNP_list[!duplicated(full_ASoC_SNP_list$ID), ]

npglut_0hr_specific_interval <-
  npglut_0hr_specific_SNPs[, c(1, 2, 3, 4, 16)]
npglut_0hr_specific_interval$POS <-
  as.numeric(npglut_0hr_specific_interval$POS) - 1
npglut_0hr_specific_interval$POS.1 <-
  as.numeric(npglut_0hr_specific_interval$POS.1)
npglut_0hr_specific_interval$strand <- "*"
colnames(npglut_0hr_specific_interval)[1:3] <-
  c("chr", "start", "end")
GRanges_npglut_0hr_SNP_interval <-
  makeGRangesFromDataFrame(df = npglut_0hr_specific_interval,
                           keep.extra.columns = T,
                           na.rm = T)

peaks_flanking_npglut_0hr_SNPs <-
  GRanges_all_peaks[GenomicRanges::nearest(x = GRanges_npglut_0hr_SNP_interval,
                                           subject = GRanges_all_peaks,
                                           ignore.strand = T,
                                           select = "arbitrary")]

peaks_flanking_npglut_0hr_SNPs <-
  subsetByOverlaps(x = peaks_flanking_npglut_0hr_SNPs,
                   ranges = blacklist_hg38_unified,
                   invert = T)

load("multiomic_obj_with_new_peaks_labeled.RData")

multiomic_obj_new <-
  UpdateSeuratObject(multiomic_obj_new)
DefaultAssay(multiomic_obj_new) <- "ATAC"
ATAC_fragment <-
  Fragments(multiomic_obj_new)

## load from here
# save(list = c("ATAC_fragment",
#               "peaks_flanking_npglut_0hr_SNPs"),
#      file = "data_2_count_npglut_hr0_SNP_flanking_peaks.RData")
load("data_2_count_npglut_hr0_SNP_flanking_peaks.RData")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/df_Lexi_cells_470K.RData")

system.time({
  npglut_hr0_peaks_count <-
    FeatureMatrix(fragments = ATAC_fragment, 
                  features = peaks_flanking_npglut_0hr_SNPs,
                  process_n = 2000, 
                  cells = df_Lexi_cells_470K,
                  verbose = T)
})
print("counting done")
save(npglut_hr0_peaks_count,
     file = "npglut_hr0_count_FeatureMatrix_76_lines_08Feb2024.RData")

load("multiomic_obj_with_new_peaks_labeled.RData")

head(colnames(npglut_hr0_peaks_count))
head(colnames(multiomic_obj_new))

dedup_npglut_hr0_peaks_count <-
  npglut_hr0_peaks_count[!duplicated(rownames(npglut_hr0_peaks_count)), ]
signac_npglut_hr0_peaks <-
  CreateChromatinAssay(counts = dedup_npglut_hr0_peaks_count,
                       sep = c("-", "-"),
                       genome = "hg38",
                       fragments = ATAC_fragment)
Seurat_npglut_hr0_peaks <-
  CreateSeuratObject(counts = signac_npglut_hr0_peaks,
                     assay = "ATAC",
                     meta.data = multiomic_obj_new@meta.data)

DefaultAssay(Seurat_npglut_hr0_peaks)
hg38_annot <-
  GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(hg38_annot) <-
  paste0('chr',
         seqlevels(hg38_annot))
genome(hg38_annot) <-
  'hg38'
Annotation(Seurat_npglut_hr0_peaks[['ATAC']]) <- hg38_annot

Seurat_npglut_hr0_peaks <-
  FindTopFeatures(Seurat_npglut_hr0_peaks,
                  min.cutoff = 10,
                  verbose = T)
Seurat_npglut_hr0_peaks <-
  RunTFIDF(Seurat_npglut_hr0_peaks,
           verbose = T)
Seurat_npglut_hr0_peaks <-
  RunSVD(Seurat_npglut_hr0_peaks,
         verbose = T)
Seurat_npglut_hr0_peaks <-
  RunUMAP(Seurat_npglut_hr0_peaks,
          reduction = 'lsi',
          dims = 2:30,
          reduction.name = 'umap.atac',
          seed.use = 42)

DimPlot(Seurat_npglut_hr0_peaks,
        reduction = 'umap.atac') +
  NoLegend()

unique(Seurat_npglut_hr0_peaks$RNA.cell.type)
subsetted_seurat_npglut_cellss_only <-
  Seurat_npglut_hr0_peaks[, Seurat_npglut_hr0_peaks$RNA.cell.type == "npglut"]

unique(subsetted_seurat_npglut_cellss_only$time.ident)
subsetted_seurat_npglut_cellss_only$time.ident <-
  factor(subsetted_seurat_npglut_cellss_only$time.ident,
         levels = sort(unique(subsetted_seurat_npglut_cellss_only$time.ident)))
unique(subsetted_seurat_npglut_cellss_only$cell.line.ident)
subsetted_seurat_npglut_cellss_only$cell.line.ident <-
  factor(subsetted_seurat_npglut_cellss_only$cell.line.ident,
         levels = sort(unique(subsetted_seurat_npglut_cellss_only$cell.line.ident)))
levels(subsetted_seurat_npglut_cellss_only$cell.line.ident)[1]

head(rowSums(subsetted_seurat_npglut_cellss_only@assays$ATAC@scale.data))

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

for (i in 1:length(levels(subsetted_seurat_npglut_cellss_only$time.ident))) {
  print(i)
  cpm_count <-
    DGEList(counts = SNP_flanking_peak_count_split_by_time[[i]],
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
  data.frame(x = rowMeans(SNP_flanking_peak_cpm[[1]],
                             na.rm = T),
             y = rowMeans(SNP_flanking_peak_cpm[[2]],
                             na.rm = T),
             stringsAsFactors = F)
perc_comparable <-
  sum(abs(df_2_plot$y - df_2_plot$x) < 0.5) / nrow(df_2_plot)


df_polygon_mapping <-
  data.frame(x = c(0, 0, 0.5,
                   max(df_2_plot$x) + 0.5,
                   max(df_2_plot$x) + 0.5,
                   max(df_2_plot$x)), 
             y = c(0.5, 0, 0, 
                   max(df_2_plot$y),
                   max(df_2_plot$y) + 0.5,
                   max(df_2_plot$y) + 0.5))


ggplot(df_2_plot,
       aes(x = x,
           y = y)) +
  geom_point(size = 0.25,
             alpha = 0.3) +
  geom_abline(slope = 1,
              intercept = 0,
              colour = "darkred") +
  xlim(c(0, max(df_2_plot$x) + 0.5)) +
  ylim(c(0, max(df_2_plot$y + 0.5))) +
  xlab("Accessibility in 0 hr npglut (log2 ATAC-seq peak CPM)") +
  ylab("Accessibility in 1 hr npglut (log2 ATAC-seq peak CPM)") +
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
  ggtitle(paste("npglut-specific ASoC SNP-flanking peaks",
                "n = 3830",
                sep = "\n"))

## 6vs0 hr #####
df_2_plot <-
  data.frame(x = rowMeans(SNP_flanking_peak_cpm[[1]],
                          na.rm = T),
             y = rowMeans(SNP_flanking_peak_cpm[[3]],
                          na.rm = T),
             stringsAsFactors = F)
perc_comparable <-
  sum(abs(df_2_plot$y - df_2_plot$x) < 0.5) / nrow(df_2_plot)


df_polygon_mapping <-
  data.frame(x = c(0, 0, 0.5,
                   max(df_2_plot$x) + 0.5,
                   max(df_2_plot$x) + 0.5,
                   max(df_2_plot$x)), 
             y = c(0.5, 0, 0, 
                   max(df_2_plot$y),
                   max(df_2_plot$y) + 0.5,
                   max(df_2_plot$y) + 0.5))


ggplot(df_2_plot,
       aes(x = x,
           y = y)) +
  geom_point(size = 0.25,
             alpha = 0.3) +
  geom_abline(slope = 1,
              intercept = 0,
              colour = "darkred") +
  xlim(c(0, max(df_2_plot$x) + 0.5)) +
  ylim(c(0, max(df_2_plot$y + 0.5))) +
  xlab("Accessibility in 0 hr npglut (log2 ATAC-seq peak CPM)") +
  ylab("Accessibility in 6 hr npglut (log2 ATAC-seq peak CPM)") +
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
  ggtitle(paste("npglut-specific ASoC SNP-flanking peaks",
                "n = 3830",
                sep = "\n"))



# load data ####
load("~/NVME/scARC_Duan_018/R_ASoC/npglut_0hr_specific_peaks_list.RData")
load("~/NVME/scARC_Duan_018/R_ASoC/npglut_0hr_specific_SNP_list.RData")

npglut_0hr_specific_peaks <-
  keepStandardChromosomes(npglut_0hr_specific_peaks,
                          pruning.mode = "coarse")
npglut_0hr_specific_peaks <-
  subsetByOverlaps(x = npglut_0hr_specific_peaks,
                   ranges = blacklist_hg38_unified,
                   invert = T)

all_peaks_list <-
  as.data.frame(npglut_0hr_specific_peaks)
all_peaks_list <-
  all_peaks_list[!(all_peaks_list$seqnames %in% c("chrX", "chrY", "chrM")), ]

all_peaks_list$start_500 <-
  all_peaks_list$start + all_peaks_list$relative_summit_position - 250
all_peaks_list$end_500 <-
  all_peaks_list$start + all_peaks_list$relative_summit_position + 250

all_peaks_by_chr <-
  split(x = all_peaks_list,
        f = all_peaks_list$seqnames)
all_peaks_by_chr <-
  all_peaks_by_chr[1:22]

## make 500 bp discontinuous non-overlapping peak sets
discontinuous_peaks_list_by_chr <-
  vector(mode = "list",
         length = length(all_peaks_by_chr))
names(discontinuous_peaks_list_by_chr) <-
  names(all_peaks_by_chr)

for (i in 1:length(all_peaks_by_chr)) {
  peaks_by_chr <-
    all_peaks_by_chr[[i]]
  peaks_by_chr <-
    peaks_by_chr[order(peaks_by_chr$start_500), ]
  for (j in 1:(nrow(peaks_by_chr))) {
    print(paste0("chr",
                 i,
                 ", peak ",
                 j))
    current_peak <-
      peaks_by_chr[j, ]
    if (j == 1) {
      previous_peak <- current_peak
      df_combined_peak <- current_peak
    } else {
      if (previous_peak$end_500 <= current_peak$start_500) {
        df_combined_peak <-
          as.data.frame(rbind(df_combined_peak,
                              current_peak))
        previous_peak <- current_peak
      } else {
        if (previous_peak$score < current_peak$score) {
          df_combined_peak[nrow(df_combined_peak), ] <-
            unlist(current_peak)
          previous_peak <- current_peak
        }
      }
    }
  }
  discontinuous_peaks_list_by_chr[[i]] <-
    df_combined_peak
}

full_discontinuous_peaks_list <-
  do.call(rbind,
          args = discontinuous_peaks_list_by_chr)
unique(full_discontinuous_peaks_list$seqnames)
sum(is.na(full_discontinuous_peaks_list$seqnames))
full_discontinuous_peaks_list[is.na(full_discontinuous_peaks_list$seqnames), ]

# write discontinuous linear peaks as bed
df_2_write <-
  full_discontinuous_peaks_list[, c(1, 13, 14, 6, 7, 5)]
df_2_write$chr <-
  str_split(string = rownames(df_2_write),
            pattern = "\\.",
            simplify = T)[, 1]
unique(df_2_write$chr)
df_2_write <-
  df_2_write[, c(7, 2:6)]
df_2_write$start_500 <-
  as.integer(df_2_write$start_500)
df_2_write$end_500 <-
  as.integer(df_2_write$end_500)
df_2_write$strand <- '*'

npglut_0hr_specific_peaks_500bp <-
  df_2_write
save(npglut_0hr_specific_peaks_500bp,
     file = "npglut_0hr_specific_peaks_500bp_08Feb2024.RData")

df_4_GRanges <-
  npglut_0hr_specific_peaks_500bp[, c(1:3, 6)]
colnames(df_4_GRanges) <-
  c("chr", "start", "end", "strand")
GRanges_npglut_0hr <-
  makeGRangesFromDataFrame(df = df_4_GRanges,
                           na.rm = T)
save(GRanges_npglut_0hr,
     file = "Granges_npglut_0hr_new_peaks_500bp.RData")

# Start from here #####
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/npglut_0hr_specific_peaks_500bp_08Feb2024.RData")
load("multiomic_obj_with_new_peaks_labeled.RData")

multiomic_obj_new <-
  UpdateSeuratObject(multiomic_obj_new)
DefaultAssay(multiomic_obj_new) <- "ATAC"
ATAC_fragment <-
  Fragments(multiomic_obj_new)

colnames(npglut_0hr_specific_peaks_500bp) <-
  c("chr", "start", "end", "seqinfo", "score", "strand")
npglut_0hr_features <-
  makeGRangesFromDataFrame(df = npglut_0hr_specific_peaks_500bp,
                           # seqinfo = npglut_0hr_specific_peaks_500bp$seqinfo,
                           keep.extra.columns = T, 
                           na.rm = T)
# save(npglut_0hr_features,
#      file = "npglut_0hr_features_500bp_4_FeatureMatrix.RData")

npglut_hr0_peaks_count <-
  FeatureMatrix(fragments = ATAC_fragment, 
                features = npglut_0hr_features,
                process_n = 4000, 
                # cells = cells, 
                verbose = T)
save(npglut_hr0_peaks_count,
     file = "npglut_hr0_count_FeatureMatrix_76_lines_08Feb2024.RData")
