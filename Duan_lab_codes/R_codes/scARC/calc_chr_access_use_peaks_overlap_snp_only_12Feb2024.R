# Siwei 12 Feb 2024
# Calc chromatin accessibility correlation using Yifan's code
# use peaks from npglut 1hr and 6 hr
# use peaks overlapping SNP intervals only

# init ####
{
  library(Seurat)
  library(Signac)
  
  library(EnsDb.Hsapiens.v86)
  library(GenomicFeatures)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  
  library(edgeR)
  
  library(RColorBrewer)
  
  library(stringr)
  library(future)
  
  library(ggplot2)
  
  library(parallel)
  library(doParallel)
  library(future)
  library(foreach)
  
  library(Rfast)
  # library()
  
  library(MASS)
}

# param #####
plan("multisession", workers = 8)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# process npglut 1 hr #####
## load intervals and FeatureMatrix reads
load("peaks_inside_npglut_1hr_SNPs_500bp.RData")
load("npglut_hr1_peaks_count_FeatureMatrix_76_lines_12Feb2024.RData")
load("multiomic_obj_new_470K_cells_frag_file.RData")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/metadata_multiomic_obj_470K_cells.RData")

## cut out npglut metadata, 1 hr and 6 hr should have the same cell number
## so should be able to create one Seurat object contains two ChromatinAssay?
metadata_multiomic_obj_npglut <-
  metadata_multiomic_obj[metadata_multiomic_obj$RNA.cell.type == "npglut", ]

npglut_hr1_peaks_count <-
  npglut_hr1_peaks_count[!duplicated(rownames(npglut_hr1_peaks_count)), ]
signac_npglut <-
  CreateChromatinAssay(counts = npglut_hr1_peaks_count,
                       sep = c("-", "-"),
                       genome = "hg38",
                       fragments = ATAC_fragment)
head(colnames(signac_npglut))
head(rownames(metadata_multiomic_obj_npglut))
Seurat_npglut <-
  CreateSeuratObject(counts = signac_npglut,
                     assay = "hr1_npglut_ATAC",
                     meta.data = metadata_multiomic_obj_npglut)

DefaultAssay(Seurat_npglut)
hg38_annot <-
  GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(hg38_annot) <-
  paste0('chr',
         seqlevels(hg38_annot))
genome(hg38_annot) <-
  'hg38'

Annotation(Seurat_npglut[['hr1_npglut_ATAC']]) <- hg38_annot

## TFIDF+SVD is not necessary if not visualising peaks by UMAP ####
## It is nevetheless recommended to run if there are not too many peaks...
Seurat_npglut <-
  FindTopFeatures(Seurat_npglut,
                  min.cutoff = 10,
                  verbose = T)
Seurat_npglut <-
  RunTFIDF(Seurat_npglut,
           verbose = T)
Seurat_npglut <-
  RunSVD(Seurat_npglut,
         verbose = T)
Seurat_npglut <-
  RunUMAP(Seurat_npglut,
          reduction = 'lsi',
          dims = 2:30,
          reduction.name = 'umap.atac',
          seed.use = 42)

DimPlot(Seurat_npglut,
        reduction = 'umap.atac') +
  NoLegend()

## confirm there is only one cell type
unique(Seurat_npglut$RNA.cell.type)

## for consistency we will directly adopt codes here
## though this will reproduce a same object anyway
subsetted_seurat_npglut_cellss_only <-
  Seurat_npglut[, Seurat_npglut$RNA.cell.type == "npglut"]

# unique(subsetted_seurat_npglut_cellss_only$time.ident)
# subsetted_seurat_npglut_cellss_only$time.ident <-
#   factor(subsetted_seurat_npglut_cellss_only$time.ident,
#          levels = sort(unique(subsetted_seurat_npglut_cellss_only$time.ident)))
unique(subsetted_seurat_npglut_cellss_only$cell.line.ident)
subsetted_seurat_npglut_cellss_only$cell.line.ident <-
  factor(subsetted_seurat_npglut_cellss_only$cell.line.ident,
         levels = sort(unique(subsetted_seurat_npglut_cellss_only$cell.line.ident)))
levels(subsetted_seurat_npglut_cellss_only$cell.line.ident)[1]

subsetted_seurat_npglut_cellss_only$time.ident <-
  factor(subsetted_seurat_npglut_cellss_only$time.ident,
         levels = c("0hr", "1hr", "6hr"))
SNP_flanking_peak_count_split_by_time <-
  vector(mode = "list",
         length = length(unique(subsetted_seurat_npglut_cellss_only$time.ident)))
names(SNP_flanking_peak_count_split_by_time) <-
  levels(subsetted_seurat_npglut_cellss_only$time.ident)
names(SNP_flanking_peak_count_split_by_time)

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
      unlist(rowSums(subsetted_count_bcmatrix@assays$hr1_npglut_ATAC@counts))
  }
  SNP_flanking_peak_count_split_by_time[[i]] <-
    raw_rowcounts_matrix
  rm(seurat_2_count)
  rm(subsetted_count_bcmatrix)
  rm(raw_rowcounts_matrix)
}

SNP_flanking_peak_cpm <-
  vector(mode = "list",
         length = length(levels(subsetted_seurat_npglut_cellss_only$time.ident)))
names(SNP_flanking_peak_cpm) <-
  levels(subsetted_seurat_npglut_cellss_only$time.ident)

## Skip sva_combat as this appears to smooth data when the counts is low...
for (i in 1:length(levels(subsetted_seurat_npglut_cellss_only$time.ident))) {
  print(i)
  raw_count_post_combat <-
    SNP_flanking_peak_count_split_by_time[[i]]
  # raw_count_post_combat <-
  #   trunc(raw_count_post_combat / 5)
  # print(colSums(raw_count_post_combat))
  raw_count_post_combat <-
    raw_count_post_combat[, (colSums(raw_count_post_combat) > 50000)]
  #                              colSums(raw_count_post_combat) < 20000)]
  print(nrow(raw_count_post_combat))
  #   ComBat_seq(counts = SNP_flanking_peak_count_split_by_time[[i]],
  #              batch = meta_line_match_group$batch,
  #              full_mod = T)
  cpm_count <-
    DGEList(counts = raw_count_post_combat,
            # samples = colnames(SNP_flanking_peak_count_split_by_time[[i]]),
            samples = colnames(raw_count_post_combat),
            genes = rownames(SNP_flanking_peak_count_split_by_time[[i]]),
            remove.zeros = F)
  cpm_count <-
    rpkm.DGEList(y = cpm_count,
                 gene.length = 500,
                 log = F,
                 prior.count = 1)
  cpm_count <-
    log2(cpm_count + 1)
  # cpm_count <-
  #   calcNormFactors(cpm_count)
  # cpm_count <-
  #   estimateDisp(cpm_count)
  # cpm_count <-
  #   cpm(cpm_count, 
  #       log = T,
  #       prior.count = 1) / 2
  SNP_flanking_peak_cpm[[i]] <-
    as.data.frame(cpm_count)
}

## find common libraries
common_lib_names <-
  sort(intersect(colnames(SNP_flanking_peak_cpm[[1]]),
                 colnames(SNP_flanking_peak_cpm[[2]])))

## 1vs0 hr, note use [[2]] ####
df_2_plot <-
  data.frame(x = rowMeans(SNP_flanking_peak_cpm[[2]],
                          na.rm = T),
             y = rowMeans(SNP_flanking_peak_cpm[[1]],
                          na.rm = T),
             stringsAsFactors = F)
# df_2_plot <-
#   data.frame(x = rowMins(SNP_flanking_peak_cpm[[2]],
#                           value = T),
#              y = rowMins(SNP_flanking_peak_cpm[[1]],
#                           value = T),
#              stringsAsFactors = F)
perc_comparable <-
  sum(abs(df_2_plot$y - df_2_plot$x) < 0.5) / nrow(df_2_plot)

df_polygon_mapping <-
  data.frame(x = c(0, 0, 1,
                   11 + 1,
                   11 + 1,
                   1), 
             y = c(1, 0, 0, 
                   1,
                   1 + 1,
                   1 + 1))
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
  ggtitle(paste("npglut 1hr-specific ASoC SNP-flanking peaks\n",
                "nPeaks =",
                nrow(df_2_plot),
                sep = " "))

View(SNP_flanking_peak_count_split_by_time[[1]])
View(SNP_flanking_peak_count_split_by_time[[2]])

write.table(SNP_flanking_peak_count_split_by_time[[1]],
            file = "npglut_0hr_interval_raw_counts.txt",
            quote = F, sep = "\t",
            row.names = T, col.names = T)
write.table(SNP_flanking_peak_count_split_by_time[[2]],
            file = "npglut_1hr_interval_raw_counts.txt",
            quote = F, sep = "\t",
            row.names = T, col.names = T)

length(as.integer(raw_count_post_combat / 10))


# split by library ####
subsetted_seurat_npglut_cellss_only <-
  Seurat_npglut[, Seurat_npglut$RNA.cell.type == "npglut"]

unique(subsetted_seurat_npglut_cellss_only$time.ident)
subsetted_seurat_npglut_cellss_only$time.ident <-
  factor(subsetted_seurat_npglut_cellss_only$time.ident,
         levels = sort(unique(subsetted_seurat_npglut_cellss_only$time.ident)))
# unique(subsetted_seurat_npglut_cellss_only$cell.line.ident)
# subsetted_seurat_npglut_cellss_only$cell.line.ident <-
#   factor(subsetted_seurat_npglut_cellss_only$cell.line.ident,
#          levels = sort(unique(subsetted_seurat_npglut_cellss_only$cell.line.ident)))
# levels(subsetted_seurat_npglut_cellss_only$cell.line.ident)[1]
# unique(subsetted_seurat_npglut_cellss_only$group.ident)
# subsetted_seurat_npglut_cellss_only$group.ident <-
#   factor(subsetted_seurat_npglut_cellss_only$group.ident,
#          levels = sort(unique(subsetted_seurat_npglut_cellss_only$group.ident)))
# levels(subsetted_seurat_npglut_cellss_only$group.ident)[1]

subsetted_seurat_npglut_cellss_only$library.ident <-
  str_split(string = subsetted_seurat_npglut_cellss_only$group.ident,
            pattern = "-",
            simplify = T)[, 1]
subsetted_seurat_npglut_cellss_only$library.ident <-
  factor(subsetted_seurat_npglut_cellss_only$library.ident,
         levels = sort(unique(subsetted_seurat_npglut_cellss_only$library.ident)))
levels(subsetted_seurat_npglut_cellss_only$library.ident)[1]
unique(subsetted_seurat_npglut_cellss_only$library.ident)

# SNP_flanking_peak_count_split_by_group <-
#   vector(mode = "list",
#          length = length(levels(subsetted_seurat_npglut_cellss_only$group.ident)))
# names(SNP_flanking_peak_count_split_by_group) <-
#   levels(subsetted_seurat_npglut_cellss_only$group.ident)
# unique(names(SNP_flanking_peak_count_split_by_group))

SNP_flanking_peak_count_split_by_time <-
  vector(mode = "list",
         length = length(levels(subsetted_seurat_npglut_cellss_only$time.ident)))
names(SNP_flanking_peak_count_split_by_time) <-
  levels(subsetted_seurat_npglut_cellss_only$time.ident)
unique(names(SNP_flanking_peak_count_split_by_time))

for (i in 1:length(levels(subsetted_seurat_npglut_cellss_only$time.ident))) {
  seurat_2_count <-
    subsetted_seurat_npglut_cellss_only[, 
                                        subsetted_seurat_npglut_cellss_only$time.ident %in%
                                          levels(subsetted_seurat_npglut_cellss_only$time.ident)[i]]
  raw_rowcounts_matrix <-
    matrix(data = 0,
           nrow = nrow(seurat_2_count),
           ncol = length(levels(seurat_2_count$library.ident)),
           dimnames = list(rownames(seurat_2_count),
                           levels(seurat_2_count$library.ident)))
  for (j in 1:length(levels(seurat_2_count$library.ident))) {
    print(paste(i, j))
    print(levels(seurat_2_count$library.ident)[j])
    if (!(j %in% c(6))) {
      subsetted_count_bcmatrix <-
        seurat_2_count[, seurat_2_count$library.ident == levels(seurat_2_count$library.ident)[j]]
      raw_rowcounts_matrix[, j] <-
        unlist(rowSums(subsetted_count_bcmatrix@assays$hr1_npglut_ATAC@counts))
    }

  }
  SNP_flanking_peak_count_split_by_time[[i]] <-
    raw_rowcounts_matrix
  rm(seurat_2_count)
  rm(subsetted_count_bcmatrix)
  rm(raw_rowcounts_matrix)
}

SNP_flanking_peak_cpm <-
  vector(mode = "list",
         length = length(levels(subsetted_seurat_npglut_cellss_only$time.ident)))
names(SNP_flanking_peak_cpm) <-
  levels(subsetted_seurat_npglut_cellss_only$time.ident)

## Skip sva_combat as this appears to smooth data when the counts is low...
for (i in 1:length(levels(subsetted_seurat_npglut_cellss_only$time.ident))) {
  print(i)
  raw_count_post_combat <-
    SNP_flanking_peak_count_split_by_time[[i]]
  # raw_count_post_combat <-
  #   trunc(raw_count_post_combat / 5)
  raw_count_post_combat <-
    raw_count_post_combat[, (colSums(raw_count_post_combat) > 100000)]
  print(colSums(raw_count_post_combat))
  print(nrow(raw_count_post_combat))
  #   ComBat_seq(counts = SNP_flanking_peak_count_split_by_time[[i]],
  #              batch = meta_line_match_group$batch,
  #              full_mod = T)
  SNP_flanking_peak_cpm[[i]] <-
    raw_count_post_combat

}

## find common libraries
common_lib_names <-
  sort(intersect(intersect(colnames(SNP_flanking_peak_cpm[[1]]),
                             colnames(SNP_flanking_peak_cpm[[2]])),
                   colnames(SNP_flanking_peak_cpm[[3]])))

lib_size_matrix <-
  matrix(data = 0,
         nrow = length(SNP_flanking_peak_cpm),
         ncol = length(common_lib_names))

## refill list with common library names and calc lib size
for (i in 1:length(levels(subsetted_seurat_npglut_cellss_only$time.ident))) {
  print(i)
  raw_count_post_combat <-
    SNP_flanking_peak_count_split_by_time[[i]]
  # raw_count_post_combat <-
  #   trunc(raw_count_post_combat / 5)
  # raw_count_post_combat <-
  #   raw_count_post_combat[, (colSums(raw_count_post_combat) > 100000)]
  raw_count_post_combat <-
    raw_count_post_combat[, colnames(raw_count_post_combat) %in% common_lib_names]
  print(colSums(raw_count_post_combat))
  lib_size_matrix[i, ] <-
    unlist(colSums(raw_count_post_combat))
  print(nrow(raw_count_post_combat))
  # cpm_count <-
  #   DGEList(counts = raw_count_post_combat,
  #           # samples = colnames(SNP_flanking_peak_count_split_by_time[[i]]),
  #           samples = colnames(raw_count_post_combat),
  #           genes = rownames(SNP_flanking_peak_count_split_by_time[[i]]),
  #           remove.zeros = F)
  # cpm_count <-
  #   rpkm.DGEList(y = cpm_count,
  #                gene.length = 500,
  #                log = T,
  #                prior.count = 1) / 2
  # cpm_count <-
  #   calcNormFactors(cpm_count)
  # cpm_count <-
  #   estimateDisp(cpm_count)
  # cpm_count <-
  #   cpm(cpm_count, 
  #       log = T,
  #       prior.count = 1) / 2
  SNP_flanking_peak_cpm[[i]] <-
    raw_count_post_combat
  
}
# View(SNP_flanking_peak_cpm[["0hr"]])

SNP_flanking_peak_manual_calc <-
  vector(mode = "list",
         length = length(levels(subsetted_seurat_npglut_cellss_only$time.ident)))
names(SNP_flanking_peak_manual_calc) <-
  levels(subsetted_seurat_npglut_cellss_only$time.ident)

for (i in 1:length(levels(subsetted_seurat_npglut_cellss_only$time.ident))) {
  print(i)
  current_count_mx <-
    SNP_flanking_peak_cpm[[i]]
  results_mx <-
    matrix(data = 0,
           nrow = nrow(current_count_mx),
           ncol = ncol(current_count_mx),
           dimnames = list(rownames(current_count_mx),
                           colnames(current_count_mx)))
  current_count_mx <-
    as.matrix(current_count_mx)
  # sum_library_size <-
  #   colSums(current_count_mx)
  # for (j in 1:ncol(current_count_mx)) {
  #   print(j)
  #   results_mx[, j] <-
  #     log2((current_count_mx[, j]) / sum_library_size[j] * 1e6 + 1) / 2
  # }
  results_mx <-
    rpkm(current_count_mx,
         gene.length = 500,
         prior.count = 0,
         log = F)
  SNP_flanking_peak_manual_calc[[i]] <-
    log2(results_mx + 1)
}
hist(SNP_flanking_peak_manual_calc[["0hr"]])
hist(rowMeans(SNP_flanking_peak_manual_calc[["0hr"]]))
hist(rowMeans(SNP_flanking_peak_manual_calc[["1hr"]]))

View(SNP_flanking_peak_manual_calc[["0hr"]])
## final calc
# for (i in 1:length(levels(subsetted_seurat_npglut_cellss_only$time.ident))) {
#   print(i)
#   raw_count_post_combat <-
#     SNP_flanking_peak_cpm[[i]]
#   # raw_count_post_combat <-
#   #   raw_count_post_combat * (lib_size_matrix[i, ] / colMeans(lib_size_matrix))
#   # raw_count_post_combat <-
#   #   raw_count_post_combat[, c(1:6, 10, 12, 13, 15)]
#   print(ncol(raw_count_post_combat))
#   # raw_count_post_combat <-
#   #   trunc(raw_count_post_combat * (lib_size_matrix[i, ] / colMeans(lib_size_matrix)))
#     # SNP_flanking_peak_count_split_by_time[[i]]
#   # raw_count_post_combat <-
#   #   trunc(raw_count_post_combat / 5)
#   # raw_count_post_combat <-
#   #   raw_count_post_combat[, (colSums(raw_count_post_combat) > 100000)]
#   # raw_count_post_combat <-
#   #   raw_count_post_combat[, colnames(raw_count_post_combat) %in% common_lib_names]
#   # print(colSums(raw_count_post_combat))
#   # lib_size_matrix[i, ] <-
#   #   unlist(colSums(raw_count_post_combat))
#   print(nrow(raw_count_post_combat))
#   # cpm_count <-
#   #   DGEList(counts = as.matrix(raw_count_post_combat),
#   #           # samples = colnames(SNP_flanking_peak_count_split_by_time[[i]]),
#   #           samples = colnames(raw_count_post_combat),
#   #           genes = rownames(SNP_flanking_peak_count_split_by_time[[i]]),
#   #           remove.zeros = F)
#   # cpm_count <-
#   #   rpkm.DGEList(y = cpm_count,
#   #                gene.length = 500,
#   #                log = T,
#   #                prior.count = 1) / 2
#   cpm_count <-
#     cpm(as.matrix(raw_count_post_combat) + 1,
#         log = T,
#         normalzed.lib.sizes = F,
#         # shrunk = F,
#         prior.count = 0) / 2
#   # cpm_count <-
#   #   cpm_count * log2((lib_size_matrix[i, ] / colMeans(lib_size_matrix)))
#   # cpm_count <-
#   #   calcNormFactors(cpm_count)
#   # cpm_count <-
#   #   estimateDisp(cpm_count)
#   # cpm_count <-
#   #   cpm(cpm_count, 
#   #       log = T,
#   #       prior.count = 1) / 2
#   # print(cpm_count)
#   # SNP_flanking_peak_cpm[[i]] <-
#   #   as.data.frame(raw_count_post_combat)
#   SNP_flanking_peak_cpm[[i]] <-
#     cpm_count
#   
# }

View(SNP_flanking_peak_cpm[["0hr"]])

## 1vs0 hr, note use [[2]] ####
df_2_plot <-
  data.frame(x = rowMeans(SNP_flanking_peak_manual_calc[[2]],
                          na.rm = T),
             y = rowMeans(SNP_flanking_peak_manual_calc[[1]],
                          na.rm = T),
             stringsAsFactors = F)
perc_comparable <-
  sum(abs(df_2_plot$y - df_2_plot$x) < 1) / nrow(df_2_plot)

# df_polygon_mapping <-
#   data.frame(x = c(0, 0, 0.5,
#                    5.5 + 0.5,
#                    5.5 + 0.5,
#                    5.5), 
#              y = c(0.5, 0, 0, 
#                    5.5,
#                    5.5 + 0.5,
#                    5.5 + 0.5))
df_polygon_mapping <-
  data.frame(x = c(0, 0, 1,
                   10.9 + 1,
                   10.9 + 1,
                   10.9), 
             y = c(1, 0, 0, 
                   10.9,
                   10.9 + 1,
                   10.9 + 1))
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
  xlab("Accessibility in 1 hr npglut (log2 ATAC-seq peak RPKM)") +
  ylab("Accessibility in 0 hr npglut (log2 ATAC-seq peak RPKM)") +
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
  ggtitle(paste("npglut 1hr-specific ASoC SNP-flanking peaks\n",
                "nPeaks =",
                nrow(df_2_plot),
                sep = " "))
