# Chuxuan Li 11/23/2021
# Using Signac to analyze group_2 and group_8 ATACseq dataset aggregated by 10x aggr
#with calling peaks on specific cell types found in 11Oct2021_ATACseq_clustering.R
#then redo the analysis with the new set of peaks


# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)

set.seed(1886)
# set threads and parallelization
plan("multisession", workers = 2)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 11474836480)

load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v1_10x_aggr_on_ATACseq/ATAC_RNA_combined_10x_aggr_labeled.RData")
# prepare the data ####
# load labeled RNAseq objects
load("/home/cli/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/integrated_labeled_updated.RData")

# project RNAseq cell type ident onto the ATACseq data
aggr_filtered$RNA.cell.type.ident <- NA
# find all barcodes corresponding to each cell type

for (i in levels(integrated_renamed_1@active.ident)){
  barcodes <- str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]][integrated_renamed_1@active.ident %in% i], 
                      end = -7L)
  
  aggr_filtered$RNA.cell.type.ident[aggr_filtered$trimmed.barcodes %in% barcodes] <- i
}

unique(aggr_filtered$RNA.cell.type.ident)

aggr_filtered$celltype.time.ident <- NA
aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "0hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% "NEFM+/CUX2- glut"] <- "NpCm_glut_0hr"
aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "1hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% "NEFM+/CUX2- glut"] <- "NpCm_glut_1hr"
aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "6hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% "NEFM+/CUX2- glut"] <- "NpCm_glut_6hr"

aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "0hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% "NEFM-/CUX2+ glut"] <- "NmCp_glut_0hr"
aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "1hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% "NEFM-/CUX2+ glut"] <- "NmCp_glut_1hr"
aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "6hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% "NEFM-/CUX2+ glut"] <- "NmCp_glut_6hr"

aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "0hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% c('GAD1+ GABA (less mature)', 
                                                                             'GAD1+/GAD2+ GABA',
                                                                             'GAD1+ GABA')] <- "GABA_0hr"
aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "1hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% c('GAD1+ GABA (less mature)', 
                                                                             'GAD1+/GAD2+ GABA',
                                                                             'GAD1+ GABA')] <- "GABA_1hr"
aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "6hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% c('GAD1+ GABA (less mature)', 
                                                                             'GAD1+/GAD2+ GABA',
                                                                             'GAD1+ GABA')] <- "GABA_6hr"

aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "0hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% c("NPC", "MAP2+ NPC")] <- "NPC_0hr"
aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "1hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% c("NPC", "MAP2+ NPC")] <- "NPC_1hr"
aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "6hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% c("NPC", "MAP2+ NPC")] <- "NPC_6hr"

aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "0hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% c("immature forebrain glut", 
                                                                             "forebrain NPC")] <- "forebrain_0hr"
aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "1hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% c("immature forebrain glut", 
                                                                             "forebrain NPC")] <- "forebrain_1hr"
aggr_filtered$celltype.time.ident[aggr_filtered$time.ident %in% "6hr" &
                                    aggr_filtered$RNA.cell.type.ident %in% c("immature forebrain glut", 
                                                                             "forebrain NPC")] <- "forebrain_6hr"
unique(aggr_filtered$celltype.time.ident)
aggr_filtered$celltype.time.ident[is.na(aggr_filtered$celltype.time.ident)] <- "others"
aggr_filtered$broad.celltype.time.ident <- NA
aggr_filtered$broad.celltype.time.ident[aggr_filtered$time.ident %in% "0hr" &
                                          aggr_filtered$celltype.time.ident %in% 
                                          c("NpCm_glut_0hr", "NmCp_glut_0hr")] <- "glut_0hr"
aggr_filtered$broad.celltype.time.ident[aggr_filtered$time.ident %in% "1hr" &
                                          aggr_filtered$celltype.time.ident %in% 
                                          c("NpCm_glut_1hr", "NmCp_glut_1hr")] <- "glut_1hr"
aggr_filtered$broad.celltype.time.ident[aggr_filtered$time.ident %in% "6hr" &
                                          aggr_filtered$celltype.time.ident %in% 
                                          c("NpCm_glut_6hr", "NmCp_glut_6hr")] <- "glut_6hr"


aggr_filtered$celltype.ident <- NA
aggr_filtered$celltype.ident[aggr_filtered$RNA.cell.type.ident %in% 'NEFM+/CUX2- glut'] <- "NpCm_glut"
aggr_filtered$celltype.ident[aggr_filtered$RNA.cell.type.ident %in% 'NEFM-/CUX2+ glut'] <- "NmCp_glut"
aggr_filtered$celltype.ident[aggr_filtered$RNA.cell.type.ident %in% c('GAD1+ GABA (less mature)', 
                                                                      'GAD1+/GAD2+ GABA',
                                                                      'GAD1+ GABA')] <- "GABA"
aggr_filtered$celltype.ident[aggr_filtered$RNA.cell.type.ident %in% c('NPC', 'MAP2+ NPC')] <- "NPC"
aggr_filtered$celltype.ident[aggr_filtered$RNA.cell.type.ident %in% c('immature forebrain glut', 
                                                                      'forebrain NPC')] <- "forebrain"

unique(aggr_filtered$celltype.ident)



# get peak set files for each cell type * time combination ####

# 14 in total: 3 cell types (NpCm glut, NmCp glut, GABA) * 3 time points +
#1 general glut cell type * 3 time points + 2 NPC cell types * 0hr time point

used_idents <- c("NpCm_glut_0hr",
                 "NpCm_glut_1hr",
                 "NpCm_glut_6hr",
                 "NmCp_glut_0hr",
                 "NmCp_glut_1hr",
                 "NmCp_glut_6hr",
                 "GABA_6hr",
                 "GABA_1hr",
                 "GABA_0hr",
                 "NPC_0hr",
                 "NPC_1hr",
                 "NPC_6hr")
# first test with one object
test_sub <- subset(x = aggr_filtered,
                   subset = celltype.time.ident == "NPC_0hr")
test_peaks <- CallPeaks(test_sub, 
                        macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2")
length(test_sub$orig.ident)

df <- data.frame(seqnames = seqnames(test_peaks),
                 starts = start(test_peaks) - 1,
                 ends = end(test_peaks),
                 strands = strand(test_peaks))
df_test <- df[grep("chr", df$seqnames), ]

write.table(df_test, 
            file="./cell_type_and_time_point_specific_peaks/test_NPC.bed", 
            quote=F, sep="\t", row.names=F, col.names=F)

# test with just 2 rows
test_two <- subset (x = aggr_filtered,
                    subset = trimmed.barcodes %in% c("AAACAGCCAAAGCGGC", "AAACAGCCAACGTGCT"))
length(test_two$orig.ident) #2

test_peaks_two <- CallPeaks(test_two, 
                            macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2")

df_two <- data.frame(seqnames = seqnames(test_peaks_two),
                     starts = start(test_peaks_two) - 1,
                     ends = end(test_peaks_two),
                     strands = strand(test_peaks_two))
df_two <- df[grep("chr", df$seqnames), ]

write.table(df_two, 
            file="./cell_type_and_time_point_specific_peaks/foo.bed", 
            quote=F, sep="\t", row.names=F, col.names=F)



length_list <- vector(mode = 'list',
                      length = 11L)
id_list <- vector(mode = 'list',
                  length = 11L)

# # loop over all needed files-part 1: three subtypes + two NPCs 
# for (i in 1:length(used_idents)){
#   # which cell type-time combination is being processed
#   id <- used_idents[i]
#   print(id)
#   
#   # create a subsetted, separate object
#   obj <- subset(x = aggr_filtered,
#                 subset = celltype.time.ident == id)
#   print(length(obj$celltype.time.ident))
#   length_list[i] <- length(obj$celltype.time.ident)
#   id_list[i] <- id
#   # call peaks with MACS2
#   obj_peaks <- CallPeaks(obj, 
#                          macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2")
#   df <- data.frame(seqnames = seqnames(obj_peaks),
#                    starts = start(obj_peaks) - 1,
#                    ends = end(obj_peaks),
#                    strands = strand(obj_peaks))
#   file_name <- paste0("./cell_type_and_time_point_specific_peaks/", 
#                       id, 
#                       "_specific_peaks.bed")
#   write.table(df, 
#               file = file_name, 
#               quote = F, 
#               sep = "\t", 
#               row.names = F, col.names = F)
#   
# }
# 
# # loop over the three files for glut combined
# used_idents <- c("glut_0hr",
#                  "glut_1hr",
#                  "glut_6hr")
# 
# for (i in 1:length(used_idents)){
#   # which cell type-time combination is being processed
#   id <- used_idents[i]
#   print(id)
#   
#   # create a subsetted, separate object
#   obj <- subset(x = aggr_filtered,
#                 subset = broad.celltype.time.ident == id)
#   print(length(obj$celltype.time.ident))
#   length_list[i] <- length(obj$celltype.time.ident)
#   id_list[i] <- id
#   
#   # call peaks with MACS2
#   obj_peaks <- CallPeaks(obj, 
#                          macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2")
#   df <- data.frame(seqnames = seqnames(obj_peaks),
#                    starts = start(obj_peaks) - 1,
#                    ends = end(obj_peaks),
#                    strands = strand(obj_peaks))
#   file_name <- paste0("./cell_type_and_time_point_specific_peaks/", 
#                       unique(peaks_by_celltype_time[[i]]$ident), 
#                       "_specific_peaks.bed")
#   write.table(df, 
#               file = file_name, 
#               quote = F, 
#               sep = "\t", 
#               row.names = F, col.names = F)
#   
# }

# call peaks on specific cell types and time points ####
# 1. call peaks specific to cell type and time point, without merging
peaks_by_celltype_time <- CallPeaks(
  object = aggr_filtered,
  group.by = "celltype.time.ident",
  macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2",
  idents = used_idents,
  combine.peaks = F
)

for (i in 1:length(peaks_by_celltype_time)){
  df <- data.frame(seqnames = seqnames(peaks_by_celltype_time[[i]]),
                   starts = start(peaks_by_celltype_time[[i]]) - 1,
                   ends = end(peaks_by_celltype_time[[i]]),
                   name = peaks_by_celltype_time[[i]]@elementMetadata@listData$name,
                   score = peaks_by_celltype_time[[i]]@elementMetadata@listData$neg_log10qvalue_summit,
                   strands = strand(peaks_by_celltype_time[[i]]))
  df <- df[grep("chr", df$seqnames), ]
  file_name <- paste0(
                      unique(peaks_by_celltype_time[[i]]$ident),
                      "_specific_peaks.bed")
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = "\t", 
              row.names = F, col.names = F)
  
}

 ordered_idents <- c(
  "NpCm_glut_1hr",
  "NpCm_glut_0hr",
  "NmCp_glut_1hr",
  "NpCm_glut_6hr",
  "GABA_6hr",
  "forebrain_0hr",
  "GABA_1hr",
  "GABA_0hr",
  "NmCp_glut_0hr",
  "NmCp_glut_6hr",
  "NPC_0hr")

for (i in 1:length(peaks_by_celltype_time)){
  id <- ordered_idents[i]
  file_name <- paste0(id, "_width_distribution_histogram.pdf")
  pdf(file = file_name)
  p <- hist(peaks_by_celltype_time[[i]]@ranges@width, 
            breaks = 200, 
            xlim = c(0, 2500),
            main = "", 
            xlab = "Peak width", 
            col = "grey40")
  print(p)
  dev.off()
}

hist(peaks_by_celltype_time[[1]]@elementMetadata@listData$score, 
     breaks = 200, 
     xlim = c(0, 1000),
     main = "", 
     xlab = "Peak intensity", 
     col = "grey40")

for (i in 1:length(peaks_by_celltype_time)){
  id <- ordered_idents[i]
  file_name <- paste0("../new_peak_set_plots/intensity_dist_histograms/", 
                      id, 
                      "_intensity_distribution_histogram.pdf")
  pdf(file = file_name)
  q <- hist(peaks_by_celltype_time[[i]]@elementMetadata@listData$score, 
            breaks = 200, 
            xlim = c(0, 1000),
            main = "", 
            xlab = "Peak intensity", 
            col = "grey40")
  print(q)
  dev.off()
}


# 2. call peaks specific for cell type only

used_idents <- c("NpCm_glut",
                 "NmCp_glut",
                 "GABA",
                 "NPC",
                 "forebrain")

peaks_by_celltype <- CallPeaks(
  object = aggr_filtered,
  group.by = "celltype.ident",
  macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2",
  idents = used_idents,
  combine.peaks = F
)

for (i in 1:length(peaks_by_celltype)){
  df <- data.frame(seqnames = seqnames(peaks_by_celltype[[i]]),
                   starts = start(peaks_by_celltype[[i]]) - 1,
                   ends = end(peaks_by_celltype[[i]]),
                   name = peaks_by_celltype[[i]]@elementMetadata@listData$name,
                   score = peaks_by_celltype[[i]]@elementMetadata@listData$neg_log10qvalue_summit,
                   strands = strand(peaks_by_celltype[[i]]))
  df <- df[grep("chr", df$seqnames), ]
  file_name <- paste0("../most_updated_obj_celltype_specific_peaks/", # pwd: /codes
                      unique(peaks_by_celltype[[i]]$ident),
                      "_specific_peaks.bed")
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = "\t", 
              row.names = F, col.names = F)
  
}


# 3. call peaks specific for time point only

used_idents <- c("0hr", "1hr", "6hr")

peaks_by_time <- CallPeaks(
  object = aggr_filtered,
  group.by = "time.ident",
  macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2",
  idents = used_idents,
  combine.peaks = F
)

for (i in 1:length(peaks_by_time)){
  df <- data.frame(seqnames = seqnames(peaks_by_time[[i]]),
                   starts = start(peaks_by_time[[i]]) - 1,
                   ends = end(peaks_by_time[[i]]),
                   name = peaks_by_time[[i]]@elementMetadata@listData$name,
                   score = peaks_by_time[[i]]@elementMetadata@listData$neg_log10qvalue_summit,
                   strands = strand(peaks_by_time[[i]]))
  df <- df[grep("chr", df$seqnames), ]
  file_name <- paste0("../most_updated_obj_time_point_specific_peaks/", # pwd: /codes
                      unique(peaks_by_time[[i]]$ident),
                      "_specific_peaks.bed")
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = "\t", 
              row.names = F, col.names = F)
  
}


# 4. merge the type&time specific peaks and use it for downstream analysis
peaks_combined <- CallPeaks(
  object = aggr_filtered,
  group.by = "celltype.time.ident",
  macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2",
  idents = used_idents,
  combine.peaks = T
)

# using the new set of peaks for downstream analysis (no subsetting)
peaks_combined <- keepStandardChromosomes(peaks_combined, 
                                          pruning.mode = "coarse")
peaks_combined <- subsetByOverlaps(x = peaks_combined, 
                                   ranges = blacklist_hg38_unified, 
                                   invert = TRUE)

#save.image(file = "1Dec2021_called_peaks_on_new_obj.RData")



df <- data.frame(seqnames = seqnames(obj_complete@assays$peaks@ranges),
                 starts = start(obj_complete@assays$peaks@ranges) - 1,
                 ends = end(obj_complete@assays$peaks@ranges),
                 strands = strand(obj_complete@assays$peaks@ranges))
file_name <- "merged_union_peak_set.bed"
write.table(df, 
            file = file_name, 
            quote = F, 
            sep = "\t", 
            row.names = F, col.names = F)


# quantify counts in each peak
macs2_counts_complete <- FeatureMatrix(
  fragments = Fragments(aggr_filtered),
  features = peaks_combined,
  cells = colnames(aggr_filtered)
)

fragpath <- "/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/atac_fragments.tsv.gz"
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/annot_aggr_signac_EnsDb_UCSC.RData")

#sum(!(aggr_filtered$celltype.time.ident %in% new_obj_celltypes)) #26064 out of 82516


obj_complete <- aggr_filtered
# create a new assay using the MACS2 peak set and add it to the Seurat object
obj_complete[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts_complete,
  fragments = fragpath,
  annotation = annot_aggr_signac_ucsc
)

# calcualate percent counts in peaks
total_fragments <- CountFragments('/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/atac_fragments.tsv.gz')

# add total fragments to the object

# test
barcodes <- obj_complete@assays$ATAC@counts@Dimnames[[2]]
obj_complete$fragments <- total_fragments$frequency_count[total_fragments$CB %in% barcodes]

# calculate percent fragments in peaks 
obj_complete <- FRiP(
  object = obj_complete,
  assay = 'ATAC',
  total.fragments = 'fragments'
)

# calculate counts of fragments in peaks 
hist(macs2_counts_complete@i)
hist(macs2_counts_complete@x, breaks = 1000, xlim = c(0, 15))
obj_complete$in_peaks_fraction <- NA

length(macs2_counts_complete@x)
in_peaks_count <- colSums(macs2_counts_complete)
head(obj_complete@assays$peaks@counts@Dimnames[[2]])
head(in_peaks_count)

recounted_fragments <- colSums(obj_complete@assays$ATAC@counts)
max(recounted_fragments)
max(in_peaks_count)

obj_complete$in_peaks_count <- in_peaks_count
obj_complete$in_peaks_fraction <- obj_complete$in_peaks_count / recounted_fragments * 100
hist(obj_complete$in_peaks_fraction, breaks = 1000)

obj_complete$fragments <- NULL
obj_complete$FRiP <- NULL

# count fraction of fragments in peaks for each cell type_time combination
for (i in used_idents){
  sub_lst <- obj_complete$in_peaks_fraction[obj_complete$celltype.time.ident %in% i]
  print(paste0(i, ": ", sum(sub_lst)/length(sub_lst)))
}



# obj_complete$in_peaks_fraction <- FractionCountsInRegion(
#   object = obj_complete, 
#   assay = 'peaks',
#   regions = peaks_combined
# )
# hist(obj_complete$in_peaks_fraction)

# downstream analysis: dimensional reduction and clustering ####
obj_complete <- FindTopFeatures(obj_complete, min.cutoff = 5)
obj_complete <- RunTFIDF(obj_complete)
obj_complete <- RunSVD(obj_complete)

# check which components are associated with technical (nonbiological) variance (sequencing depth)
DepthCor(obj_complete)

# clustering and downstream analysis
obj_complete <- RunUMAP(object = obj_complete, 
                        reduction = 'lsi', 
                        dims = 2:30, # remove the first component
                        seed.use = 999) 

obj_complete <- FindNeighbors(object = obj_complete, 
                              reduction = 'lsi', 
                              dims = 2:30) # remove the first component

obj_complete <- FindClusters(object = obj_complete, 
                             verbose = FALSE, 
                             algorithm = 3, 
                             resolution = 0.5,
                             random.seed = 99)

Idents(obj_complete) <- "seurat_clusters"

DimPlot(object = obj_complete, label = F)


DefaultAssay(obj_complete) <- "RNA"

# check important marker gene distributions
FeaturePlot(obj_complete, 
            features = c("GAD1"), 
            max.cutoff = 7.5)

FeaturePlot(obj_complete, 
            features = c("SLC17A6"), 
            max.cutoff = 10)

FeaturePlot(obj_complete, 
            features = c("VIM"), 
            max.cutoff = 20)

FeaturePlot(obj_complete, 
            features = c("NES"), 
            max.cutoff = 15)


# plot projection of RNAseq onto ATACseq
DefaultAssay(obj_complete) <- "ATAC"

DimPlot(object = obj_complete, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.cluster.ident") + 
  ggtitle("ATACseq data \n projected by RNAseq cluster numbers")

DimPlot(object = obj_complete, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.cell.type.ident") + 
  ggtitle("ATACseq data \n projected by RNAseq cell types") +
  NoLegend()

DimPlot(object = obj_complete, 
        label = TRUE, 
        repel = T,
        group.by = "fine.cell.type") + 
  ggtitle("ATACseq data \n projected by RNAseq cell types") +
  NoLegend()
# check library composition
DefaultAssay(obj_complete) <- 'peaks'

DimPlot(obj_complete, group.by = 'time.group.ident') + 
  ggtitle("distribution of group_time identity")

# project ATACseq onto RNAseq
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/codes/integrated_with_all_labels.RData")
integrated_renamed_1$ATAC.cluster.ident <- NA
integrated_renamed_1$trimmed.barcodes <- str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]], end = -7L)

for (i in 1:length(levels(obj_complete$seurat_clusters))){
  barcodes <- str_sub(obj_complete@assays$peaks@counts@Dimnames[[2]][obj_complete$seurat_clusters %in% 
                                                                       as.character(i)], end = -3L)
  
  integrated_renamed_1$ATAC.cluster.ident[integrated_renamed_1$trimmed.barcodes 
                                          %in% barcodes] <- as.character(i) 
}

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = integrated_renamed_1, 
        label = TRUE, 
        repel = T,
        group.by = "ATAC.cluster.ident") + 
  ggtitle("RNAseq data \n projected by ATACseq cluster numbers")


# plot ATACseq data by time
dimplot_by_time <- DimPlot(object = obj_complete,
                           group.by = "time.ident",
                           cols = c("magenta", "yellow", "cyan")) + 
  ggtitle("ATACseq data \n colored by time points")
dimplot_by_time[[1]]$layers[[1]]$aes_params$alpha = .2
dimplot_by_time


DimPlot(object = obj_complete, 
        group.by = "time.ident",
        cols = c("magenta", "gray", "gray")) + 
  ggtitle("0hr")

DimPlot(object = obj_complete, 
        group.by = "time.ident",
        cols = c("gray", "yellow", "gray")) + 
  ggtitle("1hr")

DimPlot(object = obj_complete, 
        group.by = "time.ident",
        cols = c("gray", "gray", "cyan")) + 
  ggtitle("6hr")


# calculate the number of cells in each bigger cell type
unique(obj_complete$fine.cell.type)
subset(obj_complete, subset = fine.cell.type == "NEFM-/CUX2+ glut") #NmCp glut: 5982
subset(obj_complete, subset = fine.cell.type == "NEFM+/CUX2- glut") #NpCm glut: 19045
subset(obj_complete, subset = fine.cell.type == "GABA") #GABA: 18156
subset(obj_complete, subset = fine.cell.type == "forebrain") #forebrain NPC: 7772
subset(obj_complete, subset = fine.cell.type == "NPC") #NPC: 5498
sum(obj_complete$fine.cell.type == "SEMA3E+ glut", na.rm = T) #SEMA3Ep glut: 5732
sum(obj_complete$fine.cell.type == "immature neuron", na.rm = T) #immature neuron: 4578
sum(obj_complete$fine.cell.type == "subcerebral immature neuron", na.rm = T) #subcerebral immature neuron: 3786



DimPlot(object = obj_complete,
        group.by = "fine.cell.type",
        label = T,
        repel = T) + 
  ggtitle("ATACseq data colored by major cell types")

# plot peak tracks for early and late response genes to check

idents_to_check <- c("GABA",
                     "NEFM+/CUX2- glut",
                     "NEFM-/CUX2+ glut",
                     "forebrain",
                     "NPC")

for (i in 1:length(idents_to_check)){
  celltype <- gsub(" ", "_", idents_to_check[i])
  celltype <- gsub("/", "_", idents_to_check[i])
  file_name <- paste0("../new_peak_set_plots/peak_track_plots/BDNF_",
                      celltype,
                      "_peak_tracks.pdf")
  pdf(file = file_name)
  p <- CoveragePlot(
    object = subset(obj_complete, subset = fine.cell.type == idents_to_check[i]),
    region = "chr11-27689113-27774038", # BDNF
    region.highlight = GRanges(seqnames = "chr11",
                               ranges = IRanges(start = 27767000, end = 27777000)),
    group.by = "time.ident"
  )
  print(p)
  dev.off()
}

for (i in 1:length(idents_to_check)){
  celltype <- gsub(" ", "_", idents_to_check[i])
  celltype <- gsub("/", "_", idents_to_check[i])
  file_name <- paste0("../new_peak_set_plots/peak_track_plots/NPAS4_",
                      celltype,
                      "_peak_tracks.pdf")
  pdf(file = file_name)
  p <- CoveragePlot(
    object = subset(obj_complete, subset = fine.cell.type == idents_to_check[i]),
    region = "chr11-66414647-66427411", # NPAS4
    region.highlight = GRanges(seqnames = "chr11",
                               ranges = IRanges(start = 66418000, end = 66418700)),
    group.by = "time.ident"
  )
  
  print(p)
  dev.off()
}

CoveragePlot(
  object = subset(obj_complete, subset = fine.cell.type == "GABA"),
  region = "VGF",
  features = "VGF",
  extend.upstream = 10000,
  extend.downstream = 10000,  
  group.by = "time.ident"
)

CoveragePlot(
  object = subset(obj_complete, subset = fine.cell.type == "NEFM+/CUX2- glut"),
  region = "VGF",
  features = "VGF",
  extend.upstream = 10000,
  extend.downstream = 10000,  
  group.by = "time.ident"
)

CoveragePlot(
  object = subset(obj_complete, subset = fine.cell.type == "NEFM-/CUX2+ glut"),
  region = "VGF",
  features = "VGF",
  extend.upstream = 10000,
  extend.downstream = 10000,  
  group.by = "time.ident"
)

# # subset the peaks and the object for downstream analysis
# new_obj_celltypes <- c("NpCm_glut_0hr",
#                        "NpCm_glut_1hr",
#                        "NpCm_glut_6hr",
#                        "NmCp_glut_0hr",
#                        "NmCp_glut_1hr",
#                        "NmCp_glut_6hr",
#                        "GABA_6hr",
#                        "GABA_1hr",
#                        "GABA_0hr",
#                        "forebrain_0hr",
#                        "forebrain_1hr",
#                        "forebrain_6hr",
#                        "NPC_6hr",
#                        "NPC_1hr",
#                        "NPC_0hr")
# 
# peaks_four_celltypes <- peaks_combined[peaks_combined$peak_called_in %in%
#                                          new_obj_celltypes]
# 
# 
# # remove peaks on nonstandard chromosomes and in genomic blacklist regions
# peaks_four_celltypes <- keepStandardChromosomes(peaks_four_celltypes, 
#                                                 pruning.mode = "coarse")
# peaks_four_celltypes <- subsetByOverlaps(x = peaks_four_celltypes, 
#                                          ranges = blacklist_hg38_unified, 
#                                          invert = TRUE)
# 
# obj_four_celltypes <- subset(obj_complete, 
#                              subset = celltype.time.ident %in% new_obj_celltypes)
# # quantify counts in each peak
# macs2_counts_four_celltypes <- FeatureMatrix(
#   fragments = Fragments(obj_four_celltypes),
#   features = peaks_four_celltypes,
#   cells = colnames(obj_four_celltypes)
# )
# 
# fragpath <- "/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/atac_fragments.tsv.gz"
# 
# 
# # create a new assay using the MACS2 peak set and add it to the Seurat object
# obj_four_celltypes[["peaks"]] <- 
#   
#   CreateChromatinAssay(
#   counts = macs2_counts,
#   fragments = fragpath,
#   annotation = annot_aggr_signac_ucsc
# )

# # downstream analysis: dimensional reduction and clustering
# obj_four_celltypes <- FindTopFeatures(obj_four_celltypes, min.cutoff = 5)
# obj_four_celltypes <- RunTFIDF(obj_four_celltypes)
# obj_four_celltypes <- RunSVD(obj_four_celltypes)
# 
# #save.image("Nov292021_called_new_peaks.RData")
# 
# 
# # check which components are associated with technical (nonbiological) variance (sequencing depth)
# DepthCor(obj_four_celltypes)
# 
# # clustering and downstream analysis
# obj_four_celltypes <- RunUMAP(object = obj_four_celltypes, 
#                               reduction = 'lsi', 
#                               dims = 2:30, # remove the first component
#                               seed.use = 999) 
# 
# obj_four_celltypes <- FindNeighbors(object = obj_four_celltypes, 
#                                     reduction = 'lsi', 
#                                     dims = 2:30) # remove the first component
# 
# obj_four_celltypes <- FindClusters(object = obj_four_celltypes, 
#                                    verbose = FALSE, 
#                                    algorithm = 3, 
#                                    resolution = 0.5,
#                                    random.seed = 99)
# 
# DimPlot(object = obj_four_celltypes, label = TRUE) + 
#   NoLegend() +
#   ggtitle("clustering of ATACseq data with new peak set")
# 
# CoveragePlot(
#   object = subset(obj_four_celltypes),
#   group.by = "time.ident",
#   region = "chr11-27719113-27734038" # BDNF
# )


# Differentially accessible peaks by time ####
ident_list <- c("GABA_1hr", "GABA_0hr",
                "GABA_6hr", "GABA_0hr",
                "GABA_6hr", "GABA_1hr",
                
                "NmCp_glut_1hr", "NmCp_glut_0hr", 
                "NmCp_glut_6hr", "NmCp_glut_0hr", 
                "NmCp_glut_6hr", "NmCp_glut_1hr",
                
                "NpCm_glut_1hr", "NpCm_glut_0hr", 
                "NpCm_glut_6hr", "NpCm_glut_0hr",
                "NpCm_glut_6hr", "NpCm_glut_1hr"
)

da_peaks_by_time_list <- vector(mode = "list",
                                length = (0.5*length(ident_list)))

Idents(obj_complete) <- "celltype.time.ident"
for (i in 1:(0.5*length(ident_list))){
  # get the correct indices
  ind1 <- 2 * i - 1
  ident1 <- ident_list[ind1]
  ind2 <- 2 * i
  ident2 <- ident_list[ind2]
  
  print(paste0("baseline: ", ident1, " comparison: ", ident2))
  
  # call FindMarkers to find differentially accessible peaks
  da_peaks <- FindMarkers(
    object = obj_complete,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = 'wilcox',
    logfc.threshold = 0,
    min.pct = 0.05)
  
  da_peaks_by_time_list[[i]] <- da_peaks
  
  # write the output into csv files for record
  df <- data.frame(peaks = rownames(da_peaks),
                   p_val = da_peaks$p_val,
                   avg_log2FC = da_peaks$avg_log2FC,
                   pct.1 = da_peaks$pct.1,
                   pct.2 = da_peaks$pct.2,
                   p_val_adj = da_peaks$p_val_adj)
  file_name <- paste0("../most_updated_obj_da_peaks_by_time_unfiltered/", # pwd: /codes
                      ident1, "_",
                      ident2, ".csv"
  )
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = ",", 
              row.names = F, 
              col.names = F)
}

da_peaks_GABA_1v0hr_test <- FindMarkers(
  object = obj_complete,
  ident.1 = "GABA_1hr", 
  ident.2 = "GABA_0hr",
  test.use = 'wilcox',
  logfc.threshold = 0.25,
  min.pct = 0,
  verbose = T)

da_peaks_GABA_1v0hr_test_2 <- FindMarkers(
  object = obj_complete,
  ident.1 = "GABA_1hr", 
  ident.2 = "GABA_0hr",
  test.use = 'wilcox',
  logfc.threshold = 0,
  min.pct = 0)

# using LR test
da_peaks_GABA_1v0hr_LR <- FindMarkers(
  object = obj_complete,
  ident.1 = "GABA_1hr", 
  ident.2 = "GABA_0hr",
  test.use = 'LR', 
  latent.vars = "in_peaks_count",
  logfc.threshold = 0,
  min.pct = 0)

da_peaks_NmCp_glut_1v0hr_LR <- FindMarkers(
  object = obj_complete,
  ident.1 = "NmCp_glut_1hr", 
  ident.2 = "NmCp_glut_0hr",
  test.use = 'LR', 
  latent.vars = "in_peaks_count",
  logfc.threshold = 0,
  min.pct = 0)

da_peaks_NpCm_glut_1v0hr_LR <- FindMarkers(
  object = obj_complete,
  ident.1 = "NpCm_glut_1hr", 
  ident.2 = "NpCm_glut_0hr",
  test.use = 'LR', 
  latent.vars = "in_peaks_count",
  logfc.threshold = 0,
  min.pct = 0)

save.image(file = "16Dec2021_glutpeaks_with_LR.RData")

df_NmCp_LR <- data.frame(p_val = da_peaks_NmCp_glut_1v0hr_LR$p_val,
                         avg_log2FC = da_peaks_NmCp_glut_1v0hr_LR$avg_log2FC,
                         pct.1 = da_peaks_NmCp_glut_1v0hr_LR$pct.1,
                         pct.2 = da_peaks_NmCp_glut_1v0hr_LR$pct.2,
                         p_val_adj = da_peaks_NmCp_glut_1v0hr_LR$p_val_adj, 
                         peak <- rownames(da_peaks_NmCp_glut_1v0hr_LR),
                         stringsAsFactors = T)

df_NpCm_LR <- data.frame(p_val = da_peaks_NpCm_glut_1v0hr_LR$p_val,
                         avg_log2FC = da_peaks_NpCm_glut_1v0hr_LR$avg_log2FC,
                         pct.1 = da_peaks_NpCm_glut_1v0hr_LR$pct.1,
                         pct.2 = da_peaks_NpCm_glut_1v0hr_LR$pct.2,
                         p_val_adj = da_peaks_NpCm_glut_1v0hr_LR$p_val_adj, 
                         peak <- rownames(da_peaks_NpCm_glut_1v0hr_LR),
                         stringsAsFactors = T)
write.table(df_NmCp_LR, 
            file = "NmCp_glut_LR_all_peaks.csv",
            quote = F, 
            sep = ",", 
            row.names = F, 
            col.names = F)
write.table(df_NpCm_LR, 
            file = "NpCm_glut_LR_all_peaks.csv",
            quote = F, 
            sep = ",", 
            row.names = F, 
            col.names = F)
df_GABA_LR <- data.frame(p_val = da_peaks_GABA_1v0hr_test_2$p_val,
                         avg_log2FC = da_peaks_GABA_1v0hr_test_2$avg_log2FC,
                         pct.1 = da_peaks_GABA_1v0hr_test_2$pct.1,
                         pct.2 = da_peaks_GABA_1v0hr_test_2$pct.2,
                         p_val_adj = da_peaks_GABA_1v0hr_test_2$p_val_adj, 
                         stringsAsFactors = T)
write.table(df_GABA_LR, 
            file = "GABA_LR_all_peaks.csv",
            quote = F, 
            sep = ",", 
            row.names = F, 
            col.names = F)

# reduce the thresholds to 0 and run again for all cell types
for (i in 1:(0.5*length(ident_list))){
  # get the correct indices
  ind1 <- 2 * i - 1
  ident1 <- ident_list[ind1]
  ind2 <- 2 * i
  ident2 <- ident_list[ind2]
  
  print(paste0("baseline: ", ident1, " comparison: ", ident2))
  
  # call FindMarkers to find differentially accessible peaks
  da_peaks <- FindMarkers(
    object = obj_complete,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = 'wilcox',
    logfc.threshold = 0,
    min.pct = 0)
  
  da_peaks_by_time_list[[i]] <- da_peaks
  
  # write the output into csv files for record
  df <- data.frame(peaks = rownames(da_peaks),
                   p_val = da_peaks$p_val,
                   avg_log2FC = da_peaks$avg_log2FC,
                   pct.1 = da_peaks$pct.1,
                   pct.2 = da_peaks$pct.2,
                   p_val_adj = da_peaks$p_val_adj)
  file_name <- paste0("../most_updated_obj_da_peaks_by_time_unfiltered/", # pwd: /codes
                      ident1, "_",
                      ident2, ".csv"
  )
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = ",", 
              row.names = F, 
              col.names = F)
}

# write bed files for annotation
# HOMER requires: Column1: chromosome, Column2: starting position, 
#Column3: ending position, Column4: Unique Peak ID, Column5: not used, 
#Column6: Strand (+/- or 0/1, where 0="+", 1="-")

for (i in 1:(0.5*length(ident_list))){
  # get the correct indices
  ind1 <- 2 * i - 1
  ident1 <- ident_list[ind1]
  ind2 <- 2 * i
  ident2 <- ident_list[ind2]
  
  print(paste0("baseline: ", ident1, " comparison: ", ident2))
  
  da_peaks <- da_peaks_by_time_list[[i]]
  
  # write the output into bed files for annotation
  df <- data.frame(chromosome = gsub(pattern = "-",
                                     replacement = "",
                                     x = str_sub(rownames(da_peaks), end = 5)),
                   start = gsub(pattern = "-",
                                replacement = "",
                                x = str_extract_all(string = rownames(da_peaks),
                                                    pattern = "-[0-9]+-"
                                )),
                   end = gsub(pattern = "-",
                              replacement = "",
                              x = str_extract_all(string = rownames(da_peaks),
                                                  pattern = "-[0-9]+$")),
                   peakID = paste0(da_peaks$avg_log2FC, '^', da_peaks$p_val, '^', da_peaks$p_val_adj),
                   p_val = da_peaks$p_val,
                   strand = rep("+", times = length(rownames(da_peaks))))
  file_name <- paste0("../most_updated_obj_da_peaks_by_time/bed/", # pwd: /codes
                      ident1, "_",
                      ident2, ".bed")
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = "\t", 
              row.names = F, 
              col.names = F)
}

# da_peaks_GABA_0v1hr <- FindMarkers(
#   object = obj_complete,
#   ident.1 = 'GABA_0hr',
#   ident.2 = 'GABA_1hr',
#   test.use = 'bimod',
#   logfc.threshold = 0.1
# )


# plot volcano plots to show da peaks
library(ggplot2)

# test plot
df <- da_peaks_by_time_list[[1]][, c("p_val", "avg_log2FC", "p_val_adj")]
df$significance <- "nonsig"
df$significance[df$p_val_adj < 0.05 & df$avg_log2FC > 0] <- "pos"
df$significance[df$p_val_adj < 0.05 & df$avg_log2FC < 0] <- "neg"
df$neg_log_p_val <- (0 - log10(df$p_val))

colors <- setNames(c("steelblue3", "grey", "red3"), c("neg", "nonsig", "pos"))

ggplot(data = df,
       aes(x = avg_log2FC, 
           y = neg_log_p_val, 
           color = significance)) + 
  geom_point() +
  xlim(-0.5, 0.5) +
  scale_color_manual(values = colors) +
  theme_minimal()

write.table(df, 
            file = "used_to_plot_volcano_plot.csv", 
            quote = F, 
            sep = ",", 
            row.names = F, 
            col.names = F)

for (i in 1:(0.5 * length(ident_list))){
  # get the correct indices
  ind1 <- 2 * i - 1
  ident1 <- ident_list[ind1]
  ind2 <- 2 * i
  ident2 <- ident_list[ind2]
  
  df <- da_peaks_by_time_list[[i]][, c("p_val", "avg_log2FC", "p_val_adj")]
  df$significance <- "nonsig"
  df$significance[df$p_val_adj < 0.05 & df$avg_log2FC > 0] <- "pos"
  df$significance[df$p_val_adj < 0.05 & df$avg_log2FC < 0] <- "neg"
  df$neg_log_p_val <- (0 - log10(df$p_val))
  
  file_name <- paste0("../new_peak_set_plots/da_peaks_by_time_volcano_plots/", # pwd: /codes
                      ident1, "_",
                      ident2, "_da_peaks_volcano_plot.pdf"
  )
  
  pdf(file = file_name)
  p <- ggplot(data = df,
              aes(x = avg_log2FC, 
                  y = neg_log_p_val, 
                  color = significance)) + 
    geom_point() +
    scale_color_manual(values = colors) +
    theme_minimal()
  print(p)
  dev.off()
}

# check the distribution of logFC and pval of positive and negative logFC da peaks
hist(df$neg_log_p_val, 
     breaks = 500, 
     xlim = c(0, 35),
     main = "GABA 0v1hr -log(p) distribution")

hist(df$avg_log2FC, 
     breaks = 500, 
     main = "GABA 0v1hr average logFC distribution")
hist(df$avg_log2FC[df$avg_log2FC < 0], 
     breaks = 500)

hist(df$neg_log_p_val[df$avg_log2FC > 0], 
     breaks = 500,
     xlim = c(0, 60))

hist(df$neg_log_p_val[df$avg_log2FC < 0], 
     breaks = 500)

colors <- setNames(c("steelblue3", "grey", "red3"), c("neg", "nonsig", "pos"))

ggplot(data = df,
       aes(x = avg_log2FC, 
           y = neg_log_p_val, 
           color = significance)) + 
  geom_point(size = 0.3, alpha = 0.5) +
  xlim(-0.5, 0.5) +
  scale_color_manual(values = colors) +
  theme_minimal()


# Differentially accessible peaks by cell type ####
ident_list <- c("GABA_0hr", "NmCp_glut_0hr",
                "GABA_1hr", "NmCp_glut_1hr",
                "GABA_6hr", "NmCp_glut_6hr",
                
                "GABA_0hr", "NpCm_glut_0hr",
                "GABA_1hr", "NpCm_glut_1hr",
                "GABA_6hr", "NpCm_glut_6hr",
                
                "NmCp_glut_0hr", "NpCm_glut_0hr",
                "NmCp_glut_1hr", "NpCm_glut_1hr",
                "NmCp_glut_6hr", "NpCm_glut_6hr")

da_peaks_by_celltype_list <- vector(mode = "list",
                                    length = (length(ident_list) / 2))

for (i in 1:length(da_peaks_by_celltype_list)){
  # get the correct indices
  ind1 <- 2 * i - 1
  ident2 <- ident_list[ind1] # reversed because pos log2FC = peak is more open in the first group
  ind2 <- 2 * i
  ident1 <- ident_list[ind2]
  
  # call FindMarkers to find differentially accessible peaks
  da_peaks <- FindMarkers(
    object = obj_complete,
    ident.1 = ident1, 
    ident.2 = ident2,
    test.use = 'wilcox',
    logfc.threshold = 0.02,
    min.pct = 0.05)
  
  da_peaks_by_celltype_list[[i]] <- da_peaks
  
  # write the output into csv files
  df <- data.frame(peaks = rownames(da_peaks),
                   p_val = da_peaks$p_val,
                   avg_log2FC = da_peaks$avg_log2FC,
                   pct.1 = da_peaks$pct.1,
                   pct.2 = da_peaks$pct.2,
                   p_val_adj = da_peaks$p_val_adj)
  file_name <- paste0("../most_updated_obj_da_peaks_by_celltype/", # pwd: /codes
                      ident1, "_",
                      ident2, ".csv"
  )
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = ",", 
              row.names = F, 
              col.names = F)
  
}

# da_peaks_GABA_0v1hr <- FindMarkers(
#   object = obj_complete,
#   ident.1 = 'GABA_0hr',
#   ident.2 = 'GABA_1hr',
#   test.use = 'bimod',
#   logfc.threshold = 0.1
# )


for (i in 1:(0.5*length(ident_list))){
  # get the correct indices
  ind1 <- 2 * i - 1
  ident2 <- ident_list[ind1]
  ind2 <- 2 * i
  ident1 <- ident_list[ind2]
  
  print(paste0("baseline: ", ident1, " comparison: ", ident2))
  
  da_peaks <- da_peaks_by_celltype_list[[i]]
  
  # write the output into bed files for annotation
  df <- data.frame(chromosome = gsub(pattern = "-",
                                     replacement = "",
                                     x = str_sub(rownames(da_peaks), end = 5)),
                   start = gsub(pattern = "-",
                                replacement = "",
                                x = str_extract_all(string = rownames(da_peaks),
                                                    pattern = "-[0-9]+-"
                                )),
                   end = gsub(pattern = "-",
                              replacement = "",
                              x = str_extract_all(string = rownames(da_peaks),
                                                  pattern = "-[0-9]+$")),
                   id = paste0(da_peaks$avg_log2FC, '^', da_peaks$p_val, '^', da_peaks$p_val_adj),
                   p_val = da_peaks$p_val,
                   strand = rep("+", times = length(rownames(da_peaks))))
  file_name <- paste0("../most_updated_obj_da_peaks_by_celltype/bed/", # pwd: /codes
                      ident1, "_",
                      ident2, ".bed")
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = "\t", 
              row.names = F, 
              col.names = F)
}

# plot volcano plots to show da peaks
library(ggplot2)

for (i in 1:length(da_peaks_by_celltype_list)){
  # get the correct indices
  ind1 <- 2 * i - 1
  ident1 <- ident_list[ind1]
  ind2 <- 2 * i
  ident2 <- ident_list[ind2]
  
  # call FindMarkers to find differentially accessible peaks
  df <- da_peaks_by_celltype_list[[i]][, c("p_val", "avg_log2FC", "p_val_adj")]
  df$significance <- "nonsig"
  df$significance[df$p_val_adj < 0.05 & df$avg_log2FC > 0] <- "pos"
  df$significance[df$p_val_adj < 0.05 & df$avg_log2FC < 0] <- "neg"
  df$neg_log_p_val <- (0 - log10(df$p_val))
  
  file_name <- paste0("../new_peak_set_plots/da_peaks_by_celltype_volcano_plots/", # pwd: /codes
                      ident1, "_",
                      ident2, "_da_peaks_volcano_plot.pdf"
  )
  
  pdf(file = file_name)
  p <- ggplot(data = df,
              aes(x = avg_log2FC, 
                  y = neg_log_p_val, 
                  color = significance)) + 
    geom_point() +
    scale_color_manual(values = colors) +
    theme_minimal()
  print(p)
  dev.off()
}

save.image(file = "02Dec2021_calc_da_peaks_motif.RData")


# Motif enrichment ####

# add motif information
library(JASPAR2020)
library(TFBSTools)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", 
              tax_group = 'vertebrates', 
              species = 9606,
              all_versions = FALSE)
)
obj_complete <- AddMotifs(
  object = obj_complete,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# separate the peaks into positive and negative change
da_peaks_by_time_pos <- vector(mode = "list", length = length(da_peaks_by_time_list))
da_peaks_by_time_neg <- vector(mode = "list", length = length(da_peaks_by_time_list))

for (i in 1:length(da_peaks_by_time_list)){
  da_peaks_by_time_pos[[i]] <- da_peaks_by_time_list[[i]][da_peaks_by_time_list[[i]]$avg_log2FC > 0, ]
  da_peaks_by_time_neg[[i]] <- da_peaks_by_time_list[[i]][da_peaks_by_time_list[[i]]$avg_log2FC < 0, ]
}

da_peaks_by_celltype_pos <- vector(mode = "list", length = length(da_peaks_by_celltype_list))
da_peaks_by_celltype_neg <- vector(mode = "list", length = length(da_peaks_by_celltype_list))

for (i in 1:length(da_peaks_by_time_list)){
  da_peaks_by_celltype_pos[[i]] <- da_peaks_by_celltype_list[[i]][da_peaks_by_celltype_list[[i]]$avg_log2FC > 0, ]
  da_peaks_by_celltype_neg[[i]] <- da_peaks_by_celltype_list[[i]][da_peaks_by_celltype_list[[i]]$avg_log2FC < 0, ]
}

# find enrichment for each pair of time points within cell type
for (i in 1:length(da_peaks_by_time_list)){
  
  enriched_motifs <- FindMotifs(
    object = obj_complete,
    features = rownames(da_peaks_by_time_pos[[i]]))
  
  ind1 <- 2 * i - 1
  ident1 <- ident_list[ind1]
  ind2 <- 2 * i
  ident2 <- ident_list[ind2]
  
  print(paste0("passed: ", ident1, "_", ident2))
  file_name <- paste0("../new_peak_set_plots/motif_by_time_plots/",
                      ident1, "_vs_",
                      ident2, "_positive_enriched_motifs.pdf")
  pdf(file = file_name)
  p <- MotifPlot(
    object = obj_complete,
    motifs = head(rownames(enriched_motifs), n = 12)
  )
  print(p)
  dev.off()
  
  # repeat for negative FC
  enriched_motifs <- FindMotifs(
    object = obj_complete,
    features = rownames(da_peaks_by_time_neg[[i]]))
  
  print(paste0("passed: ", ident1, "_", ident2))
  file_name <- paste0("../new_peak_set_plots/motif_by_time_plots/",
                      ident1, "_vs_",
                      ident2, "_negative_enriched_motifs.pdf")
  pdf(file = file_name)
  p <- MotifPlot(
    object = obj_complete,
    motifs = head(rownames(enriched_motifs), n = 12)
  )
  print(p)
  dev.off()
}

# find enrichment for each pair of cell types
for (i in 1:length(da_peaks_by_celltype_list)){
  
  enriched_motifs <- FindMotifs(
    object = obj_complete,
    features = rownames(da_peaks_by_celltype_pos[[i]]))
  
  ind1 <- 2 * i - 1
  ident1 <- ident_list[ind1]
  ind2 <- 2 * i
  ident2 <- ident_list[ind2]
  
  print(paste0("passed: ", ident1, "_", ident2))
  file_name <- paste0("../new_peak_set_plots/motif_by_celltype_plots/",
                      ident1, "_vs_",
                      ident2, "_positive_enriched_motifs.pdf")
  pdf(file = file_name)
  p <- MotifPlot(
    object = obj_complete,
    motifs = head(rownames(enriched_motifs), n = 12)
  )
  print(p)
  dev.off()
  
  # repeat for negative FC
  enriched_motifs <- FindMotifs(
    object = obj_complete,
    features = rownames(da_peaks_by_celltype_neg[[i]]))
  
  print(paste0("passed: ", ident1, "_", ident2))
  file_name <- paste0("../new_peak_set_plots/motif_by_celltype_plots/",
                      ident1, "_vs_",
                      ident2, "_negative_enriched_motifs.pdf")
  pdf(file = file_name)
  p <- MotifPlot(
    object = obj_complete,
    motifs = head(rownames(enriched_motifs), n = 12)
  )
  print(p)
  dev.off()
}

#save.image("new_peak_set_after_motif.RData")


# link peaks to gene ####
# separate object_complete into 4 cell types
# linkage plots separating time points
obj_complete$celltype.time.ident[is.na(obj_complete$celltype.time.ident)] <- "NA"
obj_complete <- subset(obj_complete, celltype.time.ident != "NA")
idents_list <- unique(obj_complete$celltype.time.ident)
idents_list <- sort(idents_list)
idents_list <- idents_list[4:15]
celltype_time_sep_list <- vector(mode = "list", 
                                 length = length(idents_list))
for (i in 1:length(celltype_time_sep_list)){
  temp_obj <- subset(obj_complete, subset = celltype.time.ident == idents_list[i])
  
  # first compute the GC content for each peak
  DefaultAssay(temp_obj) <- "peaks"
  temp_obj <- RegionStats(temp_obj, 
                          genome = BSgenome.Hsapiens.UCSC.hg38)
  print("LinkPeaks() running")
  # link peaks to genes
  temp_obj <- LinkPeaks(
    object = temp_obj,
    peak.assay = "peaks",
    expression.assay = "SCT")
  write.table(x = as.data.frame(Links(temp_obj)), 
              file = paste0(idents_list[i], ".csv"),
              quote = F, sep = ",", col.names = T, row.names = F)
  print("LinkPeaks() finished")
  celltype_time_sep_list[[i]] <- as.data.frame(Links(temp_obj))
}

# test
CoveragePlot(
  object = celltype_time_sep_list[[1]],
  region = "BDNF",
  features = "BDNF",
  expression.assay = "SCT",
  extend.upstream = 10000,
  extend.downstream = 100000
)
CoveragePlot(
  object = celltype_time_sep_list[[3]],
  region = "BDNF",
  features = "BDNF",
  expression.assay = "SCT",
  extend.upstream = 10000,
  extend.downstream = 100000
)

CoveragePlot(
  object = celltype_time_sep_list[[4]],
  region = "NPAS4",
  features = "NPAS4",
  expression.assay = "SCT",
  extend.upstream = 10000,
  extend.downstream = 10000
)


for (i in 1:9){
  id <- ident_list[i]
  file_name <- paste0("../new_peak_set_plots/linkage/BDNF_", id, ".pdf")
  
  pdf(file_name, width = 7.25, height = 3.35)
  p <- CoveragePlot(
    object = celltype_time_sep_list[[i]],
    region = "BDNF",
    features = "BDNF",
    expression.assay = "SCT",
    extend.upstream = 10000,
    extend.downstream = 100000
  )
  print(p)
  dev.off()
}

# try with only separating the cell types
celltype_only_sep_list <- vector(mode = "list", 
                                 length = (length(unique(obj_complete$celltype.time.ident)) - 1) / 3)

ident_list <- c("NpCm_glut_0hr", "NpCm_glut_1hr", "NpCm_glut_6hr",
                "NmCp_glut_0hr", "NmCp_glut_1hr", "NmCp_glut_6hr",
                "GABA_0hr", "GABA_1hr", "GABA_6hr",
                "forebrain_0hr", "forebrain_1hr", "forebrain_6hr",
                "NPC_0hr", "NPC_1hr", "NPC_6hr")

for (i in 1:length(celltype_only_sep_list)){
  temp_obj2 <- subset(obj_complete, subset = celltype.time.ident %in% c(ident_list[3 * i - 2],
                                                                       ident_list[3 * i - 1],
                                                                       ident_list[3 * i]))
  
  # first compute the GC content for each peak
  DefaultAssay(temp_obj2) <- "peaks"
  temp_obj2 <- RegionStats(temp_obj, 
                          genome = BSgenome.Hsapiens.UCSC.hg38)
  
  print("LinkPeaks() running")
  # link peaks to genes
  temp_obj2 <- LinkPeaks(
    object = temp_obj,
    peak.assay = "peaks",
    expression.assay = "SCT")
  print("LinkPeaks() finished")
  
  celltype_only_sep_list[[i]] <- temp_obj2
}


# loop over all objects to plot
for (i in 1:length(celltype_time_sep_list)){
  
  file_name <- paste0("../new_peak_set_plots/linkage/", ident_list[i], "linkage_plot.pdf")
  
  pdf(file = file_name)
  p1 <- CoveragePlot(
    object = obj,
    region = "BDNF",
    features = "BDNF",
    expression.assay = "SCT",
    group.by = "time.ident",
    extend.upstream = 10000,
    extend.downstream = 1000
  )
  print(p1)
  dev.off()
  
}


# gene activity analysis ####
gene.activities <- GeneActivity(
  obj_complete,
  assay = "peaks",
  extend.upstream = 500,
  extend.downstream = 0,
  biotypes = "protein_coding",
  max.width = 5e+05,
  verbose = TRUE)

obj_complete[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)

obj_complete <- NormalizeData(
  object = obj_complete,
  assay = 'GeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(obj_complete$nCount_RNA)
)
DefaultAssay(obj_complete) <- 'GeneActivity'

obj_complete_geneact <- RunUMAP(obj_complete, 
                                reduction = "lsi", 
                                dims = 1:30,
                                seed.use = 42)
obj_complete_geneact <- FindNeighbors(obj_complete_geneact, 
                                      reduction = "lsi", 
                                      dims = 1:30)
obj_complete_geneact <- FindClusters(obj_complete_geneact, 
                                     resolution = 0.5,
                                     graph.name = "ATAC_snn",
                                     random.seed = 42)

DimPlot(object = obj_complete_geneact, label = TRUE) + NoLegend()

FeaturePlot(object = integrated_renamed_1, 
            features = c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
                         "SEMA3E", "CUX2", "NEFM", 
                         "VIM", "SOX2", "NES"))


Idents(obj_complete_geneact) <- "seurat_clusters"

subtype_markers <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
                     "EBF1", "SERTAD4", # striatal
                     "FOXG1", # forebrain 
                     "TBR1", "TLE4", # pallial glutamatergic
                     "FEZF2", "ADCYAP1", "TCERG1L", "SEMA3E", # subcerebral
                     "POU3F2", "CUX1", "BCL11B",  # cortical
                     "LHX2", # general cortex
                     #"EOMES", # hindbrain
                     "NPY", "SST", "DLX2", "DLX5", # inhibitory
                     "SATB2", "CBLN2", "CUX2", "NEFM", # excitatory
                     "VIM", "SOX2", "NES" #NPC
)

StackedVlnPlot(obj = obj_complete_geneact, features = subtype_markers) +
  coord_flip()

save.image(file = "added_gene_activity_and_linkage.RData")

# pseudotime plot ####
heatmap_markers <- c()

df_4_heatmap <- data.fram()


heatmap.2(as.matrix(df_4_heatmap),
          Rowv = F,
          Colv = F,
          dendrogram = "none",
          scale = "none",
          col = colorRampPalette(colors = c("steelblue3", "steelblue1", "white", "red3"))(100),
          trace = "none",
          margins = c(10,10),
          cellnote = round(as.matrix(df_4_heatmap), digits = 2),
          notecol = "black")
dev.off()

