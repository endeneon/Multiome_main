# Chuxuan li 03/29/2022
# load raw, unfiltered 18-line ATACseq data and determine which are rat astrocytes
#and map the remaining cells to human demuxed barcodes

# init ####
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)

library(ggplot2)
library(patchwork)

library(stringr)
library(future)

set.seed(2022)
plan("multisession", workers = 2)

# load data ####
pathlist <- list.files(path = "/data/FASTQ/Duan_Project_024/hybrid_output/",
                       pattern = "filtered_feature_bc_matrix.h5", full.names = T,
                       recursive = T)
tsvpathlist <- list.files(path = "/data/FASTQ/Duan_Project_024/hybrid_output/",
                          pattern = "atac_fragments.tsv.gz$", recursive = T,
                          full.names = T)
objlist <- vector(mode = "list", length = length(pathlist))

# create Seurat object from fragment files and h5 files
for (i in 1:length(pathlist)){
  counts <- Read10X_h5(pathlist[[i]])
  counts <- counts$Peaks
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = tsvpathlist[i]
  )
  
  objlist[[i]] <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC"
  )
}

save(objlist, file = "raw_obj_list.RData")

# Identify rat astrocytes (unsuccessful) ####
dimreduclust <- function(obj){
  obj <- RunUMAP(obj,
                 reduction = "lsi",
                 dims = 1:30,
                 seed.use = 11)
  obj <- FindNeighbors(obj,
                       reduction = "lsi",
                       dims = 1:30)
  obj <- FindClusters(obj,
                      resolution = 0.8,
                      random.seed = 11)
  return(obj)
  
}

processed_lst <- vector(mode = "list", length = length(objlist))
for (i in 1:length(objlist)){
  obj <- RunTFIDF(objlist[[i]])
  obj <- FindTopFeatures(obj, min.cutoff = 'q0')
  obj <- RunSVD(obj)
  processed_lst[[i]] <- dimreduclust(obj)
}

save(processed_lst, file = "unfiltered_normalized_list.RData")

DimPlot(processed_lst[[1]])

intervals <- objlist[[3]]@assays$ATAC@counts@Dimnames[[1]]
intervals[str_detect(intervals, "^10-9100")]

setwd("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/QC_for_select_rat_plots/")
for (i in 3){
  jpeg(paste0(str_extract(pathlist[i], "[0-9]+-[0-6]"), "_dimplot_by_cluster.jpeg"))
  p <- DimPlot(processed_lst[[i]],
               label = T,
               group.by = "seurat_clusters") +
    NoLegend() +
    ggtitle(str_extract(pathlist[i], "[0-9]+-[0-6]")) +
    theme(text = element_text(size = 12))
  print(p)
  dev.off()
  
  jpeg(paste0(str_extract(pathlist[i], "[0-9]+-[0-6]"), "_Gfap_1.jpeg"))
  p <- FeaturePlot(processed_lst[[i]], features = "10-91000437-91001327") +
    ggtitle(str_extract(pathlist[i], "[0-9]+-[0-6]"))
  print(p)
  dev.off()
  jpeg(paste0(str_extract(pathlist[i], "[0-9]+-[0-6]"), "_Gfap_2.jpeg"))
  p <- FeaturePlot(processed_lst[[i]], features = "10-91002878-91003745") +
    ggtitle(str_extract(pathlist[i], "[0-9]+-[0-6]"))
  print(p)
  dev.off()
  
  jpeg(paste0(str_extract(pathlist[i], "[0-9]+-[0-6]"), "_Gfap_vln_1.jpeg"))
  p <- VlnPlot(processed_lst[[i]], "10-91002878-91003750") +
    ggtitle(paste0(str_extract(pathlist[i], "[0-9]+-[0-6]"), " Gfap expression by cluster"))
  print(p)
  dev.off()
  jpeg(paste0(str_extract(pathlist[i], "[0-9]+-[0-6]"), "_Gfap_vln_2.jpeg"))
  p <- VlnPlot(processed_lst[[i]], "10-91000429-91001327") +
    ggtitle(paste0(str_extract(pathlist[i], "[0-9]+-[0-6]"), " Gfap expression by cluster"))
  print(p)
  dev.off()
}


# separate rat from data first ####
# sum up the number of reads for each human/rat interval for mapped human cells,
#then look at the distribution and determine the average exp level for human/rat
#intervals for human cells, use this average to separate the rest of the cells

for (i in 1:length(objlist)){
  # find human vs other intervals
  all.intv <- objlist[[i]]@assays$ATAC@counts@Dimnames[[1]]
  human.intv <- all.intv[str_detect(all.intv, "^chr")]
  other.intv <- setdiff(all.intv, human.intv)
  # sum up number of reads in human/other intervals for each cell
  objlist[[i]]$total_human_interval_reads <- 
    colSums(objlist[[i]]@assays$ATAC@counts[which(rownames(objlist[[i]]@assays$ATAC@counts) %in% human.intv), ])
  print(mean(objlist[[i]]$total_human_interval_reads))
  objlist[[i]]$total_rat_interval_reads <- 
    colSums(objlist[[i]]@assays$ATAC@counts[which(rownames(objlist[[i]]@assays$ATAC@counts) %in% other.intv), ])
  print(mean(objlist[[i]]$total_rat_interval_reads))
}

# take a look at the distribution of total counts
human_reads_combined <- c()
rat_reads_combined <- c()
for (i in 1:length(objlist)){
  human_reads_combined <- c(human_reads_combined, objlist[[i]]$total_human_interval_reads)
  rat_reads_combined <- c(rat_reads_combined, objlist[[i]]$total_rat_interval_reads)
}

hist(human_reads_combined, breaks = 10000, xlim = c(0, 100000))
hist(rat_reads_combined, breaks = 10000, xlim = c(0, 10000))
median(human_reads_combined) #12049
median(rat_reads_combined) #13

# determine how to filter the rat cells out
# sum(objlist[[1]]$total_rat_interval_reads > median(objlist[[1]]$total_rat_interval_reads) & 
#       objlist[[1]]$total_human_interval_reads < quantile(objlist[[1]]$total_human_interval_reads,
#                                                          0.9))

# remove rat, save number of rat cells per lib in df
human_rat_cell_count_df <- data.frame(human_mapped = rep_len(0, length(objlist)),
                                      human_unmapped = rep_len(0, length(objlist)),
                                      rat = rep_len(0, length(objlist)))
norat_lst <- vector("list", length(objlist))
for (i in 1:length(objlist)){
  objlist[[i]]$rat.ident <- "human"
  objlist[[i]]$rat.ident[objlist[[i]]$total_rat_interval_reads > 
                           0.25 * mean(objlist[[i]]$total_rat_interval_reads)] <- "rat"
  print(unique(objlist[[i]]$rat.ident))
  human_rat_cell_count_df$rat[i] <- sum(objlist[[i]]$rat.ident == "rat")
  norat_lst[[i]] <- subset(objlist[[i]], rat.ident == "human")
}

save(norat_lst, file = "no_rat_astro_list.RData")
# now map the cells to demuxed barcodes ####

# read barcode .best files
setwd("/data/FASTQ/Duan_Project_024/hybrid_output")
pathlist <- list.files(path = ".", full.names = T, pattern = ".best.tsv", recursive = T)
pathlist <- pathlist[str_detect(pathlist, "demux_0.01", T)]
barcode_list <- vector(mode = "list", length = length(pathlist))
for (i in 1:length(pathlist)){
  barcode_list[[i]] <- read.delim(pathlist[i], header = F, row.names = NULL)
  colnames(barcode_list[[i]]) <- c("barcode", "line")
}

# save mapped human cells in list
cleanobj_lst <- vector(mode = "list", length(norat_lst))

# sep lines in each group's barcode list, then assign line by barcode
for (i in 1:length(barcode_list)){
  obj <- norat_lst[[i]]
  # total 15 barcode lists for 15 libraries, each lib has 3-4 lines
  lines <- unique(barcode_list[[i]]$line)
  print(lines)
  # initialize
  obj$cell.line.ident <- "unmatched"
  for (j in 1:length(lines)){
    # each lib has a list of lines and corresponding barcodes
    line_spec_barcodes <- barcode_list[[i]]$barcode[barcode_list[[i]]$line == lines[j]]
    obj$cell.line.ident[obj@assays$ATAC@counts@Dimnames[[2]] %in% 
                                         line_spec_barcodes] <- lines[j]
  }
  print(unique(obj$cell.line.ident))
  # save numbers of unmapped and mapped to df
  human_rat_cell_count_df$human_mapped[i] <- sum(obj$cell.line.ident != "unmatched")
  human_rat_cell_count_df$human_unmapped[i] <- sum(obj$cell.line.ident == "unmatched")
  # keep only mapped cells
  obj <- subset(obj, cell.line.ident != "unmatched")
  # pick only human intervals
  counts <- GetAssayData(obj, assay = "ATAC")
  use.intv <- obj@assays$ATAC@counts@Dimnames[[1]][str_detect(obj@assays$ATAC@counts@Dimnames[[1]], "^chr")]
  counts <- counts[(which(rownames(counts) %in% use.intv)), ]
  cleanobj_lst[[i]] <- subset(obj, features = rownames(counts))
}

setwd("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/")
write.table(human_rat_cell_count_df, file = "human_rat_cell_count_table.csv",
            quote = F, row.names = str_extract(tsvpathlist, "[0-9]+-[0|1|6]"), 
            col.names = T, sep = ",")
save(cleanobj_lst, file = "no_rat_no_unmapped_no_rat_gene_clean_list.RData")
