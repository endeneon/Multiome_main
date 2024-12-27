# Chuxuan Li 11/07/2022
# demultiplex, remove nonhuman, QC

# init ####
library(Seurat)
library(Signac)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
library(readr)
library(readxl)

# read h5, remove rat cells and genes ####
h5list <- sort(list.files(path = "/nvmefs/scARC_Duan_025_GRCh38",
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T,
                          include.dirs = F, 
                          full.names = T))
h5list <- h5list[2:16]
fraglist <- sort(list.files(path = "/nvmefs/scARC_Duan_025_GRCh38",
                          pattern = "atac_fragments.tsv.gz$",
                          recursive = T,
                          include.dirs = F, 
                          full.names = T))
fraglist <- fraglist[2:16]
ATACobjlist <- vector(mode = "list", length = length(h5list))
libs <- str_extract(h5list, "[0-9]+-[0-6]")
cell_counts <- data.frame(library = libs,
                          total = rep(0, by = length(libs)),
                          rat = rep(0, by = length(libs)),
                          nonrat = rep(0, by = length(libs)),
                          human = rep(0, by = length(libs)),
                          nonhuman = rep(0, by = length(libs)))
load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/GRCh38_mapped_demux_matched_human_only_list.RData")
load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/GRCh38_mapped_removed_rat_cells_assigned_demux_bc_list.RData")
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")

for (i in 1:length(ATACobjlist)){
  h5file <- Read10X_h5(filename = h5list[i])
  cat("h5: ", str_extract(string = h5list[i],
                          pattern = "[0-9]+-[0-6]"))
  frag <- fraglist[i]
  chromAssay <- CreateChromatinAssay(
    counts = h5file$Peaks,
    sep = c(":", "-"),
    fragments = frag,
    annotation = ens_use
  )
  obj <- CreateSeuratObject(counts = chromAssay, assay = "ATAC", 
                             project = str_extract(string = h5list[i],
                                                   pattern = "[0-9]+-[0-6]"))
  use.intv <- obj@assays$ATAC@counts@Dimnames[[1]][str_detect(obj@assays$ATAC@counts@Dimnames[[1]], "^chr")]
  counts <- obj@assays$ATAC@counts[(which(rownames(obj@assays$ATAC@counts) %in% use.intv)), ]
  obj <- subset(obj, features = rownames(counts))
  
  cell_counts$total[i] <- ncol(obj)
  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list[i],
                                                  pattern = "-[0-6]"), "-"), "hr")
  print(unique(obj$time.ident))
  # assign rat identity
  humanbc <- colnames(objlist[[i]])
  obj$rat.ident <- "rat"
  obj$rat.ident[colnames(obj) %in% humanbc] <- "human"
  cell_counts$rat[i] <- sum(obj$rat.ident == "rat")
  cell_counts$nonrat[i] <- sum(obj$rat.ident == "human")
  obj$cell.line.ident <- "unmatched"
  lines <- unique(objlist[[i]]$cell.line.ident)
  for (l in lines) {
    obj$cell.line.ident[colnames(obj) %in% colnames(objlist[[i]])[objlist[[i]]$cell.line.ident == l]] <- l
  }
  cell_counts$human[i] <- sum(obj$cell.line.ident != "unmatched")
  cell_counts$nonhuman[i] <- sum(obj$cell.line.ident == "unmatched")
  ATACobjlist[[i]] <- obj
}

write.table(cell_counts, file = "ATAC_cell_counts_matched_to_RNAseq.csv", 
            sep = ",", quote = F, col.names = T, row.names = F)


# 1. summarize number of cells and counts per library ####
plotCellCountsRectangles <- function(df) {
  p <- ggplot(df, aes(x = library, y = rep_len(2, length(libnames)), fill = cell.counts)) +
    geom_col(width = 1, position = position_dodge(width = 0.5)) +
    geom_text(aes(label = cell.counts), nudge_y = -1, hjust = "middle") + 
    scale_fill_continuous(low = "#a2c1f2", high = "#5c96f2") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          axis.title = element_blank(), axis.text.x = element_blank(),
          axis.ticks = element_blank()) + 
    ggtitle("# Nuclei") +
    coord_flip() +
    NoLegend()
  print(p)
}
plotUMIBars <- function(df) {
  p <- ggplot(df, aes(x = library, y = reads.per.cell)) +
    ylab("# UMIs per nuclei") +
    geom_col(position = "dodge", color = "gray27", fill = "gray90") +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 8)) + 
    coord_flip()
  print(p)
}

plotTotalUMIBars <- function(df) {
  p <- ggplot(df, aes(x = library, y = n.reads)) +
    ylab("# total UMIs") +
    geom_col(position = "dodge", color = "gray27", fill = "gray90") +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 8)) + 
    coord_flip()
  print(p)
}

libnames <- str_extract(string = h5list,
                        pattern = "[0-9]+-[0-6]")
df_raw <- data.frame(library = libnames,
                     cell.counts = rep_len(0, length(libnames)),
                     n.reads = rep_len(0, length(libnames)),
                     reads.per.cell = rep_len(0, length(libnames)))
df_human <- data.frame(library = libnames,
                       cell.counts = rep_len(0, length(libnames)),
                       n.reads = rep_len(0, length(libnames)),
                       reads.per.cell = rep_len(0, length(libnames)))
for (i in 1:length(libnames)){
  print(libnames[i])
  rawobj <- ATACobjlist[[i]]
  humanonlyobj <- subset(rawobj, cell.line.ident != "unmatched")
  df_raw$cell.counts[i] <- length(rawobj@assays$Peaks@counts@Dimnames[[2]])
  df_raw$n.reads[i] <- sum(rawobj@assays$Peaks@counts@x)
  df_raw$reads.per.cell[i] <- sum(rawobj@assays$Peaks@counts@x)/df_raw$cell.counts[i]
  df_human$cell.counts[i] <- length(humanonlyobj@assays$Peaks@counts@Dimnames[[2]])
  df_human$n.reads[i] <- sum(humanonlyobj@assays$Peaks@counts@x)
  df_human$reads.per.cell[i] <- sum(humanonlyobj@assays$Peaks@counts@x)/df_raw$cell.counts[i]
}
write.table(df_raw, 
            file = "./QC/raw_cell_counts_UMI_counts_per_lib.csv",
            quote = F, sep = ",", row.names = F,
            col.names = T)
write.table(df_human, 
            file = "./QC/human_cell_counts_UMI_counts_per_lib.csv",
            quote = F, sep = ",", row.names = F,
            col.names = T)
df_raw$library <- factor(df_raw$library, levels = rev(libnames))
df_human$library <- factor(df_human$library, levels = rev(libnames))
plotCellCountsRectangles(df_human)
plotUMIBars(df_human)
plotTotalUMIBars(df_human)
save(ATACobjlist, file = "ATAC_uncombined_obj_list_for_QC.RData")

# plot TSS and nucleosome signal for each library ####
# add annotation
for (i in 1:length(ATACobjlist)){
  Signac::Annotation(ATACobjlist[[i]]) <- ens_use
}
# compute QC metrics
for (i in 1:length(ATACobjlist)){
  # ATACobjlist[[i]] <- NucleosomeSignal(ATACobjlist[[i]])
  # ATACobjlist[[i]] <- TSSEnrichment(object = ATACobjlist[[i]], fast = FALSE)
  # 
  pdf(file = paste0("./QC/TSS_plots/", str_extract(h5list[[i]], "[0-9]+-[0|1|6]"),
                         "_TSS_plot.pdf"))
  p <- TSSPlot(ATACobjlist[[i]]) + 
    NoLegend() +
    theme(text = element_text(size = 14))
  print(p)
  dev.off()
  
  pdf(file = paste0("./QC/nuc_signal_histograms/", str_extract(h5list[[i]], "[0-9]+-[0|1|6]"),
                         "_nucleosome_signal_hist.pdf"))
  p <- FragmentHistogram(ATACobjlist[[i]]) +
    theme(text = element_text(size = 14))
  print(p)
  dev.off()
}
TSS_sum <- 0
num_cells <- 0
for (i in cleanobj_lst){
  TSS_sum <- TSS_sum + sum(i$TSS.enrichment)
  num_cells <- num_cells + length(colnames(i))
}
TSS_sum / num_cells

nuc_sig_sum <- 0
num_cells <- 0
for (i in cleanobj_lst){
  print(sum(is.infinite(i$nucleosome_signal)))
  print(sum(is.nan(i$nucleosome_signal)))
  ns <- i$nucleosome_signal[!is.infinite(i$nucleosome_signal)]
  ns <- ns[!is.nan(ns)]
  nuc_sig_sum <- nuc_sig_sum + sum(ns)
  num_cells <- num_cells + length(colnames(i))
}
nuc_sig_sum / num_cells
