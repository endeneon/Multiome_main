# Chuxuan Li 02/14/2023
# Plot TSS enrichment, nucleosome signal, and other QC plots for 029 ATACseq

# init ####
library(Seurat)
library(Signac)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
library(readr)
library(readxl)

# read h5 raw data ####
h5list <- sort(list.files(path = "/data/FASTQ/Duan_Project_029/GRCh38_only",
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T,
                          include.dirs = F, 
                          full.names = T))
h5list <- h5list[1:15]
fraglist <- sort(list.files(path = "/data/FASTQ/Duan_Project_029/GRCh38_only",
                            pattern = "atac_fragments.tsv.gz$",
                            recursive = T,
                            include.dirs = F, 
                            full.names = T))
fraglist <- fraglist[1:15]
raw_obj_list <- vector(mode = "list", length = length(h5list))
libs <- str_extract(h5list, "[0-9]+-[0-6]")

load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")

for (i in 1:length(raw_obj_list)){
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
  
  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list[i],
                                                  pattern = "-[0-6]"), "-"), "hr")
  print(unique(obj$time.ident))
  raw_obj_list[[i]] <- obj
}
names(raw_obj_list) <- libs

# plot TSS and nucleosome signal by library ####
# add annotation
for (i in 1:length(raw_obj_list)){
  Signac::Annotation(raw_obj_list[[i]]) <- ens_use
}
# compute QC metrics
for (i in 1:length(raw_obj_list)){
  raw_obj_list[[i]] <- NucleosomeSignal(raw_obj_list[[i]])
  raw_obj_list[[i]] <- TSSEnrichment(object = raw_obj_list[[i]], fast = FALSE)
  
  pdf(file = paste0("./QC/TSS_plots/", str_extract(h5list[[i]], "[0-9]+-[0|1|6]"),
                    "_TSS_plot.pdf"))
  p <- TSSPlot(raw_obj_list[[i]]) + 
    NoLegend() +
    theme(text = element_text(size = 14))
  print(p)
  dev.off()
  png(paste0("./QC/TSS_plots/", str_extract(h5list[[i]], "[0-9]+-[0|1|6]"),
                    "_TSS_plot.png"))
  print(p)
  dev.off()
  
  pdf(file = paste0("./QC/nuc_signal_histograms/", str_extract(h5list[[i]], "[0-9]+-[0|1|6]"),
                    "_nucleosome_signal_hist.pdf"))
  p <- FragmentHistogram(raw_obj_list[[i]]) +
    theme(text = element_text(size = 14)) +
    ggtitle(str_extract(h5list[[i]], "[0-9]+-[0|1|6]"))
  print(p)
  dev.off()
  png(paste0("./QC/nuc_signal_histograms/", str_extract(h5list[[i]], "[0-9]+-[0|1|6]"),
                    "_nucleosome_signal_hist.png"))
  print(p)
  dev.off()
}
TSS_sum <- 0
num_cells <- 0
for (i in raw_obj_list){
  TSS_sum <- TSS_sum + sum(i$TSS.enrichment)
  num_cells <- num_cells + length(colnames(i))
}
TSS_sum / num_cells #average TSS: 5.074527

nuc_sig_sum <- 0
num_cells <- 0
for (i in raw_obj_list){
  print(sum(is.infinite(i$nucleosome_signal)))
  print(sum(is.nan(i$nucleosome_signal)))
  ns <- i$nucleosome_signal[!is.infinite(i$nucleosome_signal)]
  ns <- ns[!is.nan(ns)]
  nuc_sig_sum <- nuc_sig_sum + sum(ns)
  num_cells <- num_cells + length(colnames(i))
}
nuc_sig_sum / num_cells #average nucleosome signal: 0.5446985

save(raw_obj_list, file = "raw_029_ATACseq_objlist_w_TSS_nuc_sig_removed_nonhuman_chr.RData")

# load and extract TSS and FRiP data from 10x output ####
pathlist <- sort(list.files(path = "/data/FASTQ/Duan_Project_029/GRCh38_only",
                          pattern = "summary\\.csv",
                          recursive = T,
                          include.dirs = F, 
                          full.names = T))
pathlist <- pathlist[1:15]
sumcsv_list <- vector("list", length(pathlist))
for (i in 1:length(sumcsv_list)) {
  sumcsv_list[[i]] <- read_csv(pathlist[i])
}

summary_table <- data.frame(Library = str_extract(pathlist, "[0-9]+-[0-6]"),
                            `TSS score` = rep_len(0, length(pathlist)),
                            `Fraction of high-quality fragments overlapping TSS` = rep_len(0, length(pathlist)),
                            `Fraction of high-quality fragments overlapping peaks` = rep_len(0, length(pathlist)))
for (i in 1:length(sumcsv_list)) {
  sumcsv <- sumcsv_list[[i]]
  summary_table$Fraction.of.high.quality.fragments.overlapping.TSS[i] <- sumcsv$`ATAC Fraction of high-quality fragments overlapping TSS`
  summary_table$Fraction.of.high.quality.fragments.overlapping.peaks[i] <- sumcsv$`ATAC Fraction of high-quality fragments overlapping peaks`
  summary_table$TSS.score[i] <- sumcsv$`ATAC TSS enrichment score`
}
write.table(summary_table, file = "./QC/TSS_FRiP_FRiTSS_summary_stats.csv", quote = F,
            sep = ",")
