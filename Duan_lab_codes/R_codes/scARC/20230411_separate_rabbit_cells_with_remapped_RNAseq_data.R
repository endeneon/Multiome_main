# Chuxuan Li 04/11/2023
# Look at percentage of reads mapped to mouse and rabbit genome from RNAseq
#remapped data from Siwei, separate rabbit cells

# init ####
library(Seurat)
library(ggplot2)
library(pals)
library(stringr)
library(readxl)
library(readr)
library(future)

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

# read h5, remove rat cells and genes ####
h5list <- sort(list.files(path = "/nvmefs/scARC_Duan_029_Oc2_Mm10/",
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T,
                          include.dirs = F, 
                          full.names = T))
objlist <- vector(mode = "list", length = length(h5list))
libs <- str_extract(h5list, "[0-9]+-[0-6]")
for (i in 1:length(objlist)){
  h5file <- Read10X_h5(filename = h5list[i])
  cat("h5: ", str_extract(string = h5list[i],
                          pattern = "[0-9]+-[0-6]"))
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(string = h5list[i],
                                                  pattern = "[0-9]+-[0-6]"))
  print(ncol(obj))
  cell_counts$total[i] <- ncol(obj)
  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list[i],
                                                  pattern = "-[0-6]"), "-"), "hr")
  print(unique(obj$time.ident))
  # assign rat identity
  objlist[[i]] <- obj
}

# map to demuxed barcodes, then remove human cells ####
pathlist <- sort(list.files(path = "./Duan_029_human_only_demux_barcodes",
                            full.names = T, pattern = ".txt", recursive = T))
pathlist <- pathlist[c(4:6, 10:12)]
barcode_list <- vector(mode = "list", length = length(pathlist))
for (i in 1:length(pathlist)){
  barcode_list[[i]] <- read_table(pathlist[i], col_names = FALSE)
  colnames(barcode_list[[i]]) <- c("barcode", "line")
}

# sep lines in each group's barcode list, then assign line by barcode
for (i in 1:length(barcode_list)){
  lines <- unique(barcode_list[[i]]$line)
  print(lines)
  obj <- objlist[[i]]
  obj$cell.line.ident <- "unmatched"
  for (j in 1:length(lines)){
    line_spec_barcodes <- barcode_list[[i]]$barcode[barcode_list[[i]]$line == lines[j]]
    obj$cell.line.ident[obj@assays$RNA@counts@Dimnames[[2]] %in% 
                          line_spec_barcodes] <- lines[j]
  }
  obj$human.ident <- "nonhuman"
  obj$human.ident[obj$cell.line.ident != "unmatched"] <- "human"
  objlist[[i]] <- obj
}

# remove human cells 
nohuman_list <- vector(mode = "list", length = length(objlist))
for (i in 1:length(objlist)){
  print(unique(objlist[[i]]$human.ident))
  nohuman_list[[i]] <- subset(objlist[[i]], human.ident == "nonhuman")
  print(unique(nohuman_list[[i]]$human.ident))
}

# count reads for single cells ####
for (i in 1:length(nohuman_list)) {
  obj <- nohuman_list[[i]]
  all_genes <- obj@assays$RNA@counts@Dimnames[[1]]
  rabbit_genes_ind <- str_detect(all_genes, "^[A-Z][A-Z0-9]+")
  mouse_genes_ind <- !rabbit_genes_ind[str_detect(all_genes, "5S-rRNA*", T)]
  obj$num_reads_mapped_to_rabbit <- 
    colSums(obj@assays$RNA@counts[rabbit_genes_ind, ])
  obj$num_reads_mapped_to_mouse <- 
    colSums(obj@assays$RNA@counts[mouse_genes_ind, ])
  obj$total_reads <- colSums(obj@assays$RNA@counts)
  obj$pct_reads_mapped_to_rabbit <- obj$num_reads_mapped_to_rabbit / obj$total_reads
  obj$pct_reads_mapped_to_mouse <- obj$num_reads_mapped_to_mouse / obj$total_reads
  nohuman_list[[i]] <- obj
}
names(nohuman_list) <- libs
save(nohuman_list, file = "RNA_seq_nohuman_list_only_libs_w_rabbit_mapped_to_Feb2023_new_genome.RData")

# check distribution ####
hist(obj$pct_reads_mapped_to_mouse, breaks = 500)
hist(obj$pct_reads_mapped_to_rabbit, breaks = 500)
for (i in 1:length(nohuman_list)) {
  obj <- nohuman_list[[i]]
  png(filename = paste0("./oc2_mm10_analysis/histograms_new_mapped_data/pct_reads_mapped_to_mouse_histogram_",
                        names(nohuman_list)[i], ".png"), width = 900, height = 500)
  hist(obj$pct_reads_mapped_to_mouse, breaks = 500, 
       main = "Percentage of reads mapped to mouse mm10 genome - ATACseq", 
       xlab = "Percent reads mapped to mm10")
  dev.off()
  png(filename = paste0("./oc2_mm10_analysis/histograms_new_mapped_data/pct_reads_mapped_to_rabbit_histogram_",
                        names(nohuman_list)[i], ".png"), width = 900, height = 500)
  hist(obj$pct_reads_mapped_to_rabbit, breaks = 500,
       main = "Percentage of reads mapped to rabbit oc2 genome- ATACseq", 
       xlab = "Percent reads mapped to oc2")
  dev.off()
}


