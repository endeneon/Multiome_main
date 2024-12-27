# Chuxuan Li 02/21/2023
# pool all 018-029 sequencing quality info for redo purposes

# init ####
library(readr)
library(Signac)
library(Seurat)

# 029 ####
load("/nvmefs/scARC_Duan_018/Duan_project_029_ATAC/raw_029_ATACseq_objlist_w_TSS_nuc_sig_removed_nonhuman_chr.RData")
h5list <- sort(list.files(path = "/data/FASTQ/Duan_Project_029/GRCh38_only/",
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T,
                          include.dirs = F, 
                          full.names = T))
h5list <- h5list[1:15]
libnames <- str_extract(string = h5list,
                        pattern = "[0-9]+-[0-6]")

pathlist_summary <- sort(list.files(path = "/data/FASTQ/Duan_Project_029/GRCh38_only",
                            pattern = "summary\\.csv",
                            recursive = T,
                            include.dirs = F, 
                            full.names = T))
pathlist_summary <- pathlist_summary[1:15]
sumcsv_list <- vector("list", length(pathlist_summary))
for (i in 1:length(sumcsv_list)) {
  sumcsv_list[[i]] <- read_csv(pathlist_summary[i])
}
df_029 <- data.frame(library = libnames,
                     cell.counts = rep_len(0, length(libnames)),
                     n.UMI = rep_len(0, length(libnames)),
                     UMI.per.cell = rep_len(0, length(libnames)),
                     n.raw.reads = rep_len(0, length(libnames)),
                     raw.reads.per.cell = rep_len(0, length(libnames)))

for (i in 1:length(libnames)){
  print(libnames[i])
  df_029$cell.counts[i] <- length(raw_obj_list[[i]]@assays$ATAC@counts@Dimnames[[2]])
  df_029$n.UMI[i] <- sum(raw_obj_list[[i]]@assays$ATAC@counts@x)
  df_029$UMI.per.cell[i] <- 
    sum(raw_obj_list[[i]]@assays$ATAC@counts@x)/df_029$cell.counts[i]
  
  sumcsv <- sumcsv_list[[i]]
  df_029$raw.reads.per.cell[i] <- sumcsv$`ATAC Mean raw read pairs per cell`
  df_029$n.raw.reads[i] <- sumcsv$`Estimated number of cells` * sumcsv$`ATAC Mean raw read pairs per cell`
}

#save(df_029, file = "../../df_029.RData")

# 025 ####
load("/nvmefs/scARC_Duan_018/Duan_project_025_ATAC/ATAC_uncombined_obj_list_for_QC.RData")

h5list <- sort(list.files(path = "/nvmefs/scARC_Duan_025_GRCh38/",
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T,
                          include.dirs = F, 
                          full.names = T))
h5list <- h5list[2:16]
libnames <- str_extract(string = h5list,
                        pattern = "[0-9]+-[0-6]")

pathlist_summary <- sort(list.files(path = "/nvmefs/scARC_Duan_025_GRCh38/",
                                    pattern = "summary\\.csv",
                                    recursive = T,
                                    include.dirs = F, 
                                    full.names = T))
pathlist_summary <- pathlist_summary[2:16]
sumcsv_list <- vector("list", length(pathlist_summary))
for (i in 1:length(sumcsv_list)) {
  sumcsv_list[[i]] <- read_csv(pathlist_summary[i])
}
df_025 <- data.frame(library = libnames,
                     cell.counts = rep_len(0, length(libnames)),
                     n.UMI = rep_len(0, length(libnames)),
                     UMI.per.cell = rep_len(0, length(libnames)),
                     n.raw.reads = rep_len(0, length(libnames)),
                     raw.reads.per.cell = rep_len(0, length(libnames)))
for (i in 1:length(libnames)){
  print(libnames[i])
  df_025$cell.counts[i] <- length(raw_obj_list[[i]]@assays$ATAC@counts@Dimnames[[2]])
  df_025$n.UMI[i] <- sum(raw_obj_list[[i]]@assays$ATAC@counts@x)
  df_025$UMI.per.cell[i] <-
    sum(raw_obj_list[[i]]@assays$ATAC@counts@x)/df_025$cell.counts[i]

  sumcsv <- sumcsv_list[[i]]
  df_025$raw.reads.per.cell[i] <- sumcsv$`ATAC Mean raw read pairs per cell`
  df_025$n.raw.reads[i] <- sumcsv$`Estimated number of cells` * sumcsv$`ATAC Mean raw read pairs per cell`
}

#save(df_025, file = "../../df_025.RData")

# 024 ####
load("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/added_TSS_no_rat_no_unmapped_cleanobj_list.RData")
h5list <- sort(list.files(path = "/nvmefs/scARC_Duan_024_raw_fastq/GRCh38_output/",
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T, include.dirs = F, full.names = T))
libnames <- str_extract(string = h5list,
                        pattern = "[0-9]+-[0-6]")

pathlist_summary <- sort(list.files(path = "/nvmefs/scARC_Duan_024_raw_fastq/GRCh38_output/",
                                    pattern = "summary\\.csv",
                                    recursive = T,
                                    include.dirs = F, 
                                    full.names = T))
sumcsv_list <- vector("list", length(pathlist_summary))
for (i in 1:length(sumcsv_list)) {
  sumcsv_list[[i]] <- read_csv(pathlist_summary[i])
}

df_024 <- data.frame(library = libnames,
                     cell.counts = rep_len(0, length(libnames)),
                     n.UMI = rep_len(0, length(libnames)),
                     UMI.per.cell = rep_len(0, length(libnames)),
                     n.raw.reads = rep_len(0, length(libnames)),
                     raw.reads.per.cell = rep_len(0, length(libnames)))

for (i in 1:length(libnames)){
  print(libnames[i])
  df_024$cell.counts[i] <- length(raw_obj_list[[i]]@assays$ATAC@counts@Dimnames[[2]])
  df_024$n.UMI[i] <- sum(raw_obj_list[[i]]@assays$ATAC@counts@x)
  df_024$UMI.per.cell[i] <-
    sum(raw_obj_list[[i]]@assays$ATAC@counts@x)/df_024$cell.counts[i]
  
  sumcsv <- sumcsv_list[[i]]
  df_024$raw.reads.per.cell[i] <- sumcsv$`ATAC Mean raw read pairs per cell`
  df_024$n.raw.reads[i] <- sumcsv$`Estimated number of cells` * sumcsv$`ATAC Mean raw read pairs per cell`
}

#save(df_024, file = "../../df_024.RData")

# 018 ####
h5list <- sort(list.files(path = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/",
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T,
                          include.dirs = F, 
                          full.names = T))
fraglist <- sort(list.files(path = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/",
                            pattern = "atac_fragments.tsv.gz$",
                            recursive = T,
                            include.dirs = F, 
                            full.names = T))
raw_obj_list <- vector(mode = "list", length = length(h5list))
libnames <- str_replace(str_extract(string = h5list,
                                    pattern = "[0-9]+_[0-6]"), "_", "-")
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")

for (i in 1:length(raw_obj_list)){
  h5file <- Read10X_h5(filename = h5list[i])
  cat("h5: ", str_extract(string = h5list[i],
                          pattern = "[0-9]+_[0-6]"))
  frag <- fraglist[i]
  chromAssay <- CreateChromatinAssay(
    counts = h5file$Peaks,
    sep = c(":", "-"),
    fragments = frag,
    annotation = ens_use
  )
  obj <- CreateSeuratObject(counts = chromAssay, assay = "ATAC", 
                            project = str_extract(string = h5list[i],
                                                  pattern = "[0-9]+_[0-6]"))
  #use.intv <- obj@assays$ATAC@counts@Dimnames[[1]][str_detect(obj@assays$ATAC@counts@Dimnames[[1]], "^chr")]
  #counts <- obj@assays$ATAC@counts[(which(rownames(obj@assays$ATAC@counts) %in% use.intv)), ]
  #obj <- subset(obj, features = rownames(counts))
  
  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list[i],
                                                  pattern = "[0-9]+_[0-6]"), "[0-9]+_"), "hr")
  print(unique(obj$time.ident))
  raw_obj_list[[i]] <- obj
}

pathlist_summary <- sort(list.files(path = "/nvmefs/scARC_Duan_018/GRCh38_mapped_only/",
                                    pattern = "summary\\.csv",
                                    recursive = T,
                                    include.dirs = F, 
                                    full.names = T))
sumcsv_list <- vector("list", length(pathlist_summary))
for (i in 1:length(sumcsv_list)) {
  sumcsv_list[[i]] <- read_csv(pathlist_summary[i])
}

df_018 <- data.frame(library = libnames,
                     cell.counts = rep_len(0, length(libnames)),
                     n.UMI = rep_len(0, length(libnames)),
                     UMI.per.cell = rep_len(0, length(libnames)),
                     n.raw.reads = rep_len(0, length(libnames)),
                     raw.reads.per.cell = rep_len(0, length(libnames)))

for (i in 1:length(libnames)){
  print(libnames[i])
  df_018$cell.counts[i] <- length(raw_obj_list[[i]]@assays$ATAC@counts@Dimnames[[2]])
  df_018$n.UMI[i] <- sum(raw_obj_list[[i]]@assays$ATAC@counts@x)
  df_018$UMI.per.cell[i] <-
    sum(raw_obj_list[[i]]@assays$ATAC@counts@x)/df_018$cell.counts[i]
  
  sumcsv <- sumcsv_list[[i]]
  df_018$raw.reads.per.cell[i] <- sumcsv$`ATAC Mean raw read pairs per cell`
  df_018$n.raw.reads[i] <- sumcsv$`Estimated number of cells` * sumcsv$`ATAC Mean raw read pairs per cell`
}

#save(df_018, file = "../../df_018.RData")

# Combine all ####
df_018$sequencing.batch <- "018"
df_024$sequencing.batch <- "024"
df_025$sequencing.batch <- "025"
df_029$sequencing.batch <- "029"
df_combined <- dplyr::bind_rows(df_024, df_025, df_029) 

write.csv(df_combined, file = "018-029_ATACseq_cell_raw_read_UMI_counts_summary.csv", 
          quote = F, row.names = F)
