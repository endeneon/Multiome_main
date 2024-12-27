# Chuxuan Li 03/31/2023
# Look at percentage of reads mapped to mouse and rabbit genome, look at their 
#distribution and whether cells can be separated


# init ####
library(Seurat)
library(Signac)
load("~/NVME/scARC_Duan_018/Duan_project_029_ATAC/raw_029_ATACseq_objlist_w_TSS_nuc_sig_removed_nonhuman_chr.RData")

# read h5 raw data ####
h5list <- sort(list.files(path = "/nvmefs/scARC_Duan_029_Oc2_Mm10/",
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T,
                          include.dirs = F, 
                          full.names = T))
fraglist <- sort(list.files(path = "/nvmefs/scARC_Duan_029_Oc2_Mm10/",
                            pattern = "atac_fragments.tsv.gz$",
                            recursive = T,
                            include.dirs = F, 
                            full.names = T))
raw_obj_list <- vector(mode = "list", length = length(h5list))
libs <- str_extract(h5list, "[0-9]+-[0-6]")

for (i in 1:length(raw_obj_list)){
  h5file <- Read10X_h5(filename = h5list[i])
  cat("h5: ", str_extract(string = h5list[i],
                          pattern = "[0-9]+-[0-6]"))
  frag <- fraglist[i]
  chromAssay <- CreateChromatinAssay(
    counts = h5file$Peaks,
    sep = c(":", "-"),
    fragments = frag
  )
  obj <- CreateSeuratObject(counts = chromAssay, assay = "ATAC", 
                            project = str_extract(string = h5list[i],
                                                  pattern = "[0-9]+-[0-6]"))

  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list[i],
                                                  pattern = "-[0-6]"), "-"), "hr")
  print(unique(obj$time.ident))
  raw_obj_list[[i]] <- obj
}
names(raw_obj_list) <- libs

# count reads for single cells ####
for (i in 1:length(raw_obj_list)) {
  obj <- raw_obj_list[[i]]
  obj$num_reads_mapped_to_mouse <- 
    colSums(obj@assays$ATAC@counts[str_detect(obj@assays$ATAC@counts@Dimnames[[1]], 
                                              "mm10"), ])
  obj$num_reads_mapped_to_rabbit <- 
    colSums(obj@assays$ATAC@counts[str_detect(obj@assays$ATAC@counts@Dimnames[[1]], 
                                              "oc2"), ])
  obj$total_reads <- colSums(obj@assays$ATAC@counts)
  obj$pct_reads_mapped_to_mouse <- obj$num_reads_mapped_to_mouse / obj$total_reads
  obj$pct_reads_mapped_to_rabbit <- obj$num_reads_mapped_to_rabbit / obj$total_reads
  raw_obj_list[[i]] <- obj
}

# check distribution ####
hist(obj$pct_reads_mapped_to_mouse, breaks = 500)
hist(obj$pct_reads_mapped_to_rabbit, breaks = 500)
for (i in 1:length(raw_obj_list)) {
  obj <- raw_obj_list[[i]]
  png(filename = paste0("./rabbit_analysis/QC/pct_mapped_to_genome_histograms/pct_mapped_to_mouse_histogram_",
                        names(raw_obj_list)[i], ".png"), width = 900, height = 500)
  hist(obj$pct_reads_mapped_to_mouse, breaks = 500, 
       main = "Percentage of reads mapped to mouse mm10 genome - ATACseq", 
       xlab = "Percent reads mapped to mm10")
  dev.off()
  png(filename = paste0("./rabbit_analysis/QC/pct_mapped_to_genome_histograms/pct_mapped_to_rabbit_histogram_",
                        names(raw_obj_list)[i], ".png"), width = 900, height = 500)
  hist(obj$pct_reads_mapped_to_rabbit, breaks = 500,
       main = "Percentage of reads mapped to rabbit oc2 genome- ATACseq", 
       xlab = "Percent reads mapped to oc2")
  dev.off()
}

# total reads distribution ####
for (i in 1:length(raw_obj_list)) {
  obj <- raw_obj_list[[i]]
  png(filename = paste0("./rabbit_analysis/QC/total_reads_histograms/total_reads_histogram_",
                        names(raw_obj_list)[i], ".png"), width = 900, height = 500)
  hist(obj$total_reads, breaks = 500, 
       main = paste0("Total number of mapped reads - ", names(raw_obj_list)[i], " ATACseq"), 
       xlab = "Number of reads")
  dev.off()
}

# separate by pct reads mapped to rabbit ####

# # count reads ####
# reads_count_df <- data.frame(library = names(raw_obj_list),
#                              num_reads_mapped_to_mouse = rep_len(0, length(raw_obj_list)),
#                              num_reads_mapped_to_rabbit = rep_len(0, length(raw_obj_list)),
#                              total_reads = rep_len(0, length(raw_obj_list)),
#                              pct_reads_mapped_to_mouse = rep_len(0, length(raw_obj_list)),
#                              pct_reads_mapped_to_rabbit = rep_len(0, length(raw_obj_list)))
# unique(str_sub(obj@assays$ATAC@counts@Dimnames[[1]], end = 7L))
# for (i in 1:length(raw_obj_list)) {
#   obj <- raw_obj_list[[i]]
#   reads_count_df$num_reads_mapped_to_mouse[i] <- 
#     sum(rowSums(obj@assays$ATAC@counts[str_detect(obj@assays$ATAC@counts@Dimnames[[1]], "mm10"), ]))
#   reads_count_df$num_reads_mapped_to_rabbit[i] <- 
#     sum(rowSums(obj@assays$ATAC@counts[str_detect(obj@assays$ATAC@counts@Dimnames[[1]], "oc2"), ]))
#   reads_count_df$total_reads[i] <- sum(rowSums(obj@assays$ATAC@counts))
#   reads_count_df$pct_reads_mapped_to_mouse[i] <- 
#     reads_count_df$num_reads_mapped_to_mouse[i] / reads_count_df$total_reads[i]
#   reads_count_df$pct_reads_mapped_to_rabbit[i] <- 
#     reads_count_df$num_reads_mapped_to_rabbit[i] / reads_count_df$total_reads[i]
# }

