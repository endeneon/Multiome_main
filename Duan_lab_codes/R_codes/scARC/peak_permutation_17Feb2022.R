# Siwei 17 Feb 2022

# permutation peaks by changing peakset identity

# init
library(Seurat)
library(Signac)
library(rtracklayer)
library(GenomicRanges)
library(readr)



# load database
# 5 lines ATAC-Seq data
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/ATACseq_5line_objonly_labeled.RData")
# 20 lines RNA-Seq data
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/demux_20line_integrated_labeled_obj.RData")

# set parallel computing
library(future)
plan("multicore", workers = 4)

# give each peak a unique identity, 
# use names from rownames(obj_complete@assays$peaks@data)
obj_complete@assays$peaks@diff.peaks <- F
DefaultAssay(obj_complete) <- "peaks"
# obj_complete@assays$peaks@var.features <-
#   rep_len(F, length.out = ncol(obj_complete@assays$peaks@data))
# obj_complete@assays$peaks@key <-
#   rep_len("unknown", length.out = ncol(obj_complete@assays$peaks@data))

obj_complete@assays$peaks@misc[["diff.peak.ident"]] <-
  rep_len(F, length.out = ncol(obj_complete@assays$peaks@data))
obj_complete@assays$peaks@misc[["peak.ident"]] <-
  rep_len("unknown", length.out = ncol(obj_complete@assays$peaks@data))
obj_complete@assays$peaks@misc[["peak.names"]] <-
  rownames(obj_complete@assays$peaks@data)

# assign metadata to peaks by read in bed files
# read in all barcode files by cell types
cell_type_barcodes <-
  list.files(path = "../R_peak_permutation/for_He_lab_Jan2022/peaks_called_by_cell_types/",
             pattern = "*.bed",
             full.names = T)
cell_types <- c("GABA", "NmCP_glut", "NPC_glut", "NpCm_glut")

i <- 1
# assign peak identity
for (i in 1:length(cell_type_barcodes)) {
  print(cell_types[i])
  temp_peaks <-
    data.frame(read_delim(cell_type_barcodes[i], 
                          delim = "\t", escape_double = FALSE, 
                          col_names = FALSE, trim_ws = TRUE),
               stringsAsFactors = F)
  temp_peaks$X2 <- as.character(temp_peaks$X2)
  temp_peaks$X3 <- as.character(temp_peaks$X3)
  temp_peak_names <- 
    apply(temp_peaks, 1, 
          function(x) paste(x[1], x[2], x[3], sep = '-'))
  
  
  obj_complete@assays$peaks@misc$peak.ident[obj_complete@assays$peaks@misc$peak.names %in% 
                                              temp_peak_names] <- cell_types[i]
}



peak_ranges_4_permutation <-
  sample(x = obj_complete@assays$peaks@ranges[obj_complete$],
         size = 10000,
         replace = F)
cells_4_permutation <-
  sample(x = obj_complete@as)


#####

