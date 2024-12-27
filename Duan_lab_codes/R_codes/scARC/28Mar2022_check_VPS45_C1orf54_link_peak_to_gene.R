# Chuxuan Li 03/28/2022
# Check the gene-peak links for peaks near VPS45 and C1orf54

library(GenomicRanges)
library(Seurat)
library(Signac)

library(readr)

# read data ####
setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/link_peak_to_genes_results")
pathlist <- list.files("./", pattern = "*.csv")
celltype_time_sep_list <- vector(mode = "list", length = length(pathlist))
for (i in 1:length(pathlist)){
  celltype_time_sep_list[[i]] <- as.data.frame(read_csv(pathlist[i]),row.names = F)
}
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/new_peak_set_after_motif.RData")
obj_complete$celltype.time.ident[is.na(obj_complete$celltype.time.ident)] <- "NA"
obj_complete <- subset(obj_complete, celltype.time.ident != "NA")
idents_list <- unique(obj_complete$celltype.time.ident)
idents_list <- sort(idents_list)
idents_list <- idents_list[4:15]

# find the genes of interest ####
selected_rows <- vector(mode = "list", length = length(celltype_time_sep_list))
for(i in 1:length(celltype_time_sep_list)){
  selected_rows[[i]] <- celltype_time_sep_list[[i]][celltype_time_sep_list[[i]]$gene %in% 
                                      c("VPS45", "C1orf54"), ]
}

idents_list <- unique(obj_complete$celltype.time.ident)
idents_list <- sort(idents_list)
idents_list <- idents_list[4:15]

for (i in 1:length(idents_list)){
  print(idents_list[i])
  tempobj <- subset(obj_complete, celltype.time.ident == idents_list[i])
  links <- makeGRangesFromDataFrame(celltype_time_sep_list[[i]], keep.extra.columns = T)
  Links(tempobj) <- links
  pdf(file = paste0(idents_list[i], "_coverage_plot.pdf"), width = 8, height = 10)
  p <- CoveragePlot(
    object = tempobj,
    region = "VPS45",
    expression.assay = "peaks",
    extend.upstream = 10000,
    extend.downstream = 150000
  )
  print(p)
  dev.off()
  
}
