# Chuxuan Li 05/13/2022
# Generate barcodes in txt files for ASoC from 18-line data

# init ####
library(Seurat)
library(Signac)
library(stringr)
load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/ATAC_obj_added_MACS2_peaks.RData")


ATAC <- subset(ATAC, cell.type != "unknown")
times <- unique(ATAC$time.ident)
lines <- unique(ATAC$cell.line.ident)
types <- unique(ATAC$cell.type)

DefaultAssay(ATAC)

# loop ####
for (i in times) {
  for (j in lines) {
    for (k in types) {
      filename <- paste(i, j, k, "barcodes.txt", sep = "_")
      print(filename)
      tempobj <- subset(ATAC, time.ident == i & cell.line.ident == j & cell.type == k)
      cells <- colnames(tempobj)
      cells <- str_replace(cells, "-[0-9]+", "-1")
      write.table(cells, file = filename, quote = F, 
                  sep = "\t", row.names = F, col.names = F)
    }
  }
}
