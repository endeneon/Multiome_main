# Chuxuan Li 03/07/2023
# output by-line-time-type barcodes for 025 data

# init ####
library(stringr)
library(Seurat)
load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/integrated_labeled.RData")

# output ####
lines <- unique(integrated_labeled$cell.line.ident)
times <- unique(integrated_labeled$time.ident)
types <- unique(integrated_labeled$cell.type)
types <- types[types != "unidentified"]

for (l in lines) {
  for (t in times) {
    for (y in types) {
      subobj <- subset(integrated_labeled, cell.line.ident == l & time.ident == t & cell.type == y)
      bc <- subobj@assays$RNA@counts@Dimnames[[2]]
      write.table(bc, file = paste0("./Analysis_part2_GRCh38/by_line_type_time_barcodes_025/025_",
                                    t, "_", l, "_", y, "_barcodes.txt"), 
                  quote = F, sep = "\t", row.names = F, col.names = F)
    }
  }
}
