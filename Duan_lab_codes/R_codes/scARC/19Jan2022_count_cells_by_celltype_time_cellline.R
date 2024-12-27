# Chuxuan Li 01/19/2022
# count number of cells in each cell type at each time point from each cell type

library(Seurat)
library(Signac)
library(ggplot2)
library(readr)
library(stringr)
library(RColorBrewer)

# check cell line composition of each cluster
# first read in all the cell line barcode files
g_2_0_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_0_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
g_2_0_CD_27_barcodes <- unlist(g_2_0_CD_27_barcodes)

g_2_0_CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_0_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
g_2_0_CD_54_barcodes <- unlist(g_2_0_CD_54_barcodes)


g_2_1_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_1_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
g_2_1_CD_27_barcodes <- unlist(g_2_1_CD_27_barcodes)
g_2_1_CD_27_barcodes <- gsub(pattern = "-1",
                             replacement = "-2",
                             x = g_2_1_CD_27_barcodes)

g_2_1_CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_1_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
g_2_1_CD_54_barcodes <- unlist(g_2_1_CD_54_barcodes)
g_2_1_CD_54_barcodes <- gsub(pattern = "-1",
                             replacement = "-2",
                             x = g_2_1_CD_54_barcodes)

g_2_6_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_6_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
g_2_6_CD_27_barcodes <- unlist(g_2_6_CD_27_barcodes)
g_2_6_CD_27_barcodes <- gsub(pattern = "-1",
                             replacement = "-3",
                             x = g_2_6_CD_27_barcodes)
g_2_6_CD_54_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_6_common_CD27_CD54.best.CD_54", 
           col_names = FALSE)
g_2_6_CD_54_barcodes <- unlist(g_2_6_CD_54_barcodes)
g_2_6_CD_54_barcodes <- gsub(pattern = "-1",
                             replacement = "-3",
                             x = g_2_6_CD_54_barcodes)

# import group 8 cell line barcodes
g_8_0_CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_0_common_CD08_CD25_CD26.best.CD_08", 
           col_names = FALSE)
g_8_0_CD_08_barcodes <- unlist(g_8_0_CD_08_barcodes)
g_8_0_CD_08_barcodes <- gsub(pattern = "-1",
                             replacement = "-4",
                             x = g_8_0_CD_08_barcodes)


g_8_0_CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_0_common_CD08_CD25_CD26.best.CD_25", 
           col_names = FALSE)
g_8_0_CD_25_barcodes <- unlist(g_8_0_CD_25_barcodes)
g_8_0_CD_25_barcodes <- gsub(pattern = "-1",
                             replacement = "-4",
                             x = g_8_0_CD_25_barcodes)
g_8_0_CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_0_common_CD08_CD25_CD26.best.CD_26", 
           col_names = FALSE)
g_8_0_CD_26_barcodes <- unlist(g_8_0_CD_26_barcodes)
g_8_0_CD_26_barcodes <- gsub(pattern = "-1",
                             replacement = "-4",
                             x = g_8_0_CD_26_barcodes)

g_8_1_CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_1_common_CD08_CD25_CD26.best.CD_08", 
           col_names = FALSE)
g_8_1_CD_08_barcodes <- unlist(g_8_1_CD_08_barcodes)
g_8_1_CD_08_barcodes <- gsub(pattern = "-1",
                             replacement = "-5",
                             x = g_8_1_CD_08_barcodes)
g_8_1_CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_1_common_CD08_CD25_CD26.best.CD_25", 
           col_names = FALSE)
g_8_1_CD_25_barcodes <- unlist(g_8_1_CD_25_barcodes)
g_8_1_CD_25_barcodes <- gsub(pattern = "-1",
                             replacement = "-5",
                             x = g_8_1_CD_25_barcodes)

g_8_1_CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_1_common_CD08_CD25_CD26.best.CD_26", 
           col_names = FALSE)
g_8_1_CD_26_barcodes <- unlist(g_8_1_CD_26_barcodes)
g_8_1_CD_26_barcodes <- gsub(pattern = "-1",
                             replacement = "-5",
                             x = g_8_1_CD_26_barcodes)

g_8_6_CD_08_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_6_common_CD08_CD25_CD26.best.CD_08", 
           col_names = FALSE)
g_8_6_CD_08_barcodes <- unlist(g_8_6_CD_08_barcodes)
g_8_6_CD_08_barcodes <- gsub(pattern = "-1",
                             replacement = "-6",
                             x = g_8_6_CD_08_barcodes)

g_8_6_CD_25_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_6_common_CD08_CD25_CD26.best.CD_25", 
           col_names = FALSE)
g_8_6_CD_25_barcodes <- unlist(g_8_6_CD_25_barcodes)
g_8_6_CD_25_barcodes <- gsub(pattern = "-1",
                             replacement = "-6",
                             x = g_8_6_CD_25_barcodes)

g_8_6_CD_26_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_6_common_CD08_CD25_CD26.best.CD_26", 
           col_names = FALSE)
g_8_6_CD_26_barcodes <- unlist(g_8_6_CD_26_barcodes)
g_8_6_CD_26_barcodes <- gsub(pattern = "-1",
                             replacement = "-6",
                             x = g_8_6_CD_26_barcodes)

# assign cell line identity
obj_complete$cell.line.ident <- NA
obj_complete$cell.line.ident[obj_complete@assays$RNA@counts@Dimnames[[2]] %in% 
                             c(g_2_0_CD_27_barcodes,
                               g_2_1_CD_27_barcodes,
                               g_2_6_CD_27_barcodes)] <- "CD_27"
obj_complete$cell.line.ident[obj_complete@assays$RNA@counts@Dimnames[[2]] %in% 
                             c(g_2_0_CD_54_barcodes,
                               g_2_1_CD_54_barcodes,
                               g_2_6_CD_54_barcodes)] <- "CD_54"
obj_complete$cell.line.ident[obj_complete@assays$RNA@counts@Dimnames[[2]] %in% 
                             c(g_8_0_CD_08_barcodes,
                               g_8_1_CD_08_barcodes,
                               g_8_6_CD_08_barcodes)] <- "CD_08"
obj_complete$cell.line.ident[obj_complete@assays$RNA@counts@Dimnames[[2]] %in% 
                             c(g_8_0_CD_25_barcodes,
                               g_8_1_CD_25_barcodes,
                               g_8_6_CD_25_barcodes)] <- "CD_25"
obj_complete$cell.line.ident[obj_complete@assays$RNA@counts@Dimnames[[2]] %in% 
                             c(g_8_0_CD_26_barcodes,
                               g_8_1_CD_26_barcodes,
                               g_8_6_CD_26_barcodes)] <- "CD_26"
# check if anything remains unassigned
sum(is.na(obj_complete$cell.line.ident))

# remove NAs
obj_complete$cell.line.ident[is.na(obj_complete$cell.line.ident)] <- "unmatched"
obj_complete$fine.cell.type[is.na(obj_complete$fine.cell.type)] <- "unknown"

# get lists of variable values
types <- unique(obj_complete$fine.cell.type)
times <- unique(obj_complete$time.ident)
lines <- unique(obj_complete$cell.line.ident)
total.len <- length(types) * length(times) * length(lines)

df <- data.frame(cell.type = rep("", total.len),
                 time.point = rep("", total.len),
                 cell.line = rep("", total.len),
                 reads = rep(0, total.len))


i = 0
for (c in types){
  print(paste0("cell type: ", c))
  for (t in times){
    print(paste0("time: ", t))
    for (l in lines){
      print(paste0("cell line: ", l))
      i = i + 1
      nread <- sum(obj_complete@assays$peaks@counts@x[obj_complete$fine.cell.type == c &
                                                      obj_complete$time.ident == t &
                                                      obj_complete$cell.line.ident == l])
      print(nread)
      df[i, 1] <- c
      df[i, 2] <- t
      df[i, 3] <- l
      df[i, 4] <- nread
    }
  }
}

# output df
write.table(df, file = "ATACseq_reads_count.csv",
            quote = F, sep = ",", col.names = c("cell.type", "time.point", "cell.line", "counts"))

# plot it
ggplot(data = df, aes(x = cell.line, 
                      y = cell.type,
                      fill = time.point,
                      size = reads,
                      group = time.point)) +
  scale_fill_manual(values = brewer.pal(3, "Set2")) +
  geom_point(shape = 21,
             na.rm = TRUE,
             position = position_dodge(width = 0.9)) +
  scale_size_area() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  coord_flip()
