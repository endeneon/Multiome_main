# Chuxuan Li 02/23/2022
# plot dotplot showing the cell line/time point composition of each cluster in
#the 6-lib normalized then integrated object


# init ####

library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)
library(future)

library(RColorBrewer)

load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/27Jan2022_5line_added_cellline_idents.RData")


length(levels(RNAseq_integrated_labeled@meta.data$seurat_clusters))
# get the list of clusters, time points, cell lines
clusters <- levels(RNAseq_integrated_labeled$seurat_clusters)
times <- unique(RNAseq_integrated_labeled$time.ident)
lines <- unique(RNAseq_integrated_labeled$cell.line.ident)
lines <- lines[lines != "unmatched"]

nc <- length(clusters)
nt <- length(times)
nl <- length(lines)

# # vector to store all data
# v <- vector(mode = "numeric", length = nc * nt * nl * 4)
# initialize matrix to store number of cells
df_4_plot <- data.frame(value = rep_len(0, length.out = (nc * nt * nl)),
                        cluster = rep_len(NA, length.out = (nc * nt * nl)),
                        time = rep_len(NA, length.out = (nc * nt * nl)),
                        cell_line = rep_len(NA, length.out = (nc * nt * nl)),
                        stringsAsFactors = F)

i <- 1
for (c in clusters){
  for (t in times){
    for (l in lines){
      # print(l)
      s <- sum((RNAseq_integrated_labeled$time.ident == t) & 
                 (RNAseq_integrated_labeled$cell.line.ident %in% l) & 
                 (RNAseq_integrated_labeled$seurat_clusters == c))
      df_4_plot$value[i] <- s
      df_4_plot$cluster[i] <- c
      df_4_plot$time[i] <- t
      df_4_plot$cell_line[i] <- l
      i = i + 1
    }
  }
}

df_4_plot$cluster <- factor(df_4_plot$cluster, 
                            levels = as.character(0:nc))
#c("0",  "1",  "2",  "3",  "4",  5  6  7  8  9 10 
#11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26))
df_4_plot$time <- factor(df_4_plot$time, 
                         levels = c("6hr", "1hr", "0hr"))
df_4_plot$cell_line <- factor(df_4_plot$cell_line)

hist(df_4_plot$value)
ggplot(data = df_4_plot, aes(x = cell_line, 
                             y = cluster,
                             fill = time,
                             size = value,
                             group = time)) +
  scale_fill_manual(values = brewer.pal(3, "Dark2")) +
  geom_point(shape = 21,
             na.rm = TRUE,
             position = position_dodge(width = 0.9)) +
  scale_size_area() +
  theme_bw() +
  coord_flip()
