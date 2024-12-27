# Chuxuan Li 09/29/2021
# Look at sample composition for each of the clusters in the combined sample

# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)
library(future)

library(RColorBrewer)

combined_gex <- readRDS("./DESeq2/combined_g2_g8_USE_THIS.rds")

seurat_obj <- RunPCA(combined_gex,
                     seed.use = 42,
                     verbose = T)
seurat_obj <- RunUMAP(seurat_obj,
                      dims = 1:30,
                      verbose = T)
seurat_obj <- FindNeighbors(seurat_obj,
                            dims = 1:30,
                            verbose = T)
seurat_obj <- FindClusters(seurat_obj,
                           verbose = T, 
                           resolution = 0.13,
                           random.seed = 42)


length(levels(seurat_obj@meta.data$seurat_clusters))
# get the list of clusters, time points, cell lines
clusters <- levels(seurat_obj$seurat_clusters)
times <- unique(seurat_obj$time.ident)
lines <- unique(seurat_obj$cell.line.ident)
lines <- lines[!(is.na(lines))]

nc <- length(clusters)
nt <- length(times)
nl <- length(lines)

# plot a Dimplot first
DimPlot(seurat_obj, 
        label = T) +
  NoLegend()
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
      s <- sum((combined_gex$time.ident == t) & (combined_gex$cell.line.ident %in% l) & (combined_gex$seurat_clusters == c))
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
