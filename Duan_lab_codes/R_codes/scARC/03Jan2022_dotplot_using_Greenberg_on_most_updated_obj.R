# Chuxuan Li 01/03/2022
# calculating fold change and plotting for BDNF based on Siwei 12/24/2021 using
#most_updated_RNA_ATAC_combined_obj_with_all_labels.RData

library(Seurat)
library(Signac)
library(future)
library(sctransform)
library(glmGamPoi)

library(ggplot2)
library(RColorBrewer)

# setup multithread
plan("multisession", workers = 1) # should not use "multicore" here
options(future.globals.maxSize = 1474836480)
set.seed(42)
options(future.seed = T)

## use "integrate" object
DefaultAssay(aggr_filtered)
DefaultAssay(aggr_filtered) <- "SCT"

sort(unique(aggr_filtered$seurat_clusters))

unique(Idents(aggr_filtered))

i <- 0

# calc hour 1 vs 0
# save.image(file = "raw_data_with_SCT.RData")

for (i in sort(unique(aggr_filtered$seurat_clusters))) {
  print(i)
  if (i == "0" ) {
    print("first cluster")
    df_temp <- FindMarkers(object = aggr_filtered,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- i
    df_1vs0 <- df_temp[df_temp$Gene_Symbol %in% "BDNF", ]
  } else {
    print("B")
    df_temp <- FindMarkers(object = aggr_filtered,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- i
    df_to_append <- df_temp[df_temp$Gene_Symbol %in% "BDNF", ]
    df_1vs0 <- rbind(df_4vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
  }
}

# calculate 6hr vs 0hr
for (i in sort(unique(aggr_filtered$seurat_clusters))) {
  print(i)
  if (i == "0" ) {
    print("A")
    df_temp <- FindMarkers(object = aggr_filtered,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- i
    df_6vs0 <- df_temp[df_temp$Gene_Symbol %in% "BDNF", ]
  } else {
    print("B")
    df_temp <- FindMarkers(object = aggr_filtered,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- i
    df_to_append <- df_temp[df_temp$Gene_Symbol %in% "BDNF", ]
    df_6vs0 <- rbind(df_6vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
    # df_1vs0
  }
}



# make a bubble plot
df_to_plot <- rbind(df_1vs0,
                    df_6vs0)
df_to_plot$source <- c(rep_len("1hr vs 0hr", length.out = 15),
                       rep_len("6hr vs 0hr", length.out = 15))


ggplot(df_to_plot,
       aes(x = as.factor(source),
           y = as.factor(cluster),
           size = pct.1 * 100,
           fill = 2 ^ avg_log2FC * pct.1 / pct.2)) +
  geom_point(shape = 21) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  theme_classic() +
  coord_flip()

cell_count_per_cluster <- summary(SCT_iCell$seurat_clusters)
cell_count_per_cluster <- as.data.frame(cell_count_per_cluster)
cell_count_per_cluster$cluster <- rownames(cell_count_per_cluster)

df_test <- merge(df_to_plot,
                 cell_count_per_cluster,
                 by = "cluster",
                 all.x = T)

ggplot(df_test,
       aes(x = as.factor(source),
           y = as.factor(cluster),
           size = pct.1 * 100,
           fill = 2 ^ avg_log2FC * pct.1 / pct.2)) +
  geom_point(shape = 21) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  theme_classic() +
  coord_flip()  

