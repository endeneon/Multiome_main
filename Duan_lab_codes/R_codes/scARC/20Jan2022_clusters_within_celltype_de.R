# Chuxuan Li 01/20/2022
# Split the four_celltype_obj into individual cell types, then look at the clusters
#within these cell types and perform differential gene expression analysis.


# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)

library(readr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(graphics)
library(stringr)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)


set.seed(120)


setwd("/home/cli/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes")

load("four_celltype_obj.RData")

# split the object
unique(four_celltype_obj$fine.cell.type)
npcm <- subset(four_celltype_obj, subset = fine.cell.type == "NpCm_glut")
nmcp <- subset(four_celltype_obj, subset = fine.cell.type == "NmCp_glut")
gaba <- subset(four_celltype_obj, subset = fine.cell.type == "GABA")

# check the clusters
unique(npcm$RNA.cluster.ident)
unique(nmcp$RNA.cluster.ident)
unique(gaba$RNA.cluster.ident)
gaba$RNA.cluster.ident[is.na(gaba$RNA.cluster.ident)] <- "unknown"
gaba <- subset(gaba, subset = RNA.cluster.ident != "unknown")

# prepare the objects
Idents(npcm) <- "RNA.cluster.ident"
Idents(nmcp) <- "RNA.cluster.ident"
Idents(gaba) <- "RNA.cluster.ident"

DimPlot(npcm, group.by = "time.ident")
sum(npcm$RNA.cluster.ident == "6")
cluster6 <- subset(npcm, subset = RNA.cluster.ident == "6")
sum(cluster6$time.ident == "6hr")
unique(npcm$RNA.cluster.ident)

# npcm ####
i = 0
for (c in c("6", "1", "8", "13", "15", "16")) {
  i = i + 1
  if (i == 1) {
    print("the first cluster")
    df_temp <- FindMarkers(object = npcm,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = c,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- c
    df_1vs0 <- df_temp[df_temp$Gene_Symbol %in% c("ATP2A2", "BDNF", "EGR1",
                                                  "ETF1", "FOS", "MMP1", "NPAS4",
                                                  "NR4A1", "PPP1R13B", "RGS2",
                                                  "ZSWIM6"), ]
  } else {
    print("the rest")
    df_temp <- FindMarkers(object = npcm,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = c,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- c
    df_to_append <- df_temp[df_temp$Gene_Symbol %in% c("ATP2A2", "BDNF", "EGR1",
                                                       "ETF1", "FOS", "MMP1", "NPAS4",
                                                       "NR4A1", "PPP1R13B", "RGS2",
                                                       "ZSWIM6"), ]
    df_1vs0 <- rbind(df_1vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
  }
}


i = 0
for (c in sort(unique(npcm$RNA.cluster.ident))) {
  i = i + 1
  if (i == 1) {
    print("the first cluster")
    df_temp <- FindMarkers(object = npcm,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = c,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- c
    df_6vs0 <- df_temp[df_temp$Gene_Symbol %in% c("ATP2A2", "BDNF", "EGR1",
                                                  "ETF1", "FOS", "MMP1", "NPAS4",
                                                  "NR4A1", "PPP1R13B", "RGS2",
                                                  "ZSWIM6"), ]
  } else {
    print("the rest")
    df_temp <- FindMarkers(object = npcm,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = c,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- c
    df_to_append <- df_temp[df_temp$Gene_Symbol %in% c("ATP2A2", "BDNF", "EGR1",
                                                       "ETF1", "FOS", "MMP1", "NPAS4",
                                                       "NR4A1", "PPP1R13B", "RGS2",
                                                       "ZSWIM6"), ]
    df_6vs0 <- rbind(df_6vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
    # df_1vs0
  }
}



# make a df to store
df_npcm <- rbind(df_1vs0,
                df_6vs0)




# nmcp ####
i = 0
for (c in unique(nmcp$RNA.cluster.ident)) {
  i = i + 1
  if (i == 1) {
    print("the first cluster")
    df_temp <- FindMarkers(object = nmcp,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = c,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- c
    df_1vs0 <- df_temp[df_temp$Gene_Symbol %in% c("ATP2A2", "BDNF", "EGR1",
                                                  "ETF1", "FOS", "MMP1", "NPAS4",
                                                  "NR4A1", "PPP1R13B", "RGS2",
                                                  "ZSWIM6"), ]
  } else {
    print("the rest")
    df_temp <- FindMarkers(object = nmcp,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = c,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- c
    df_to_append <- df_temp[df_temp$Gene_Symbol %in% c("ATP2A2", "BDNF", "EGR1",
                                                       "ETF1", "FOS", "MMP1", "NPAS4",
                                                       "NR4A1", "PPP1R13B", "RGS2",
                                                       "ZSWIM6"), ]
    df_1vs0 <- rbind(df_1vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
  }
}


i = 0
for (c in sort(unique(nmcp$RNA.cluster.ident))) {
  i = i + 1
  if (i == 1) {
    print("the first cluster")
    df_temp <- FindMarkers(object = nmcp,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = c,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- c
    df_6vs0 <- df_temp[df_temp$Gene_Symbol %in% c("ATP2A2", "BDNF", "EGR1",
                                                  "ETF1", "FOS", "MMP1", "NPAS4",
                                                  "NR4A1", "PPP1R13B", "RGS2",
                                                  "ZSWIM6"), ]
  } else {
    print("the rest")
    df_temp <- FindMarkers(object = nmcp,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = c,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- c
    df_to_append <- df_temp[df_temp$Gene_Symbol %in% c("ATP2A2", "BDNF", "EGR1",
                                                       "ETF1", "FOS", "MMP1", "NPAS4",
                                                       "NR4A1", "PPP1R13B", "RGS2",
                                                       "ZSWIM6"), ]
    df_6vs0 <- rbind(df_6vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
    # df_1vs0
  }
}

# make a df to store
df_nmcp <- rbind(df_1vs0,
                 df_6vs0)



# gaba ####
i = 0
for (c in unique(gaba$RNA.cluster.ident)) {
  i = i + 1
  if (i == 1) {
    print("the first cluster")
    df_temp <- FindMarkers(object = gaba,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = c,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- c
    df_1vs0 <- df_temp[df_temp$Gene_Symbol %in% c("ATP2A2", "BDNF", "EGR1",
                                                  "ETF1", "FOS", "MMP1", "NPAS4",
                                                  "NR4A1", "PPP1R13B", "RGS2",
                                                  "ZSWIM6"), ]
  } else {
    print("the rest")
    df_temp <- FindMarkers(object = gaba,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = c,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- c
    df_to_append <- df_temp[df_temp$Gene_Symbol %in% c("ATP2A2", "BDNF", "EGR1",
                                                       "ETF1", "FOS", "MMP1", "NPAS4",
                                                       "NR4A1", "PPP1R13B", "RGS2",
                                                       "ZSWIM6"), ]
    df_1vs0 <- rbind(df_1vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
  }
}


i = 0
for (c in sort(unique(gaba$RNA.cluster.ident))) {
  i = i + 1
  if (i == 1) {
    print("the first cluster")
    df_temp <- FindMarkers(object = gaba,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = c,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- c
    df_6vs0 <- df_temp[df_temp$Gene_Symbol %in% c("ATP2A2", "BDNF", "EGR1",
                                                  "ETF1", "FOS", "MMP1", "NPAS4",
                                                  "NR4A1", "PPP1R13B", "RGS2",
                                                  "ZSWIM6"), ]
  } else {
    print("the rest")
    df_temp <- FindMarkers(object = gaba,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = c,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_temp$Gene_Symbol <- rownames(df_temp)
    df_temp$cluster <- c
    df_to_append <- df_temp[df_temp$Gene_Symbol %in% c("ATP2A2", "BDNF", "EGR1",
                                                       "ETF1", "FOS", "MMP1", "NPAS4",
                                                       "NR4A1", "PPP1R13B", "RGS2",
                                                       "ZSWIM6"), ]
    df_6vs0 <- rbind(df_6vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
    # df_1vs0
  }
}


# make a df to store
df_gaba <- rbind(df_1vs0,
                 df_6vs0)


# plotting ####
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

