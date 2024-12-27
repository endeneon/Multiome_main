# Chuxuan Li 01/06/2022
# using 10x-aggr produced dataset - RNAseq portion, then check differential
#expression

# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)

library(readr)
library(stringr)
library(dplyr)

library(ggplot2)
library(RColorBrewer)
library(viridis)
library(graphics)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)

library(future)

set.seed(1886)
# set threads and parallelization
plan("multisession", workers = 1)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 11474836480)



# read files and normalization ####

# read in the combined object obtained by the -aggr method from 10x
aggr_raw <- Read10X_h5(filename = "/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/filtered_feature_bc_matrix.h5")

aggr_gex <- CreateSeuratObject(counts = aggr_raw$`Gene Expression`)
rm(aggr_raw)

aggr_gex$time.group.ident <- str_sub(aggr_gex@assays$RNA@data@Dimnames[[2]], 
                                 start = -1L)

aggr_gex$time.ident <- NA
aggr_gex$time.ident[aggr_gex$time.group.ident %in% c("1", "4")] <- "0hr"
aggr_gex$time.ident[aggr_gex$time.group.ident %in% c("2", "5")] <- "1hr"  
aggr_gex$time.ident[aggr_gex$time.group.ident %in% c("3", "6")] <- "6hr"

aggr_gex <- PercentageFeatureSet(aggr_gex,
                                 pattern = c("^MT-"),
                                 col.name = "percent.mt")
aggr_gex <- SCTransform(aggr_gex, 
                        vars.to.regress = "percent.mt",
                        method = "glmGamPoi",
                        variable.features.n = 8000,
                        seed.use = 115, 
                        return.only.var.genes = F,
                        verbose = T)

# split by time then integrate ####
t0 <- subset(aggr_gex, subset = time.ident == "0hr")
t1 <- subset(aggr_gex, subset = time.ident == "1hr")
t6 <- subset(aggr_gex, subset = time.ident == "6hr")
obj.lst <- c(t0, t1, t6)

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = obj.lst, 
                                      nfeatures = 15000,
                                      fvf.nfeatures = 15000)
anchors <- FindIntegrationAnchors(object.list = obj.lst, anchor.features = features)
aggr_ts_combined <- IntegrateData(anchorset = anchors)


# prelim downstream analysis ####
DefaultAssay(aggr_ts_combined) <- "integrated"

aggr_ts_combined <- FindVariableFeatures(object = aggr_ts_combined, nfeatures = )
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(aggr_ts_combined), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(aggr_ts_combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(aggr_gex)
aggr_ts_combined <- ScaleData(aggr_ts_combined, features = all.genes)

aggr_ts_combined <- RunPCA(aggr_ts_combined, 
                           assay = "integrated",
                           verbose = T, 
                           seed.use = 42)

aggr_ts_combined <- RunUMAP(aggr_ts_combined, 
                            reduction = "pca", 
                            dims = 1:30,
                            seed.use = 42)

aggr_ts_combined <- FindNeighbors(aggr_ts_combined, 
                                  reduction = "pca", 
                                  dims = 1:30)

aggr_ts_combined <- FindClusters(aggr_ts_combined, 
                                 resolution = 0.5,
                                 random.seed = 42)

save("aggr_ts_combined", file = "aggr_split_by_time_then_integrated_obj_inc_varfeat_to_15000.RData")

DimPlot(aggr_ts_combined,
        repel = F, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of 10x-aggr data")

# check time point distribution 
by_time_dimplot <- DimPlot(aggr_ts_combined,
                              #cols = c("steelblue3", "mediumpurple4", "indianred3"),
                              cols = c("cyan", "magenta", "yellow"),
                              # cols = "Set2",
                              # cols = g_2_6_gex@meta.data$cell.line.ident,
                              group.by = "time.ident",
                              label = F) +
  ggtitle("10x-aggr data split by time \n then integrated by Seurat \n colored by time")
by_time_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
by_time_dimplot

DimPlot(aggr_ts_combined,
        #cols = c("steelblue3", "mediumpurple4", "indianred3"),
        cols = c("cyan", "transparent", "transparent"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("10x-aggr data split by time \n then integrated by Seurat \n colored by time - 0hr")

DimPlot(aggr_ts_combined,
        #cols = c("steelblue3", "mediumpurple4", "indianred3"),
        cols = c("transparent", "magenta", "transparent"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("10x-aggr data split by time \n then integrated by Seurat \n colored by time - 1hr")

DimPlot(aggr_ts_combined,
        #cols = c("steelblue3", "mediumpurple4", "indianred3"),
        cols = c("transparent", "transparent", "gold"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("10x-aggr data split by time \n then integrated by Seurat \n colored by time - 6hr")

# check library distribution
by_lib_dimplot <- DimPlot(aggr_ts_combined,
                          #cols = c("steelblue3", "mediumpurple4", "indianred3"),
                          cols = "Set1",
                          #cols = c("transparent", "transparent", "transparent", "transparent", "transparent", "steelblue"),
                          # cols = g_2_6_gex@meta.data$cell.line.ident,
                          group.by = "time.group.ident",
                          label = F) +
  ggtitle("10x-aggr data split by time \n then integrated by Seurat \n colored by library")
by_lib_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
by_lib_dimplot

sum(aggr_ts_combined$time.group.ident == "1") #g20: 14572
sum(aggr_ts_combined$time.group.ident == "2") #g21: 13560
sum(aggr_ts_combined$time.group.ident == "3") #g26: 14721
sum(aggr_ts_combined$time.group.ident == "4") #g80: 13261
sum(aggr_ts_combined$time.group.ident == "5") #g81: 12534
sum(aggr_ts_combined$time.group.ident == "6") #g86: 14445




# assign cell types ####

DefaultAssay(aggr_ts_combined) <- "integrated"
#plot stacked vlnplot
subtype_markers <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
                     "EBF1", "SERTAD4", # striatal
                     "FOXG1", # forebrain
                     "TBR1", "FOXP2", "TLE4", # pallial glutamatergic
                     "FEZF2", "ADCYAP1", "TCERG1L", "SEMA3E", # subcerebral
                     "POU3F2", "CUX1", "BCL11B",  # cortical
                     "LHX2", # general cortex
                     #"EOMES", # hindbrain
                     "NPY", "SST", "DLX2", "DLX5", # inhibitory
                     "SATB2", "CBLN2", "CUX2", "NEFM", # excitatory
                     "VIM", "SOX2", "NES" #NPC
)

StackedVlnPlot(obj = aggr_ts_combined,
               features = subtype_markers) +
  coord_flip()

# more stacked vlnplots for late/early response genes
trimmed_markers <- c("GAD1", "GAD2", "SLC17A6",
                     "EBF1", # striatal
                     "ADCYAP1", "TCERG1L", "SEMA3E", # subcerebral
                     "CUX1", "BCL11B",  # cortical
                     "SST", # inhibitory
                     "SATB2", "CUX2", "NEFM", # excitatory
                     "VIM", "SOX2", "NES", #NPC
                     "MAP2", "DCX"
)

other_markers <- c("SLC17A7", "SERTAD4", "FOXG1", # forebrain
                   "POU3F2", "LHX2" # general cortex
)
StackedVlnPlot(obj = aggr_ts_combined, features = trimmed_markers) +
  coord_flip()
StackedVlnPlot(obj = aggr_ts_combined, features = other_markers) +
  coord_flip()


FeaturePlot(aggr_ts_combined, features = c("GAD1", "GAD2"), max.cutoff = 10)
FeaturePlot(aggr_ts_combined, features = c("SLC17A6", "SLC17A7"), max.cutoff = 15)
FeaturePlot(aggr_ts_combined, features = c("CUX2", "NEFM"), max.cutoff = 15)
FeaturePlot(aggr_ts_combined, features = c("SEMA3E"), max.cutoff = 20)
FeaturePlot(aggr_ts_combined, features = c("FOXG1"), max.cutoff = 10)
FeaturePlot(aggr_ts_combined, features = c("VIM"), max.cutoff = 30)
FeaturePlot(aggr_ts_combined, features = c("SOX2"), max.cutoff = 10) 
FeaturePlot(aggr_ts_combined, features = c("NES"), max.cutoff = 15) 


FeaturePlot(aggr_ts_combined, features = "BDNF", cols = brewer.pal(n = 5, name = "Purples"))
# early_late_markers <- c("EGR1", "FOS", "FOSB", "NPAS4", # early response
#                         "VGF", "BDNF", "NRN1", "PNOC" # late response
# )
# StackedVlnPlot(obj = aggr_ts_combined, features = early_late_markers) +
#   coord_flip() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


# try a dendrogram plot
aggr_ts_combined <- BuildClusterTree(object = aggr_ts_combined, 
                               assay = "SCT",
                               features = trimmed_markers,
                               verbose = T)

PlotClusterTree(aggr_ts_combined, 
                direction = "rightwards",
                use.edge.length = F)

# assign cluster identities based on the violin plots
new.cluster.ids <- 
  c("NEFM_neg_glut", "NEFM_pos_glut", "GABA", "forebrain_NPC", "GABA", 
    "GABA", "NEFM_pos_glut", "SEMA3E_pos_glut", "NEFM_neg_glut", "unknown", 
    "unknown", "NPC", "NEFM_pos_glut", "NPC", "NEFM_pos_glut", 
    "SEMA3E_pos_glut", "NEFM_pos_glut", "forebrain_NPC", "NEFM_neg_glut", "NEFM_neg_glut",
    "GABA", "unknown", "forebrain_NPC", "NEFM_pos_glut", "NEFM_pos_glut")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(aggr_labeled$seurat_clusters))
names(new.cluster.ids) <- levels(aggr_labeled)

aggr_labeled <- RenameIdents(aggr_labeled, new.cluster.ids)

DimPlot(aggr_labeled, 
        reduction = "umap", 
        label = TRUE, 
        repel = T, 
        pt.size = 0.3) +
  ggtitle("independently normalized RNAseq data \n integrated by Seurat, labeled by cell type")

sum(aggr_labeled@active.ident == "forebrain_NPC") #7766
sum(aggr_labeled@active.ident == "NEFM_neg_glut") #16769
sum(aggr_labeled@active.ident == "NEFM_pos_glut") #21248
sum(aggr_labeled@active.ident == "SEMA3E_pos_glut") #5537
sum(aggr_labeled@active.ident == "GABA") #19714
sum(aggr_labeled@active.ident == "NPC") #5333


save("aggr_labeled", file = "aggr_labeled_inc_varfeat.RData")


# differentially expressed genes ####

aggr_labeled$broad.cell.type <- NA
aggr_labeled$broad.cell.type[aggr_labeled@active.ident == "subcerebral_neuron"] <- "subcerebral_neuron"
aggr_labeled$broad.cell.type[aggr_labeled@active.ident == "forebrain_NPC"] <- "forebrain_NPC"
aggr_labeled$broad.cell.type[aggr_labeled@active.ident == "NEFM_pos_glut"] <- "NEFM_pos_glut"
aggr_labeled$broad.cell.type[aggr_labeled@active.ident == "NEFM_neg_glut"] <- "NEFM_neg_glut"
aggr_labeled$broad.cell.type[aggr_labeled@active.ident == "NPC"] <- "NPC"
aggr_labeled$broad.cell.type[aggr_labeled@active.ident == "GABA"] <- "GABA"
unique(aggr_labeled$broad.cell.type)

# find variable markers for each time point
GABA_subset <- subset(aggr_labeled, subset = broad.cell.type == "GABA")
NPCM_pos_glut_subset <- subset(aggr_labeled, subset = broad.cell.type == "NEFM_pos_glut")
NPCM_neg_glut_subset <- subset(aggr_labeled, subset = broad.cell.type == "NEFM_neg_glut")
NPC_subset <- subset(aggr_labeled, subset = broad.cell.type == "NPC")
forebrain_NPC_subset <- subset(aggr_labeled, subset = broad.cell.type == "forebrain_NPC")


for (t in c(GABA_subset, NPCM_pos_glut_subset, NPCM_neg_glut_subset, NPC_subset, forebrain_NPC_subset)){
  t_name <- unique(t$broad.cell.type)
  print(t_name)
  
  # differentially expressed genes comparing one big cell type with all others
  Idents(t) <- "time.ident"
  print(DefaultAssay(t))
  zero_one <- FindMarkers(t, 
                          ident.1 = "1hr", 
                          logfc.threshold = 0, 
                          min.pct = 0.1,
                          ident.2 = "0hr", 
                          group.by = "time.ident",
                          test.use = "LR",
                          verbose = TRUE)
  zero_six <- FindMarkers(t, 
                          ident.1 = "6hr", 
                          logfc.threshold = 0,
                          min.pct = 0.1,
                          ident.2 = "0hr", 
                          group.by = "time.ident",
                          test.use = "LR",
                          verbose = TRUE)
  
  # output the list of variable markers 
  write.table(zero_one, 
              file = paste0(t_name, "_DE_markers_1v0hr.csv"),
              col.names = T,
              row.names = T,
              sep = ",",
              quote = F)
  print(nrow(zero_one))
  write.table(zero_six, 
              file = paste0(t_name, "_DE_markers_6v0hr.csv"),
              col.names = T,
              row.names = T,
              sep = ",",
              quote = F)
  
}

# DefaultAssay(aggr_labeled) <- "SCT"
# DPnew2(aggr_labeled, 
#        features = early_late_markers, 
#        cols = c("royalblue3", "mediumpurple3", "red3"), 
#        dot.scale = 6, 
#        split.by = "time.ident") +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 60, 
#                                    hjust = 1), validate = TRUE) +
#   coord_flip() 
# 
# # get an object with only glut and GABA cells
# glut_GABA_only_obj <- subset(aggr_labeled, subset = broad.cell.type %in% c("NEFM_pos_glut", "NEFM_neg_glut", "GABA"))
# DefaultAssay(glut_GABA_only_obj) <- "SCT"
# DPnew2(glut_GABA_only_obj, 
#        features = early_late_markers, 
#        cols = c("royalblue3", "mediumpurple3", "red3"), 
#        dot.scale = 6, 
#        split.by = "time.ident") +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 60, 
#                                    hjust = 1), validate = TRUE) +
#   coord_flip() 


# use new method to plot dotplot
## use "integrate" object
DefaultAssay(aggr_labeled) <- "integrated"
sort(unique(aggr_labeled$broad.cell.type))
unique(Idents(aggr_labeled))
Idents(aggr_labeled) <- aggr_labeled$broad.cell.type
unique(aggr_labeled$time.ident)

FeaturePlot(object = aggr_labeled, features = "FOS")

i <- 0
# calc hour 1 vs 0
for (i in sort(unique(aggr_labeled$broad.cell.type))) {
  print(i)
  if (i == "forebrain_NPC" ) {
    print("first cell type")
    df_1vs0 <- FindMarkers(object = aggr_labeled,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_1vs0$Gene_Symbol <- rownames(df_1vs0)
    df_1vs0$cluster <- i
    #df_1vs0 <- df_temp[df_temp$Gene_Symbol %in% "BDNF", ]
  } else {
    print("others")
    df_to_append <- FindMarkers(object = aggr_labeled,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_to_append$Gene_Symbol <- rownames(df_to_append)
    df_to_append$cluster <- i
    #df_to_append <- df_1vs0[df_1vs0$Gene_Symbol %in% "BDNF", ]
    df_1vs0 <- rbind(df_1vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
  }
}

i <- 0
# calculate 6hr vs 0hr
for (i in sort(unique(aggr_labeled$broad.cell.type))) {
  print(i)
  if (i == "forebrain_NPC" ) {
    print("first cell type")
    df_6vs0 <- FindMarkers(object = aggr_labeled,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_6vs0$Gene_Symbol <- rownames(df_6vs0)
    df_6vs0$cluster <- i
    #df_6vs0 <- df_temp[df_temp$Gene_Symbol %in% "BDNF", ]
  } else {
    df_to_append <- FindMarkers(object = aggr_labeled,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_to_append$Gene_Symbol <- rownames(df_to_append)
    df_to_append$cluster <- i
    #df_to_append <- df_temp[df_temp$Gene_Symbol %in% "BDNF", ]
    df_6vs0 <- rbind(df_6vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
  }
}



# make a bubble plot

DefaultAssay(aggr_labeled)
hist(df_1vs0$pct.2, breaks = 10)

df_to_plot <- rbind(df_1vs0[df_1vs0$Gene_Symbol == "BDNF", ],
                    df_6vs0[df_6vs0$Gene_Symbol == "BDNF", ])
df_to_plot$source <- c(rep_len("1v0hr", length.out = 5),
                       rep_len("6v0hr", length.out = 5))

df_to_plot$cluster <- unlist(as.character(df_to_plot$cluster))


df_to_plot$avg_log2FC[df_to_plot$avg_log2FC > 5] <- 0
ggplot(df_to_plot,
       aes(x = as.factor(source),
           y = as.factor(cluster),
           size = pct.1 * 100,
           #fill = 2 ^ avg_log2FC * pct.1 / pct.2)
           fill = avg_log2FC)) +
  geom_point(shape = 21) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  theme_light() #+
#coord_flip()

df_to_plot_filtered <- df_to_plot
df_to_plot_filtered$avg_log2FC[df_to_plot_filtered$p_val_adj > 0.05] <- 0
df_to_plot_filtered$pct.1[df_to_plot_filtered$p_val_adj > 0.05] <- 0
df_to_plot_filtered$pct.2[df_to_plot_filtered$pct.2 == 0] <- 1

df_to_plot_filtered$cluster <- factor(df_to_plot_filtered$cluster, 
                                      levels = rev(unique(as.character(df_to_plot_filtered$cluster))))

ggplot(df_to_plot_filtered,
       aes(x = as.factor(source),
           y = as.factor(cluster),
           size = pct.1 * 100,
           #fill = 2 ^ avg_log2FC * pct.1 / pct.2)
           fill = avg_log2FC)) +
  geom_point(shape = 21) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  theme_light()

