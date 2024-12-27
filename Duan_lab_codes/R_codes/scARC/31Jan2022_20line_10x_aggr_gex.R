# Chuxuan Li 01/31/2022
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
aggr_raw <- 
  Read10X_h5(filename = "/home/cli/Data/FASTQ/Duan_Project_022_Reseq/Duan_022_reseq_aggregated/outs/filtered_feature_bc_matrix.h5")
aggr_gex <- CreateSeuratObject(counts = aggr_raw$`Gene Expression`)
rm(aggr_raw)

aggr_gex$time.group.ident <- str_extract(aggr_gex@assays$RNA@data@Dimnames[[2]], 
                                         "[0-9]+$")

aggr_gex$time.ident <- NA
aggr_gex$time.ident[aggr_gex$time.group.ident %in% as.character(seq(1, 13, 3))] <- "0hr"
aggr_gex$time.ident[aggr_gex$time.group.ident %in% as.character(seq(2, 14, 3))] <- "1hr"  
aggr_gex$time.ident[aggr_gex$time.group.ident %in% as.character(seq(3, 15, 3))] <- "6hr"
unique(aggr_gex$time.ident)

aggr_gex <- PercentageFeatureSet(aggr_gex,
                                 pattern = c("^MT-"),
                                 col.name = "percent.mt")
# aggr_gex <- SCTransform(aggr_gex, 
#                         vars.to.regress = "percent.mt",
#                         method = "glmGamPoi",
#                         variable.features.n = 8000,
#                         seed.use = 115, 
#                         return.only.var.genes = F,
#                         verbose = T)

# Dowmstream analysis ####
aggr_gex <- FindVariableFeatures(object = aggr_gex, nfeatures = 8000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(aggr_gex), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(aggr_gex)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

aggr_gex <- NormalizeData(aggr_gex)
aggr_gex <- ScaleData(aggr_gex, assay = "RNA", vars.to.regress = "percent.mt", verbose = T)
aggr_gex <- RunPCA(aggr_gex, 
                   assay = "RNA",
                   verbose = T,  
                   seed.use = 42)
ElbowPlot(aggr_gex)
aggr_gex <- RunUMAP(aggr_gex, 
                    reduction = "pca", 
                    dims = 1:30,
                    seed.use = 42)

aggr_gex <- FindNeighbors(aggr_gex, 
                          reduction = "pca", 
                          dims = 1:30)

aggr_gex <- FindClusters(aggr_gex, 
                         resolution = 0.5,
                         random.seed = 42)

save("aggr_gex", file = "aggr_split_by_time_then_integrated_obj_inc_varfeat_to_15000.RData")

DimPlot(aggr_gex,
        repel = F, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of 10x-aggr data")

# check time point distribution 
by_time_dimplot <- DimPlot(aggr_gex,
                           #cols = c("steelblue3", "mediumpurple4", "indianred3"),
                           cols = c("cyan", "magenta", "yellow"),
                           # cols = "Set2",
                           # cols = g_2_6_gex@meta.data$cell.line.ident,
                           group.by = "time.ident",
                           label = F) +
  ggtitle("10x-aggr data split by time \n then integrated by Seurat \n colored by time")
by_time_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
by_time_dimplot

DimPlot(aggr_gex,
        #cols = c("steelblue3", "mediumpurple4", "indianred3"),
        cols = c("cyan", "transparent", "transparent"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("10x-aggr data split by time \n then integrated by Seurat \n colored by time - 0hr")

DimPlot(aggr_gex,
        #cols = c("steelblue3", "mediumpurple4", "indianred3"),
        cols = c("transparent", "magenta", "transparent"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("10x-aggr data split by time \n then integrated by Seurat \n colored by time - 1hr")

DimPlot(aggr_gex,
        #cols = c("steelblue3", "mediumpurple4", "indianred3"),
        cols = c("transparent", "transparent", "gold"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("10x-aggr data split by time \n then integrated by Seurat \n colored by time - 6hr")

# check library distribution
by_lib_dimplot <- DimPlot(aggr_gex,
                          #cols = c("steelblue3", "mediumpurple4", "indianred3"),
                          cols = "Set1",
                          #cols = c("transparent", "transparent", "transparent", "transparent", "transparent", "steelblue"),
                          # cols = g_2_6_gex@meta.data$cell.line.ident,
                          group.by = "time.group.ident",
                          label = F) +
  ggtitle("10x-aggr data split by time \n then integrated by Seurat \n colored by library")
by_lib_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
by_lib_dimplot

sum(aggr_gex$time.group.ident == "1") #g20: 14572
sum(aggr_gex$time.group.ident == "2") #g21: 13560
sum(aggr_gex$time.group.ident == "3") #g26: 14721
sum(aggr_gex$time.group.ident == "4") #g80: 13261
sum(aggr_gex$time.group.ident == "5") #g81: 12534
sum(aggr_gex$time.group.ident == "6") #g86: 14445




# assign cell types ####

DefaultAssay(aggr_gex) <- "integrated"
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

StackedVlnPlot(obj = aggr_gex,
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
StackedVlnPlot(obj = aggr_gex, features = trimmed_markers) +
  coord_flip()
StackedVlnPlot(obj = aggr_gex, features = other_markers) +
  coord_flip()


FeaturePlot(aggr_gex, features = c("GAD1", "GAD2"), max.cutoff = 10)
FeaturePlot(aggr_gex, features = c("SLC17A6", "SLC17A7"), max.cutoff = 15)
FeaturePlot(aggr_gex, features = c("CUX2", "NEFM"), max.cutoff = 15)
FeaturePlot(aggr_gex, features = c("SEMA3E"), max.cutoff = 20)
FeaturePlot(aggr_gex, features = c("FOXG1"), max.cutoff = 10)
FeaturePlot(aggr_gex, features = c("VIM"), max.cutoff = 30)
FeaturePlot(aggr_gex, features = c("SOX2"), max.cutoff = 10) 
FeaturePlot(aggr_gex, features = c("NES"), max.cutoff = 15) 


FeaturePlot(aggr_gex, features = "BDNF", cols = brewer.pal(n = 5, name = "Purples"))
# early_late_markers <- c("EGR1", "FOS", "FOSB", "NPAS4", # early response
#                         "VGF", "BDNF", "NRN1", "PNOC" # late response
# )
# StackedVlnPlot(obj = aggr_gex, features = early_late_markers) +
#   coord_flip() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


# try a dendrogram plot
aggr_gex <- BuildClusterTree(object = aggr_gex, 
                                     assay = "SCT",
                                     features = trimmed_markers,
                                     verbose = T)

PlotClusterTree(aggr_gex, 
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
