# Siwei 18 Nov 2024
# Re-analyse Lexi's results using the grouping info investigated

# Grouping info: ####
# 20080-0: sample 1 (Danii), anesthesia;
# 20088-1: sample 8 (Danii), control;
# 20088-6: sample 3 (Danii), anesthesia;
# 60060-0: sample 2 (Danii), anesthesia;
# 60060-1: sample 0 (Danii), control;
# 60060-6: sample 9 (Danii), control.

# init ####
{
  library(Seurat)
  library(ggplot2)
  library(pals)
  library(stringr)
  library(future)
  
  library(readr)
  
  library(harmony)
  
  library(EnhancedVolcano)
}

plan("multisession", workers = 4)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

load("integrated_labeled_after_filtering2.RData")


integrated_labelled_filtered_raw_RNA_only <-
  readRDS("integrated_labelled_filtered_raw_RNA_only_oc2_18Nov2024.RDs")


# practise with Harmony
# data("cell_lines")
head(integrated_labelled_filtered_raw_RNA_only@meta.data)

# Assign control/anesthesia conditions
DefaultAssay(integrated_labelled_filtered_raw_RNA_only) <- "RNA"
integrated_labelled_filtered_raw_RNA_only$treatment <- NA
integrated_labelled_filtered_raw_RNA_only$treatment[integrated_labelled_filtered_raw_RNA_only$orig.ident %in%
                                                      c("20088-1",
                                                        "60060-1",
                                                        "60060-6")] <-
  "control"
integrated_labelled_filtered_raw_RNA_only$treatment[integrated_labelled_filtered_raw_RNA_only$orig.ident %in%
                                                      c("20088-0",
                                                        "20088-6",
                                                        "60060-0")] <-
  "anesthesia"
sum(is.nan(integrated_labelled_filtered_raw_RNA_only$treatment))

# Re-calculate embeddings
integrated_oc2_analysis <-
  NormalizeData(integrated_labelled_filtered_raw_RNA_only,
                normalization.method = "LogNormalize",
                verbose = T)
integrated_oc2_analysis <-
  FindVariableFeatures(integrated_oc2_analysis,
                       selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)
integrated_oc2_analysis <-
  integrated_oc2_analysis %>%
  ScaleData(verbose = T) %>%
  RunPCA(verbose = T)

integrated_oc2_analysis$cell.type <- NULL

Idents(integrated_oc2_analysis) <- "treatment"


Harmony_integrated_oc2 <-
  RunHarmony(integrated_oc2_analysis,
             group.by.vars = "group.ident",
             plot_convergence = T)

Harmony_integrated_oc2 <-
  RunUMAP(Harmony_integrated_oc2,
          reduction = "harmony",
          dims = 1:20)

Harmony_integrated_oc2 <-
  Harmony_integrated_oc2 %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.5,
               random.seed = 42)

DimPlot(Harmony_integrated_oc2,
        reduction = "umap",
        group.by = "treatment")
DimPlot(Harmony_integrated_oc2,
        reduction = "umap",
        group.by = "RNA_snn_res.0.5")

Harmony_integrated_oc2 <-
  FindClusters(Harmony_integrated_oc2,
               resolution = 0.1,
               random.seed = 42)
DimPlot(Harmony_integrated_oc2,
        reduction = "umap",
        group.by = "RNA_snn_res.0.1")

# make a lookup table between oc2 ensemblid and gene_symbol #####
oc2_geneid_symbol_lookup <-
  read_delim("~/Data/Databases/Genomes/CellRanger_10x/raw_GRCh38_mm10_oc2/Orycun2-2020-A-build/oryCun2_geneid_gene_symbol_list.txt",
             delim = "\t")
oc2_geneid_symbol_lookup <-
  oc2_geneid_symbol_lookup[!(is.na(oc2_geneid_symbol_lookup$`Gene name`)), ]


oc2_geneid_symbol_lookup[oc2_geneid_symbol_lookup$`Gene name` == "GFAP", ]


# ENSOCUG00000016389 SLC17A6 
# ENSOCUG00000012799 SLC17A7
# ENSOCUG00000012644 GAD1 
# ENSOCUG00000017480 GAD2
# ENSOCUG00000006895 GFAP

FeaturePlot(Harmony_integrated_oc2,
            features = c("ENSOCUG00000016389"))
FeaturePlot(Harmony_integrated_oc2,
            features = c("ENSOCUG00000012644"))

# cut out neurons (cluster 0:2) and run FindVariables again

Harmony_neurons_only <-
  Harmony_integrated_oc2[, Harmony_integrated_oc2$RNA_snn_res.0.1 %in% c(0:2)]

Harmony_neurons_only <-
  Harmony_neurons_only %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 100,
                       verbose = T) %>%
  ScaleData(verbose = T) %>%
  RunPCA(verbose = T)

Harmony_neurons_only <-
  RunUMAP(Harmony_neurons_only,
          reduction = "harmony",
          seed.use = 42, umap.method = "uwot",
          repulsion.strength = 0.5,
          dims = 1:30)
DimPlot(Harmony_neurons_only,
        reduction = "umap",
        group.by = "treatment")
FeaturePlot(Harmony_neurons_only,
            features = c("ENSOCUG00000016389"))
FeaturePlot(Harmony_neurons_only,
            features = c("ENSOCUG00000012644"))

Idents(Harmony_neurons_only) <- "RNA_snn_res.0.5"
DimPlot(Harmony_neurons_only)

hist(Harmony_neurons_only$nFeature_RNA)


Harmony_neurons_only <-
  Harmony_neurons_only[, !(Harmony_neurons_only$RNA_snn_res.0.5 %in% c(9, 14))]
Harmony_neurons_only <-
  Harmony_neurons_only[, !(Harmony_neurons_only$RNA_snn_res.0.5 %in% c(7))]


Harmony_neurons_only <-
  Harmony_neurons_only %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 500,
                       verbose = T) %>%
  ScaleData(verbose = T) %>%
  RunPCA(verbose = T)

Harmony_neurons_only <-
  RunUMAP(Harmony_neurons_only,
          reduction = "harmony",
          seed.use = 42, umap.method = "uwot",
          repulsion.strength = 0.5,
          dims = 1:30)

Harmony_neurons_only <-
  Harmony_neurons_only %>%
  FindNeighbors() %>%
  FindClusters(resolution = 1,
               random.seed = 42)
Harmony_neurons_only <-
  Harmony_neurons_only %>%
  FindNeighbors() %>%
  FindClusters(resolution = 0.4,
               random.seed = 42)

Idents(Harmony_neurons_only) <- "RNA_snn_res.1"
DimPlot(Harmony_neurons_only)

Idents(Harmony_neurons_only) <- "RNA_snn_res.0.4"
DimPlot(Harmony_neurons_only)

Harmony_neurons_glut <-
  Harmony_neurons_only[, Harmony_neurons_only$RNA_snn_res.0.4 %in% c(0, 1, 4)]
DimPlot(Harmony_neurons_glut)


unique(Harmony_neurons_glut$treatment)
Idents(Harmony_neurons_glut) <- "treatment"

# up-regulated in anesthesia
Glut_anes_markers <-
  FindMarkers(Harmony_neurons_glut,
              ident.1 = "anesthesia", 
              logfc.threshold = 0,
              test.use = "MAST",
              verbose = T,
              random.seed = 42)

EnhancedVolcano(toptable = Glut_anes_markers,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                lab = rownames(Glut_anes_markers),
                pCutoff = 0.05,
                FCcutoff = 0.5)

write.table(Glut_anes_markers,
            file = "Glut_anes_markers_MAST.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)

FeaturePlot(Harmony_neurons_only,
            features = c("ENSOCUG00000016389"))
FeaturePlot(Harmony_neurons_only,
            features = c("ENSOCUG00000012644"))


Harmony_neurons_GABA <-
  Harmony_neurons_only[, !(Harmony_neurons_only$RNA_snn_res.0.4 %in% c(0, 1, 4, 6, 7))]
DimPlot(Harmony_neurons_GABA)
Idents(Harmony_neurons_GABA) <- "treatment"

# up-regulated in anesthesia
GABA_anes_markers <-
  FindMarkers(Harmony_neurons_GABA,
              ident.1 = "anesthesia", 
              logfc.threshold = 0,
              test.use = "MAST",
              verbose = T,
              random.seed = 42)

EnhancedVolcano(toptable = GABA_anes_markers,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                lab = rownames(GABA_anes_markers),
                pCutoff = 0.05,
                FCcutoff = 0.5)

write.table(GABA_anes_markers,
            file = "GABA_anes_markers_MAST.tsv",
            quote = F, sep = "\t",
            row.names = T, col.names = T)
