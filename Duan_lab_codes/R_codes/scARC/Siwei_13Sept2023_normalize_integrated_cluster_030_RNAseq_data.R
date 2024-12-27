# Chuxuan Li 05/05/2023
# Normalize -> integrate -> cluster -> label 030 data

# init ####
library(Seurat)
library(ggplot2)
library(pals)
library(stringr)
library(future)
library(dplyr)

plan("multisession", workers = 12)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

source("../../Duan_project_029_RNA/stacked_vln_copy.R")

QCed_lst <-
  readRDS("GRCh38_mapped_after_QC_list.RData")

# normalize ####
transformed_lst <-
  vector(mode = "list",
         length(QCed_lst))

for (i in 1:length(QCed_lst)){
  print(paste0("now at ", i))
  # proceed with normalization
  obj <-
    PercentageFeatureSet(QCed_lst[[i]],
                         pattern = c("^MT-"),
                         col.name = "percent.mt")

  obj <-
    FindVariableFeatures(obj,
                         nfeatures = 8000)
  transformed_lst[[i]] <-
    SCTransform(obj,
                vars.to.regress = "percent.mt",
                method = "glmGamPoi",
                return.only.var.genes = F,
                variable.features.n = 8000,
                seed.use = 42,
                verbose = T)
}

# integration using reciprocal PCA ####
# find anchors
features <-
  SelectIntegrationFeatures(object.list = transformed_lst,
                                      nfeatures = 3000,
                                      fvf.nfeatures = 3000)
transformed_lst <-
  PrepSCTIntegration(transformed_lst,
                     anchor.features = features)
transformed_lst <-
  lapply(X = transformed_lst,
         FUN = function(x) {
           x <-
             ScaleData(x,
                       features = features)
           x <-
             RunPCA(x,
                    features = features)
         })

anchors <-
  FindIntegrationAnchors(object.list = transformed_lst,
                         anchor.features = features,
                         reference = c(1, 2, 3),
                         reduction = "rpca",
                         normalization.method = "SCT",
                         scale = F,
                         dims = 1:50)
#Found 10806 anchors
#all.genes <- transformed_lst[[1]]@assays$RNA@counts@Dimnames[[1]]

# integrate
integrated <-
  IntegrateData(anchorset = anchors,
                verbose = T)
saveRDS(transformed_lst,
     file = "transformed_and_prepared_for_integration_list.RData")
saveRDS(integrated,
     file = "integrated_obj_nfeature_8000.RData")

# QC ####
integrated[["percent.mt"]] <-
  PercentageFeatureSet(integrated,
                       assay = "RNA",
                       pattern = "^MT-")
DefaultAssay(integrated) <- "RNA"
Idents(integrated) <- "orig.ident"
VlnPlot(integrated,
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt"),
        ncol = 3,
        fill.by = "feature",
        pt.size = 0)


# clustering ####
unique(integrated$orig.ident)
unique(integrated$time.ident)

DefaultAssay(integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated)
integrated <-
  RunPCA(integrated,
         verbose = T,
         seed.use = 42)
integrated <-
  RunUMAP(integrated,
          reduction = "pca",
          dims = 1:30,
          seed.use = 42)
integrated <-
  FindNeighbors(integrated,
                reduction = "pca",
                dims = 1:30)
integrated <-
  FindClusters(integrated,
               resolution = 0.5,
               random.seed = 42)

DimPlot(integrated,
        label = T,
        repel = T,
        group.by = "seurat_clusters") +
  # NoLegend() +
  ggtitle("Clustering 025 17/46 RNAseq data") +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))

lib_colors <-
  DiscretePalette(n = length(unique(Idents(integrated))),
                  "alphabet")
p <-
  DimPlot(integrated,
          label = F,
          cols = lib_colors,
          group.by = "orig.ident") +
  ggtitle("by library") +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))
p$layers[[1]]$aes_params$alpha <- 0.1
p

p <-
  DimPlot(integrated,
          label = F,
          group.by = "time.ident") +
  ggtitle("by time point") +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))
p$layers[[1]]$aes_params$alpha <- 0.2
p

# add sequencing batch info and plot again
integrated$seq.batch.ident <- "025_17_46"
unique(integrated$orig.ident)
integrated$seq.batch.ident[integrated$orig.ident %in% c("17-0",
                                                        "17-1",
                                                        "17-6")] <- "025-17"
integrated$seq.batch.ident[integrated$orig.ident %in% c("46-0",
                                                        "46-1",
                                                        "46-6")] <- "025-46"

sbatches <- c("025-17", "025-46")
p <-
  DimPlot(integrated,
          label = F,
          group.by = "seq.batch.ident") +
  ggtitle("by sequencing batch") +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))

g <- ggplot_build(p)
colors <- unique(g$data[[1]]$colour)
color_order <-
  rep_len("transparent",
          length.out = length(sbatches))

for (i in 1:length(sbatches)) {
  color_order[i] <- colors[i]
  fname <-
    paste0("dimplot_by_seq_batch_",
           sbatches[i],
           ".png")
  png(fname, width = 600, height = 600)
  p <- DimPlot(integrated,
               label = F,
               group.by = "seq.batch.ident",
               cols = color_order) +
    ggtitle(paste0("sequencing batch ",
                   sbatches[i])) +
    theme(text = element_text(size = 10),
          axis.text = element_text(size = 10))
  print(p)
  dev.off()
  color_order <-
    rep_len("transparent",
            length.out = length(sbatches))
}

# cell type identification ####
DefaultAssay(integrated) <- "SCT"
FeaturePlot(integrated,
            features = c("SOX2",
                         "VIM"))
FeaturePlot(integrated,
            features = c("GAD1",
                         "GAD2",
                         "SLC17A6",
                         "SLC17A7"), ncol = 2)
FeaturePlot(integrated,
            features = c("POU5F1", "NANOG"))
FeaturePlot(integrated,
            features = c("NEFM", "MAP2"))

FeaturePlot(integrated,
            features = c("SST", "BCL11B",
                         "GAD1",
                         "GAD2",
                         "SLC17A6",
                         "SLC17A7",
                         "NEFM", "MAP2",
                         "SOX2",
                         "VIM",
                         "POU5F1", "NANOG"),
            ncol = 4)


trimmed_markers <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
                     "EBF1", # striatal
                     "SEMA3E", # subcerebral
                     "BCL11B",  # cortical
                     "SST", # inhibitory
                     "SATB2",  "NEFM", # excitatory
                     "VIM", "SOX2",  #NPC
                     "SLC17A7", "SERTAD4", "FOXG1",  # forebrain
                     "POU3F2", "LHX2", # general cortex
                     "ADCYAP1", "CUX1", "CUX2",
                     "MAP2", "DCX")

StackedVlnPlot(obj = integrated, features = trimmed_markers) +
  coord_flip()

Idents(integrated) <- "seurat_clusters"

# assign cell types
new.cluster.ids <-
  c("NEFM+ glut", "NEFM- glut", "NEFM- glut", "NEFM+ glut", "GABA", "GABA", # [0-5]
    "NEFM- glut", "NEFM- glut", "unknown", "GABA", "NEFM- glut", # [6-10]
    "NEFM- glut", "GABA", "NEFM- glut", "GABA", "GABA", # [11-15]
    "unknown")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated$seurat_clusters))

names(new.cluster.ids) <- levels(integrated)
integrated_labeled <- RenameIdents(integrated, new.cluster.ids)

DimPlot(integrated_labeled,
        reduction = "umap",
        label = TRUE,
        repel = T,
        pt.size = 0.3,
        cols = c("#E6D2B844", #npglut,
                 "#4C005C44", #nmglut
                 "#B33E5244", #GABA
                 "#0075DC44", #SEMA3E GABA
                 "#0075DC", #unknown
                 "#E6B8BF", #SST GABA
                 "#4C005C", #BCL11B GABA
                 "#0075DC" #unknown neuron
        )) +
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10))

integrated_labeled$fine.cell.type <- integrated_labeled@active.ident
integrated_labeled$cell.type <- as.character(integrated_labeled$fine.cell.type)
integrated_labeled$cell.type[str_detect(integrated_labeled$fine.cell.type, "GABA")] <- "GABA"
integrated_labeled$cell.type[str_detect(integrated_labeled$fine.cell.type, "unknown")] <- "unidentified"
integrated_labeled$cell.type.forplot <- integrated_labeled$cell.type
integrated_labeled$cell.type[integrated_labeled$fine.cell.type == "NEFM- glut"] <- "nmglut"
integrated_labeled$cell.type[integrated_labeled$fine.cell.type == "NEFM+ glut"] <- "npglut"
unique(integrated_labeled$cell.type)
unique(integrated_labeled$cell.type.forplot)
integrated_labeled$cell.type.counts <- integrated_labeled$cell.type.forplot
types <- unique(integrated_labeled$cell.type.forplot)
for (i in 1:length(types)) {
  count <- sum(integrated_labeled$cell.type.forplot == types[i])
  print(count)
  integrated_labeled$cell.type.counts[integrated_labeled$cell.type.counts == types[i]] <-
    paste0(types[i], "\n", count)
}
unique(integrated_labeled$cell.type.counts)
DimPlot(integrated_labeled, reduction = "umap", group.by = "cell.type.counts",
        label = TRUE, repel = F, pt.size = 0.3,
        cols = c("#B33E52", "#E6D2B8", "#CCAA7A", "#54990F")) +
  ggtitle("labeled by cell type and number of cells in each type") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))


integrated_labeled$cell.type <- Idents(integrated_labeled)
save(integrated_labeled, file = "integrated_labeled_025_17_46.RData")


unique(integrated_labeled$orig.ident)
unique(integrated_labeled$cell.line.ident)
unique(integrated_labeled$seurat_clusters)
unique(integrated_labeled$cell.type)

dir.create("all_barcodes")

for (time_ident in unique(integrated_labeled$time.ident)) {
  for (cell_line_ident in unique(integrated_labeled$cell.line.ident)) {
    for (cell_type in unique(integrated_labeled$cell.type)) {
      current_output_set <-
        subset(integrated_labeled,
               subset = (time.ident == time_ident &
                           cell.line.ident == cell_line_ident &
                           cell.type == cell_type))
      writeout_barcodes <-
        colnames(current_output_set)
      writeout_barcodes <-
        str_split(string = writeout_barcodes,
                  pattern = "_",
                  simplify = T)[, 1]
      print(length(writeout_barcodes))

      file_name_string <-
        str_c(time_ident,
              cell_line_ident,
              cell_type,
              sep = "_") %>%
        str_replace_all(c('\\- glut' = '_neg_glut',
                          '\\+ glut' = '_pos_glut'))
      file_name_string <-
        str_c(file_name_string,
              "_barcodes.txt")

      print(file_name_string)

      if (!str_detect(string = file_name_string,
                      pattern = "unknown")) {
        write.table(writeout_barcodes,
                    file = paste0("all_barcodes/",
                                  file_name_string),
                    quote = F, sep = "\t",
                    row.names = F, col.names = F)
      }
    }
  }
}
