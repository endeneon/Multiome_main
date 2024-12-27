# Chuxuan Li 05/05/2023
# Normalize -> integrate -> cluster -> label 030 data

# init ####
library(Seurat)
library(ggplot2)
library(pals)
library(stringr)
library(future)

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

load("/nvmefs/scARC_Duan_018/Duan_project_030_RNA/Analysis_part2_mapped_to_Hg38_only/GRCh38_mapped_after_QC_list.RData")

# normalize ####
transformed_lst <- vector(mode = "list", length(QCed_lst))

for (i in 1:length(QCed_lst)){
  print(paste0("now at ", i))
  # proceed with normalization
  obj <- PercentageFeatureSet(QCed_lst[[i]], pattern = c("^MT-"), col.name = "percent.mt")
  
  obj <- FindVariableFeatures(obj, nfeatures = 8000)
  transformed_lst[[i]] <- SCTransform(obj,
                                      vars.to.regress = "percent.mt",
                                      method = "glmGamPoi", 
                                      return.only.var.genes = F,
                                      variable.features.n = 8000,
                                      seed.use = 42,
                                      verbose = T)
}

# integration using reciprocal PCA ####
# find anchors
features <- SelectIntegrationFeatures(object.list = transformed_lst,
                                      nfeatures = 3000,
                                      fvf.nfeatures = 3000)
transformed_lst <- PrepSCTIntegration(transformed_lst, anchor.features = features)
transformed_lst <- lapply(X = transformed_lst, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(object.list = transformed_lst,
                                  anchor.features = features,
                                  reference = c(1, 2, 3),
                                  reduction = "rpca",
                                  normalization.method = "SCT",
                                  scale = F, 
                                  dims = 1:50)
#Found 9036 anchors
#all.genes <- transformed_lst[[1]]@assays$RNA@counts@Dimnames[[1]]

# integrate
integrated <- IntegrateData(anchorset = anchors,
                            verbose = T)
save(transformed_lst, file = "transformed_and_prepared_for_integration_list.RData")
save(integrated, file = "integrated_obj_nfeature_8000.RData")

# QC ####
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, assay = "RNA",
                                                   pattern = "^MT-")
DefaultAssay(integrated) <- "RNA"
Idents(integrated) <- "orig.ident"
VlnPlot(integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, fill.by = "feature", pt.size = 0)


# clustering ####
unique(integrated$orig.ident)
unique(integrated$time.ident)

DefaultAssay(integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, verbose = T, seed.use = 11)
integrated <- RunUMAP(integrated, reduction = "pca",
                      dims = 1:30, seed.use = 11)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5, random.seed = 11)

DimPlot(integrated, label = T, group.by = "seurat_clusters") +
  NoLegend() +
  ggtitle("Clustering 030 RNAseq data") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

lib_colors <- DiscretePalette(n = length(unique(Idents(integrated))), "alphabet")
p <- DimPlot(integrated, label = F, cols = lib_colors,
        group.by = "orig.ident") +
  ggtitle("by library") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
p$layers[[1]]$aes_params$alpha <- 0.1
p

p <- DimPlot(integrated, label = F, group.by = "time.ident") +
  ggtitle("by time point") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
p$layers[[1]]$aes_params$alpha <- 0.2
p

# add sequencing batch info and plot again
integrated$seq.batch.ident <- "030"
unique(integrated$lib.ident)
integrated$seq.batch.ident[integrated$lib.ident %in% c("2-0", "2-1", "2-6",
                                                       "8-0", "8-1", "8-6")] <- "018"
integrated$seq.batch.ident[integrated$lib.ident %in% c("22-0", "22-1", "22-6",
                                                       "36-0", "36-1", "36-6",
                                                       "39-0", "39-1","39-6",
                                                       "44-0", "44-1", "44-6", 
                                                       "49-0", "49-1", "49-6")] <- "022"
integrated$seq.batch.ident[integrated$lib.ident %in% c("05-0", "05-1", "05-6",
                                                       "33-0", "33-1", "33-6",
                                                       "35-0", "35-1", "35-6",
                                                       "51-0", "51-1", "51-6", 
                                                       "63-0", "63-1", "63-6")] <- "024"
integrated$seq.batch.ident[integrated$lib.ident %in% c("09-0", "09-1", "09-6", 
                                                       "13-0", "13-1", "13-6",
                                                       "21-0", "21-1", "21-6", 
                                                       "23-0", "23-1", "23-6",
                                                       "53-0", "53-1", "53-6",
                                                       "17-0", "17-1", "17-6",
                                                       "46-0", "46-1", "46-6")] <- "025"
integrated$seq.batch.ident[integrated$lib.ident %in% c("20087-0", "20087-1", "20087-6",
                                                       "20088-0", "20088-1", "20088-6",
                                                       "50040-0", "50040-1", "50040-6",
                                                       "60060-0", "60060-1", "60060-6",
                                                       "70179-0", "70179-1", "70179-6",
                                                       "40201-0", "40201-1", "40201-6")] <- "029"

sbatches <- c("018", "022", "024", "025", "029", "030")
p <- DimPlot(integrated_labeled, label = F, group.by = "seq.batch.ident") +
  ggtitle("by sequencing batch") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))
g <- ggplot_build(p)
colors <- unique(g$data[[1]]$colour)
color_order <- rep_len("transparent", length(sbatches))
for (i in 1:length(sbatches)) {
  color_order[i] <- colors[i]
  fname <- paste0("dimplot_by_seq_batch_", sbatches[i], ".png")
  png(fname, width = 600, height = 450)
  p <- DimPlot(integrated_labeled, label = F, group.by = "seq.batch.ident", cols = color_order) +
    ggtitle(paste0("sequencing batch ", sbatches[i])) +
    theme(text = element_text(size = 10), axis.text = element_text(size = 10))
  print(p)
  dev.off()
  color_order <- rep_len("transparent", length(sbatches))
}

# cell type identification ####
DefaultAssay(integrated) <- "SCT"
FeaturePlot(integrated, features = c("SOX2", "VIM"))
FeaturePlot(integrated, features = c("GAD1", "GAD2", "SLC17A6", "SLC17A7"), ncol = 2)
FeaturePlot(integrated, features = c("POU5F1", "NANOG"))
FeaturePlot(integrated, features = c("NEFM", "MAP2"))


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
  c("NEFM- glut", "NEFM+ glut", "NEFM- glut", "GABA", "NEFM- glut", "GABA",
    "GABA", "GABA", "SEMA3E+ GABA", "GABA", "NEFM+ glut",
    "SEMA3E+ GABA", "GABA", "unknown", "GABA", "SST+ GABA",
    "NEFM+ glut", "unknown", "BCL11B+ GABA", "unknown neuron", "unknown neuron")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated$seurat_clusters))

names(new.cluster.ids) <- levels(integrated)
integrated_labeled <- RenameIdents(integrated, new.cluster.ids)

DimPlot(integrated_labeled, reduction = "umap", label = TRUE, repel = F,
        pt.size = 0.3, cols = c("#E6D2B8", #nmglut,
                                "#CCAA7A", #npglut
                                "#B33E52", #GABA
                                "#CC7A88", #SEMA3E GABA
                                "#993F00", #unknown
                                "#E6B8BF", #SST GABA
                                "#4C005C", #BCL11B GABA
                                "#0075DC" #unknown neuron
        )) +
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

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

save(integrated_labeled, file = "integrated_labeled.RData")
