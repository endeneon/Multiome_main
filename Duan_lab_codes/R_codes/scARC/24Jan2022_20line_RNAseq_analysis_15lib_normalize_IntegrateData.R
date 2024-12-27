# Chuxuan Li 12/23/2021
# RNAseq data analysis by SCT-normalizing 6 libraries individually first, then
# use IntegrateData() in Seurat to integrate 6 libraries


# init ####
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(monocle3)

library(stringr)
library(ggplot2)
library(readr)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)

library(future)

# set threads and parallelization
plan("multisession", workers = 1)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 107374182400)


# load data ####
setwd("/data/FASTQ/Duan_Project_022_Reseq")
# note this h5 file contains both atac-seq and gex information
# Read files, separate the ATACseq data from the .h5 matrix
h5list <- list.files(path = ".", 
                     pattern = "filtered_feature_bc_matrix.h5",
                     recursive = T,
                     include.dirs = T)
h5list <- h5list[2:length(h5list)]

objlist <- vector(mode = "list", length = length(h5list))
transformed_lst <- vector(mode = "list", length(h5list))

# normalization
for (i in 1:length(objlist)){
  h5file <- Read10X_h5(filename = h5list[i])
  print(str_extract(string = h5list[i], 
                    pattern = "[0-9][0-9]_[0-6]"))
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(string = h5list[i], 
                                                  pattern = "[0-9][0-9]_[0-6]"))
  obj$lib.ident <- str_extract(string = h5list[i], 
                               pattern = "[0-9][0-9]_[0-6]")
  print(paste0(i, " number of genes: ", nrow(obj), 
               ", number of cells: ", ncol(obj)))
  obj <- 
    PercentageFeatureSet(obj,
                         pattern = c("^MT-"),
                         col.name = "percent.mt")
  objlist[[i]] <- FindVariableFeatures(obj, nfeatures = 3000)
  transformed_lst[[i]] <- SCTransform(objlist[[i]], 
                                      vars.to.regress = "percent.mt",
                                      method = "glmGamPoi",
                                      variable.features.n = 8000,
                                      seed.use = 115,
                                      verbose = T)
}

setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat")
# integration using reciprocal PCA ####
# find anchors
features <- SelectIntegrationFeatures(object.list = transformed_lst, 
                                      nfeatures = 3000, 
                                      fvf.nfeatures = 3000)

transformed_lst <- lapply(X = transformed_lst, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})
# save.image(file = "24Jan2022_20lines_before_integration.RData")

anchors <- FindIntegrationAnchors(object.list = transformed_lst, 
                                  reference = c(1, 2, 3), 
                                  reduction = "rpca",
                                  dims = 1:50)

# integrate
integrated <- IntegrateData(anchorset = anchors, 
                            verbose = T)
# save("integrated", file = "obj_after_integration.RData")


# downstream analysis
# integrated$group.ident <- NA
# integrated$group.ident[integrated$lib.ident %in% c(1, 2, 3)] <- "lib_22"
# integrated$group.ident[integrated$lib.ident %in% c(4, 5, 6)] <- "lib_36"
# integrated$group.ident[integrated$lib.ident %in% c(4, 5, 6)] <- "lib_36"
unique(integrated$orig.ident)
unique(integrated$lib.ident)
integrated$time.ident <- NA
zerohr_origidents <- vector(mode = "list", length = length(h5list) / 3)
onehr_origidents <- vector(mode = "list", length = length(h5list) / 3)
sixhr_origidents <- vector(mode = "list", length = length(h5list) / 3)

for (i in 1:length(h5list)){
  zerohr_origidents[i] <- str_extract(string = h5list[i], 
                                      pattern = "[0-9][0-9]_0")
  onehr_origidents[i] <- str_extract(string = h5list[i], 
                                      pattern = "[0-9][0-9]_1")
  sixhr_origidents[i] <- str_extract(string = h5list[i], 
                                      pattern = "[0-9][0-9]_6")
}
integrated$time.ident[integrated$orig.ident %in% zerohr_origidents] <- "0hr"
integrated$time.ident[integrated$lib.ident %in% onehr_origidents] <- "1hr"
integrated$time.ident[integrated$lib.ident %in% sixhr_origidents] <- "6hr"
unique(integrated$time.ident)


# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, 
                     verbose = T, 
                     seed.use = 11)
integrated <- RunUMAP(integrated, 
                      reduction = "pca", 
                      dims = 1:30,
                      seed.use = 11)
integrated <- FindNeighbors(integrated, 
                            reduction = "pca", 
                            dims = 1:30)
integrated <- FindClusters(integrated, 
                           resolution = 0.5,
                           random.seed = 11)

DimPlot(integrated, 
        label = T) +
  NoLegend() +
  ggtitle("clustering of independently normalized, \n Seurat-integrated data")


# cell type labeling ####
# check time point distribution 
integrated_dimplot <- DimPlot(integrated,
                              #cols = c("steelblue3", "mediumpurple4", "indianred3"),
                              cols = c("cyan", "magenta", "yellow"),
                              # cols = "Set2",
                              # cols = g_2_6_gex@meta.data$cell.line.ident,
                              group.by = "time.ident",
                              label = F) +
  ggtitle("independently normalized RNAseq data \n integrated by Seurat \n colored by time")
integrated_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
integrated_dimplot

# check library distribution
by_lib_dimplot <- DimPlot(integrated,
                          #cols = c("steelblue3", "mediumpurple4", "indianred3"),
                          
                          #cols = c("transparent", "transparent", "transparent", "transparent", "transparent", "steelblue"),
                          # cols = g_2_6_gex@meta.data$cell.line.ident,
                          group.by = "lib.ident",
                          label = F) +
  ggtitle("independently normalized RNAseq data \n integrated by Seurat \n colored by library")
by_lib_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
by_lib_dimplot

FeaturePlot(integrated, 
            features = c("SOX2", "GAD1", "SLC17A6", "VIM"))
DefaultAssay(integrated) <- "SCT"
FeaturePlot(integrated, 
            features = c("POU5F1", "NANOG"))
FeaturePlot(integrated, 
            features = c("MAP2"))
FeaturePlot(integrated, 
            features = c("NEFM", "CUX2"))

trimmed_markers <- c("GAD1", "GAD2", "SLC17A6",
                     "EBF1", # striatal
                     "SEMA3E", # subcerebral
                     "BCL11B",  # cortical
                     "SST", # inhibitory
                     "SATB2",  "NEFM", # excitatory
                     "VIM", "SOX2", "FOXG1"  #NPC
)

other_markers <- c("SLC17A7", "SERTAD4", "FOXG1", # forebrain
                   "POU3F2", "LHX2" # general cortex
)

DefaultAssay(integrated) <- "integrated"
StackedVlnPlot(obj = integrated, features = trimmed_markers) +
  coord_flip()
StackedVlnPlot(obj = integrated, features = other_markers) +
  coord_flip()

early_late_markers <- c("EGR1", "FOS", "FOSB", "NPAS4", # early response
                        "VGF", "BDNF", "NRN1", "PNOC" # late response
)
StackedVlnPlot(obj = integrated, features = early_late_markers) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

FeaturePlot(integrated, features = c("SEMA3E"))
FeaturePlot(integrated, features = c("FOXG1"))
FeaturePlot(integrated, features = c("GFAP"))

# assign cell types
new.cluster.ids <- 
  c("GABA", "NEFM_pos_glut", "NPC", "GABA", "NEFM_neg_glut",
    "GABA", "GABA", "GABA", "GABA", "unknown",
    "NEFM_pos_glut", "GABA", "unknown", "unknown", "unknown",
    "unknown")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated$seurat_clusters))


names(new.cluster.ids) <- levels(integrated)

RNAseq_integrated_labeled <- RenameIdents(integrated, new.cluster.ids)
DimPlot(RNAseq_integrated_labeled, 
        reduction = "umap", 
        label = TRUE, 
        repel = T, 
        pt.size = 0.3) +
  ggtitle("independently normalized RNAseq data \n integrated by Seurat, labeled by cell type")

sum(RNAseq_integrated_labeled@active.ident == "NPC") #11838
sum(RNAseq_integrated_labeled@active.ident == "NEFM_neg_glut") #7482
sum(RNAseq_integrated_labeled@active.ident == "NEFM_pos_glut") #17602
sum(RNAseq_integrated_labeled@active.ident == "GABA") #47775
sum(RNAseq_integrated_labeled@active.ident == "unknown") #10693

RNAseq_integrated_labeled$broad.cell.type <- NA
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled@active.ident == "NEFM_pos_glut"] <- "NEFM_pos_glut"
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled@active.ident == "NEFM_neg_glut"] <- "NEFM_neg_glut"
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled@active.ident == "NPC"] <- "NPC"
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled@active.ident == "GABA"] <- "GABA"
RNAseq_integrated_labeled$broad.cell.type[RNAseq_integrated_labeled@active.ident == "unknown"] <- "unknown"
unique(RNAseq_integrated_labeled$broad.cell.type)

save("RNAseq_integrated_labeled", file = "nfeature5000_labeled_obj.RData")


# differential expression ####
unique(Idents(RNAseq_integrated_labeled))
Idents(RNAseq_integrated_labeled) <- "lib.ident"

GABA <- subset(RNAseq_integrated_labeled, broad.cell.type == "GABA")
glut <- subset(RNAseq_integrated_labeled, broad.cell.type %in% c("NEFM_pos_glut", 
                                                                 "NEFM_neg_glut"))

j <- 0
# calc hour 1 vs 0
for (i in sort(unique(GABA$lib.ident))) {
  print(i)
  j = j + 1
  if (j == 1) {
    print("first lib")
    df_1vs0 <- FindMarkers(object = GABA,
                           ident.1 = "1hr",
                           group.by = 'time.ident',
                           features = c("FOS", "VGF"),
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_1vs0$gene_symbol <- rownames(df_1vs0)
    df_1vs0$lib <- i
  } else {
    print("others")
    df_to_append <- FindMarkers(object = GABA,
                                ident.1 = "1hr",
                                group.by = 'time.ident',
                                features = c("FOS", "VGF"),
                                subset.ident = i,
                                ident.2 = "0hr",
                                min.pct = 0.0,
                                logfc.threshold = 0.0)
    df_to_append$gene_symbol <- rownames(df_to_append)
    df_to_append$lib <- i
    df_1vs0 <- rbind(df_1vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
  }
}

j <- 0
# calculate 6hr vs 0hr
for (i in sort(unique(GABA$lib.ident))) {
  print(i)
  j = j + 1
  if (j == 1) {
    print("first cell type")
    df_6vs0 <- FindMarkers(object = GABA,
                           ident.1 = "6hr",
                           group.by = 'time.ident',
                           features = c("FOS", "VGF"),
                           subset.ident = i,
                           ident.2 = "0hr",
                           min.pct = 0.0,
                           logfc.threshold = 0.0)
    df_6vs0$gene_symbol <- rownames(df_6vs0)
    df_6vs0$lib <- i
  } else {
    print("others")
    df_to_append <- FindMarkers(object = GABA,
                                ident.1 = "6hr",
                                group.by = 'time.ident',
                                features = c("FOS", "VGF"),
                                subset.ident = i,
                                ident.2 = "0hr",
                                min.pct = 0.0,
                                logfc.threshold = 0.0)
    df_to_append$gene_symbol <- rownames(df_to_append)
    df_to_append$lib <- i
    df_6vs0 <- rbind(df_6vs0,
                     df_to_append, 
                     make.row.names = F,
                     stringsAsFactors = F)
  }
}

colnames(df_6vs0) <- colnames(df_1vs0)
# make a bubble plot
df_to_plot <- rbind(df_1vs0,
                    df_6vs0)
df_to_plot$source <- c(rep_len("1v0hr", length.out = 10),
                       rep_len("6v0hr", length.out = 10))


ggplot(df_to_plot,
       aes(x = as.factor(source),
           y = as.factor(lib),
           size = pct.1 * 100,
           fill = avg_log2FC)) +
  geom_point(shape = 21) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  facet_grid(cols = vars(gene_symbol)) + 
  theme_bw() 

  
sum(RNAseq_integrated_labeled@assays$integrated@data@Dimnames[[1]] == "NPAS4")


unique(RNAseq_integrated_labeled$lib.ident)
RNAseq_integrated_labeled$lib.ident <- str_sub(RNAseq_integrated_labeled$lib.ident,
                                               end = -3L)
DimPlot(RNAseq_integrated_labeled, 
        group.by = "lib.ident",
        cols = c("transparent", "transparent", "transparent", "transparent", "mediumpurple2")) +
  ggtitle("colored by one library")



# Gene set enrichment ####

# first calculate time de markers
Idents(RNAseq_integrated_labeled) <- "time.ident"
df_1vs0 <- FindMarkers(object = RNAseq_integrated_labeled,
                      ident.1 = "1hr",
                      ident.2 = "0hr",
                      min.pct = 0.0,
                      logfc.threshold = 0.0)
df_1vs0$gene_symbol <- rownames(df_1vs0)

df_6vs0 <- FindMarkers(object = RNAseq_integrated_labeled,
                       ident.1 = "6hr",
                       ident.2 = "0hr",
                       min.pct = 0.0,
                       logfc.threshold = 0.0)
df_6vs0$gene_symbol <- rownames(df_6vs0)

df_1vs0_filtered <- df_1vs0[df_1vs0$p_val_adj < 0.05, ]
df_6vs0_filtered <- df_6vs0[df_6vs0$p_val_adj < 0.05, ]

# output for GO term analysis
write.table(df_1vs0_filtered, file = "de_1v0.csv", quote = F, sep = ",", col.names = T, row.names = F)
write.table(df_6vs0_filtered, file = "de_6v0.csv", quote = F, sep = ",", col.names = T, row.names = F)
write.table(df_1vs0_filtered$gene_symbol, file = "de_1v0_genes.csv", quote = F, 
            sep = ",", col.names = F, row.names = F)
write.table(df_6vs0_filtered$gene_symbol, file = "de_6v0_genes.csv", quote = F, 
            sep = ",", col.names = F, row.names = F)

# read GO term results
analysis_1v0 <- 
  read.delim("analysis_1v0.txt")
analysis_6v0 <- 
  read.delim("analysis_6v0.txt")
col_names <- c("name", "REF", "#", "expected", "over_under", 
               "fold_enrichment", "raw_p_value", "FDR")
colnames(analysis_1v0) <- col_names
colnames(analysis_6v0) <- col_names
up_1v0 <- analysis_1v0[analysis_1v0$over_under == "+", ]
up_6v0 <- analysis_6v0[analysis_6v0$over_under == "+", ]
do_1v0 <- analysis_1v0[analysis_1v0$over_under == "-", ]
do_6v0 <- analysis_6v0[analysis_6v0$over_under == "-", ]
up <- rbind(up_1v0[978:997, ], up_6v0[1184:1203, ])
up$time <- c(rep_len("1v0", 20), rep_len("6v0", 20))
up$name <- gsub(pattern = '\\(GO:[0-9]+\\)', replacement = "", x = up$name)

do <- rbind(do_1v0, do_6v0)
do <- do[!is.na(do$name), ]
do$time <- c(rep_len("1v0", 18), rep_len("6v0", 19))
do$name <- gsub(pattern = "\\(GO:[0-9]+\\)", replacement = "", x = do$name)

# plot with dotplot
ggplot(up,
       aes(x = as.factor(time),
           y = as.factor(name),
           size = (0 - log10(raw_p_value)),
           fill = fold_enrichment)) +
  geom_point(shape = 21) +
  labs(title = "Gene set enrichment (upregulated)\n22 out of 1219 GO terms", x = "", y = "", 
       fill = "fold\nenrichment \n", size = "p value") +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))

save("df_1vs0_filtered", "df_6vs0_filtered", file = "filtered_de_time_markers.RData")
ggplot(do,
       aes(x = as.factor(time),
           y = as.factor(name),
           size = (0 - log10(raw_p_value)),
           fill = fold_enrichment)) +
  geom_point(shape = 21) +
  labs(title = "Gene set\nenrichment\n(downregulated)", x = "", y = "", 
       fill = "fold\nenrichment \n", size = "p value") +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))


# pseudotime trajectory ####
expmatrix <- RNAseq_integrated_labeled@assays$RNA@counts
cellmeta <- RNAseq_integrated_labeled@meta.data
cellmeta$cell <- rownames(cellmeta)
genemeta <- data.frame(gene_short_name = RNAseq_integrated_labeled@assays$RNA@counts@Dimnames[[1]],
                      numreads = rowSums(RNAseq_integrated_labeled@assays$RNA@counts))
sum(rownames(cellmeta) != colnames(expmatrix))
sum(rownames(genemeta) != rownames(expmatrix))
# create cds
cds <- new_cell_data_set(expression_data = expmatrix,
                         cell_metadata = cellmeta,
                         gene_metadata = genemeta)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, 
                 alignment_group = "lib.ident")

cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster = FALSE,  color_cells_by = "broad.cell.type")
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "broad.cell.type",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = T)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

