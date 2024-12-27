# Chuxuan Li 2/9/2022
# Using demultiplexed barcodes for each cell line, identify cells in rpca-integrated
#data by cell line, then remove cells not matched to any cell line, and proceed
#with analysis


# init ####
library(Seurat)
library(glmGamPoi)
library(sctransform)

library(ggplot2)
library(RColorBrewer)
library(viridis)
library(graphics)
library(ggrepel)

library(dplyr)
library(readr)
library(stringr)
library(readxl)

library(future)
# set threads and parallelization
plan("multisession", workers = 2)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

# load data ####
setwd("/data/FASTQ/Duan_Project_022_Reseq")
# note this h5 file contains both atac-seq and gex information
# Read files, separate the ATACseq data from the .h5 matrix
h5list <- list.files(path = ".",
                     pattern = "filtered_feature_bc_matrix.h5",
                     recursive = T,
                     include.dirs = T)
h5list <- h5list[!h5list %in% grep("hybrid", h5list, value = T)]
h5list <- h5list[2:length(h5list)]
objlist <- vector(mode = "list", length = length(h5list))
transformed_lst <- vector(mode = "list", length(h5list))

for (i in 1:length(objlist)){
  
  h5file <- Read10X_h5(filename = h5list[i])
  print(str_extract(string = h5list[i],
                    pattern = "[0-9][0-9]_[0-6]"))
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(string = h5list[i],
                                                  pattern = "[0-9][0-9]_[0-6]"))
  # assign library identity
  obj$lib.ident <- str_extract(string = h5list[i],
                               pattern = "[0-9][0-9]_[0-6]")
  print(obj$lib.ident)
  # check number of genes and cells
  print(paste0(i, " number of genes: ", nrow(obj),
               ", number of cells: ", ncol(obj)))
  objlist[[i]] <- obj
}

# read in demultiplexed barcodes


cdlist <- data.frame(lib22 = c("CD_22", "CD_23", "CD_37", "CD_40"),
                     lib36 = c("CD_36", "CD_38", "CD_43", "CD_45"),
                     lib39 = c("CD_39", "CD_42", "CD_47", "CD_48"),
                     lib44 = c("CD_31", "CD_44", "CD_46", "CD_57"),
                     lib49 = c("CD_32","CD_49",  "CD_58", "CD_61"))

id <- c("22_0", "22_1", "22_6",
"36_0", "36_1", "36_6", 
"39_0", "39_1", "39_6", 
"44_0", "44_1", "44_6",
"49_0", "49_1","49_6")

for (i in 1:5){
  # read in barcodes for each lib
  for(j in 1:3){
    print(id[i * 3 - 3 + j])
    df <-
      read.delim(paste0("./demux_barcode_output/output_barcodes/libraries_",
                        id[i * 3 - 3 + j],
                        "_SNG_by_cell_line.tsv"),
                 quote="")
    obj <- objlist[[i * 3 - 3 + j]]
    print(unique(obj$orig.ident))
    
    # assign cell line looping over 4 cell lines per lib
    obj$cell.line.ident <- "unknown"
    for (k in 1:4){
      # assign cell line in each lib
      ident <- cdlist[k, i]
      print(ident)
      subdf <- df$BARCODE[df$SNG.1ST == ident]
      obj$cell.line.ident[obj@assays$RNA@counts@Dimnames[[2]] %in% subdf] <- ident
    }
    print(unique(obj$cell.line.ident))
    objlist[[i * 3 - 3 + j]] <- obj
  }
  
}


# check counts of unknowns
unknown_counts <- rep_len(0, length(objlist))
for(i in 1:length(objlist)){
  unknown_counts[i] <- (sum(objlist[[i]]$cell.line.ident == "unknown"))
}
sum(unknown_counts)

sum(obj$cell.line.ident == "unknown")
subobjlist <- vector(mode = "list", length = length(objlist))
for (i in 1:length(objlist)){
  # proceed with normalization
  obj <- objlist[[i]]
  subobjlist[[i]] <- subset(obj, cell.line.ident != "unknown")
  print(unique(subobjlist[[i]]$cell.line.ident))
}

# normalization
for (i in 1:length(objlist)){
  # proceed with normalization
  subobjlist[[i]] <-
    PercentageFeatureSet(subobjlist[[i]],
                         pattern = c("^MT-"),
                         col.name = "percent.mt")
  
  subobjlist[[i]] <- FindVariableFeatures(subobjlist[[i]], nfeatures = 3000)
  transformed_lst[[i]] <- SCTransform(subobjlist[[i]],
                                      vars.to.regress = "percent.mt",
                                      method = "glmGamPoi",
                                      variable.features.n = 8000,
                                      seed.use = 115,
                                      verbose = T)
}


setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat")
save(transformed_lst, file = "20line_demuxed_removed_nonhuman_transformed_lst.RData")

# try integration without rPCA ####
features <- SelectIntegrationFeatures(object.list = transformed_lst,
                                      nfeatures = 8000, fvf.nfeatures = 8000)
intg_lst <- PrepSCTIntegration(object.list = transformed_lst,
                               anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = intg_lst,
                                  normalization.method = "SCT",
                                  anchor.features = features)
demuxed_integrated <- IntegrateData(anchorset = anchors,
                                    normalization.method = "SCT")

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
                                  anchor.features = features,
                                  reference = c(1, 2, 3),
                                  reduction = "rpca",
                                  dims = 1:50)
#all.genes <- objlist[[1]]@assays$RNA@counts@Dimnames[[1]]

# integrate
integrated <- IntegrateData(anchorset = anchors,
                            features.to.integrate = features,
                            verbose = T)
save("integrated", file = "demuxed_obj_after_integration_all_genes.RData")


# clustering ####
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

DefaultAssay(integrated) <- "integrated"
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
        label = T,
        group.by = "seurat_clusters") +
  NoLegend() +
  ggtitle("clustering of 20-cell line RNAseq data\n(removed cells not matched to any cell line)") +
  theme(text = element_text(size = 12))



# remove clusters 5 and 13 which are probably rat cells, then rerun graphs
integrated$rat.ident <- "human"
integrated$rat.ident[integrated$seurat_clusters %in% c(13, 17, 20)] <- "rat"
sum(integrated$seurat_clusters %in% c(13, 17, 20)) #5660
sum(integrated$rat.ident ==  "rat") #5660
clean_integrated <- subset(integrated, subset = rat.ident != "rat")
length(unique(clean_integrated$seurat_clusters))


clean_integrated <- RunPCA(clean_integrated,
                           verbose = T,
                           seed.use = 11)
clean_integrated <- RunUMAP(clean_integrated,
                            reduction = "pca",
                            dims = 1:30,
                            seed.use = 11)
clean_integrated <- FindNeighbors(clean_integrated,
                                  reduction = "pca",
                                  dims = 1:30)
clean_integrated <- FindClusters(clean_integrated,
                                 resolution = 1.3,
                                 random.seed = 11)

DimPlot(clean_integrated,
        label = T,
        group.by = "seurat_clusters") +
  NoLegend() +
  ggtitle("clustering of 20-cell line RNAseq data\n(removed cells not matched to any cell line)") +
  theme(text = element_text(size = 12))


Idents(clean_integrated) <- "seurat_clusters"
# check time point distribution
integrated_dimplot <- DimPlot(clean_integrated,
                              #cols = c("steelblue3", "mediumpurple4", "indianred3"),
                              cols = c("cyan", "magenta", "yellow"),
                              # cols = "Set2",
                              # cols = g_2_6_gex@meta.data$cell.line.ident,
                              group.by = "time.ident",
                              label = F) +
  ggtitle("clustering of 20-cell line RNAseq data\n(removed cells not matched to any cell line)\n colored by time")
integrated_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
integrated_dimplot

DimPlot(clean_integrated,
        #cols = c("steelblue3", "mediumpurple4", "indianred3"),
        cols = c("cyan", "transparent", "transparent"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("0hr")
DimPlot(clean_integrated,
        #cols = c("steelblue3", "mediumpurple4", "indianred3"),
        cols = c("transparent", "magenta", "transparent"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("1hr")
DimPlot(clean_integrated,
        #cols = c("steelblue3", "mediumpurple4", "indianred3"),
        cols = c("transparent", "transparent", "gold"),
        # cols = "Set2",
        # cols = g_2_6_gex@meta.data$cell.line.ident,
        group.by = "time.ident",
        label = F) +
  ggtitle("6hr")

# check library distribution
by_lib_dimplot <- DimPlot(clean_integrated,
                          #cols = c("steelblue3", "mediumpurple4", "indianred3"),
                          
                          #cols = c("transparent", "transparent", "transparent", "transparent", "transparent", "steelblue"),
                          # cols = g_2_6_gex@meta.data$cell.line.ident,
                          group.by = "lib.ident",
                          label = F) +
  ggtitle("clustering of 20-cell line RNAseq data\n(removed cells not matched to any cell line)\ncolored by library")
by_lib_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
by_lib_dimplot

by_lib_dimplot <- DimPlot(clean_integrated,
                          #cols = c("steelblue3", "mediumpurple4", "indianred3"),
                          
                          #cols = c("transparent", "transparent", "transparent", "transparent", "transparent", "steelblue"),
                          # cols = g_2_6_gex@meta.data$cell.line.ident,
                          group.by = "cell.line.ident",
                          label = F) +
  ggtitle("clustering of 20-cell line RNAseq data\n(removed cells not matched to any cell line)\ncolored by cell line")
by_lib_dimplot[[1]]$layers[[1]]$aes_params$alpha = .2
by_lib_dimplot



# cell type identification ####
DefaultAssay(clean_integrated) <- "SCT"
FeaturePlot(clean_integrated,
            features = c("SOX2", "VIM"))
FeaturePlot(clean_integrated,
            features = c("GAD1", "SLC17A6"))
FeaturePlot(clean_integrated,
            features = c("POU5F1", "NANOG"))
FeaturePlot(clean_integrated,
            features = c("MAP2"))
FeaturePlot(clean_integrated,
            features = c("NEFM", "CUX2"))


trimmed_markers <- c("GAD1", "GAD2", "SLC17A6",
                     "EBF1", # striatal
                     "SEMA3E", # subcerebral
                     "BCL11B",  # cortical
                     "SST", # inhibitory
                     "SATB2",  "NEFM", # excitatory
                     "VIM", "SOX2",  #NPC
                     "SLC17A7", "SERTAD4", "FOXG1",  # forebrain
                     "POU3F2", "LHX2", # general cortex
                     "ADCYAP1", "CUX1", "CUX2", 
                     "MAP2", "DCX"
)


DefaultAssay(integrated_labeled) <- "SCT"
Idents(integrated_labeled) <- "seurat_clusters"
StackedVlnPlot(obj = integrated_labeled, features = trimmed_markers) +
  coord_flip()


early_late_markers <- c("EGR1", "FOS", "FOSB", "NPAS4", # early response
                        "VGF", "BDNF", "NRN1", "PNOC" # late response
)
StackedVlnPlot(obj = integrated, features = early_late_markers) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

FeaturePlot(clean_integrated, features = c("SEMA3E"))
FeaturePlot(clean_integrated, features = c("FOXG1"))

# assign cell types
new.cluster.ids <-
  c("GABA", "NEFM_pos_glut", "NEFM_neg_glut", "GABA", "GABA",
    "GABA", "GABA", "SEMA3E_pos_GABA", "GABA", "GABA",
    "SEMA3E_pos_GABA", "GABA", "NEFM_pos_glut", "NEFM_neg_glut", "NEFM_neg_glut",
    "GABA", "NEFM_pos_glut", "unknown", "unknown", "NEFM_neg_glut",
    "unknown")
#unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(clean_integrated$seurat_clusters))


names(new.cluster.ids) <- levels(clean_integrated)
integrated_labeled <- RenameIdents(clean_integrated, new.cluster.ids)
DimPlot(integrated_labeled,
        reduction = "umap",
        label = TRUE,
        repel = F,
        pt.size = 0.3) +
  ggtitle("clustering of 20-cell line RNAseq data\nlabeled by cell type") +
  theme(text = element_text(size = 11))



integrated_labeled$cell.type <- NA
integrated_labeled$cell.type[integrated_labeled@active.ident == "NEFM_pos_glut"] <- "NEFM_pos_glut"
integrated_labeled$cell.type[integrated_labeled@active.ident == "NEFM_neg_glut"] <- "NEFM_neg_glut"
integrated_labeled$cell.type[integrated_labeled@active.ident == "GABA"] <- "GABA"
integrated_labeled$cell.type[integrated_labeled@active.ident == "unknown"] <- "unknown"
integrated_labeled$cell.type[integrated_labeled@active.ident == "SEMA3E_pos_GABA"] <- "SEMA3E_pos_GABA"
unique(integrated_labeled$cell.type)

# count number in each cell type
integrated_labeled <- subset(integrated_labeled, seurat_clusters != "14")
sum(integrated_labeled@active.ident == "NEFM_neg_glut") #8210
sum(integrated_labeled@active.ident == "NEFM_pos_glut") #11947
sum(integrated_labeled@active.ident == "GABA") #36885
sum(integrated_labeled@active.ident == "SEMA3E_pos_GABA") #6989
sum(integrated_labeled@active.ident == "unknown") #3354

integrated_labeled$cell.type.counts <- 0
integrated_labeled$cell.type.counts[integrated_labeled$cell.type == "NEFM_pos_glut"] <- "NEFM_pos_glut\n11947"
integrated_labeled$cell.type.counts[integrated_labeled$cell.type == "NEFM_neg_glut"] <- "NEFM_neg_glut\n8210"
integrated_labeled$cell.type.counts[integrated_labeled$cell.type == "GABA"] <- "GABA\n36885"
integrated_labeled$cell.type.counts[integrated_labeled$cell.type == "unknown"] <- "unknown\n3354"
integrated_labeled$cell.type.counts[integrated_labeled$cell.type == "SEMA3E_pos_GABA"] <- "SEMA3E_pos_GABA\n6989"
unique(integrated_labeled$cell.type.counts)

DimPlot(integrated_labeled,
        reduction = "umap",
        label = TRUE,
        repel = F,
        pt.size = 0.3,
        group.by = "cell.type.counts") +
  NoLegend() +
  ggtitle("clustering of 20-cell line RNAseq data\nlabeled by cell type") +
  theme(text = element_text(size = 11))

# count number of cells in each cluster
integrated_labeled$cluster.counts <- 0
for (i in unique(integrated_labeled$seurat_clusters)){
  lab <- sum(integrated_labeled$seurat_clusters == i)
  print(lab)
  integrated_labeled$cluster.counts[integrated_labeled$seurat_clusters == i] <- lab
  
}
unique(integrated_labeled$cluster.counts)

DimPlot(integrated_labeled,
        reduction = "umap",
        label = TRUE,
        repel = F,
        pt.size = 0.3,
        group.by = "cluster.counts"
        #order = c(7957,7319,4780,4745,4558,3955,3800,3770,3677,3332,3219,2988,
        #2902, 2381, 1873, 1726, 1675, 1321, 1049, 358)
) +
  scale_color_discrete(label = rev(0:19)) +
  ggtitle("clustering of 20-cell line RNAseq data\nwith cell counts by cluster") +
  theme(text = element_text(size = 11))


# count number of cells in each lib
lib.idents <- unique(integrated_labeled$lib.ident)
lib_count <- rep_len(0, length.out = length(lib.idents))
for (i in 1:length(lib.idents)){
  print(i)
  lib_count[i] <- sum(integrated_labeled$lib.ident == lib.idents[i])
}
names(lib_count) <- lib.idents
lib_count_table <- rbind(lib22 = lib_count[1:3],
                         lib36 = lib_count[4:6],
                         lib39 = lib_count[7:9],
                         lib44 = lib_count[10:12],
                         lib49 = lib_count[13:15])
colnames(lib_count_table) <- c("0hr", "1hr", "6hr")


# count number of cells in each cell line
line.idents <- unique(integrated_labeled$cell.line.ident)
line_count <- rep_len(0, length.out = length(line.idents))
for (i in 1:length(line.idents)){
  print(i)
  line_count[i] <- sum(integrated_labeled$cell.line.ident == line.idents[i])
}
names(line_count) <- line.idents
lib_count_table <- rbind(lib22 = lib_count[1:3],
                         lib36 = lib_count[4:6],
                         lib39 = lib_count[7:9],
                         lib44 = lib_count[10:12],
                         lib49 = lib_count[13:15])

save("clean_integrated", file = "demux_integrated_cleaned_obj.RData")
save("integrated_labeled", file = "demux_integrated_labeled_obj.RData")


# differential expression comparing cell types ####
DefaultAssay(integrated_labeled) <- "RNA"
# Run the standard workflow for visualization and clustering
integrated_labeled <- ScaleData(integrated_labeled)

unique(Idents(integrated_labeled))


# write function to call time specific genes

find_time_DEG <- function(obj, metadata, metadata_name, test2use, min.pct){
  
  # Generic function to calculate time-specific DEG for an obj within subgroups
  #diveded according to another ident, e.g. cell type. 
  # Inputs:
  # obj: a Seurat object, containing counts and metadata for time and another identity
  # metadata: a named vector, one metadata column from obj
  # metadata_name: a string, the name of the metadata column to be used in FindMarkers()
  # test2use: a string, indicates which test to send to FindMarkers() as input 
  # min.pct: a numerical, minimum percentage of cells expressing the gene, send to FindMarkers()
  # Returns: 
  # df_list_2_return: a list, contains two data frames df_1vs0 and df_6vs0, which
  #contains DEG comparing 1hr vs 0hr and DEG comparing 6hr vs 0hr
  
  Idents(obj) <- metadata_name
  j <- 0
  # calc hour 1 vs 0
  for (i in sort(unique(metadata))) {
    print(i)
    j = j + 1
    if (j == 1) {
      print("first cell type")
      df_1vs0 <- FindMarkers(object = obj,
                             slot = "scale.data",
                             ident.1 = "1hr",
                             group.by = 'time.ident', 
                             subset.ident = i,
                             ident.2 = "0hr",
                             min.pct = min.pct,
                             logfc.threshold = 0.0,
                             test.use = test2use)
      df_1vs0$gene_symbol <- rownames(df_1vs0)
      df_1vs0$cell_type <- i
    } else {
      print("others")
      df_to_append <- FindMarkers(object = obj,
                                  slot = "scale.data",
                                  ident.1 = "1hr",
                                  group.by = 'time.ident',
                                  subset.ident = i,
                                  ident.2 = "0hr",
                                  min.pct = min.pct,
                                  logfc.threshold = 0.0,
                                  test.use = test2use)
      df_to_append$gene_symbol <- rownames(df_to_append)
      df_to_append$cell_type <- i
      df_1vs0 <- rbind(df_1vs0,
                       df_to_append,
                       make.row.names = F,
                       stringsAsFactors = F)
    }
  }
  
  j <- 0
  # calculate 6hr vs 0hr
  for (i in sort(unique(metadata))) {
    print(i)
    j = j + 1
    if (j == 1) {
      print("first cell type")
      df_6vs0 <- FindMarkers(object = obj,
                             slot = "scale.data",
                             ident.1 = "6hr",
                             group.by = 'time.ident',
                             subset.ident = i,
                             ident.2 = "0hr",
                             min.pct = min.pct,
                             logfc.threshold = 0.0,
                             test.use = test2use)
      df_6vs0$gene_symbol <- rownames(df_6vs0)
      df_6vs0$cell_type <- i
    } else {
      print("others")
      df_to_append <- FindMarkers(object = obj,
                                  slot = "scale.data",
                                  ident.1 = "6hr",
                                  group.by = 'time.ident',
                                  subset.ident = i,
                                  ident.2 = "0hr",
                                  min.pct = min.pct,
                                  logfc.threshold = 0.0,
                                  test.use = test2use)
      df_to_append$gene_symbol <- rownames(df_to_append)
      df_to_append$cell_type <- i
      df_6vs0 <- rbind(df_6vs0,
                       df_to_append,
                       make.row.names = F,
                       stringsAsFactors = F)
    }
  }
  df_list_2_return <- vector(mode = "list", length = 2L)
  df_list_2_return[[1]] <- df_1vs0
  df_list_2_return[[2]] <- df_6vs0
  
  return(df_list_2_return)
}

# use MAST, min.pct = 0.0 - 0.05 - 0.1
dfs_pct0 <- find_time_DEG(obj = integrated_labeled,
                          metadata = integrated_labeled$cell.type, 
                          metadata_name = "cell.type",
                          test2use = "MAST", 
                          min.pct = 0.0)
dfs_pct0.05 <- find_time_DEG(obj = integrated_labeled,
                             metadata = integrated_labeled$cell.type, 
                             metadata_name = "cell.type",
                             test2use = "MAST", 
                             min.pct = 0.05)
dfs_pct0.1 <- find_time_DEG(obj = integrated_labeled,
                            metadata = integrated_labeled$cell.type, 
                            metadata_name = "cell.type",
                            test2use = "MAST", 
                            min.pct = 0.1)
dfs_pct0.15 <- find_time_DEG(obj = integrated_labeled,
                            metadata = integrated_labeled$cell.type, 
                            metadata_name = "cell.type",
                            test2use = "MAST", 
                            min.pct = 0.15)

#save.image("DEG_pct0-0.15.RData")
dfs_pct0.01 <- find_time_DEG(obj = integrated_labeled,
                             metadata = integrated_labeled$cell.type, 
                             metadata_name = "cell.type",
                             test2use = "MAST",
                             min.pct = 0.01)
dfs_pct0.02 <- find_time_DEG(obj = integrated_labeled,
                             metadata = integrated_labeled$cell.type, 
                             metadata_name = "cell.type",
                             test2use = "MAST",
                             min.pct = 0.02)
dfs_pct0.03 <- find_time_DEG(obj = integrated_labeled,
                             metadata = integrated_labeled$cell.type, 
                             metadata_name = "cell.type",
                             test2use = "MAST",
                             min.pct = 0.03)
rm("integrated_labeled")
save.image(file = "mDEG_pct0.01-0.03.RData")

# make dotplot
df_to_plot <- rbind(df_1vs0[df_1vs0$gene_symbol %in% c("BDNF", "FOS", "VGF", "NPAS4"), ],
                    df_6vs0[df_6vs0$gene_symbol %in% c("BDNF", "FOS", "VGF", "NPAS4"), ])
df_to_plot$source <- c(rep_len("1v0hr", length.out = 15),
                       rep_len("6v0hr", length.out = 15))


ggplot(df_to_plot,
       aes(x = as.factor(source),
           y = as.factor(cell_type),
           size = pct.1 * 100,
           fill = avg_log2FC)) +
  geom_point(shape = 21) +
  ylab("cell type") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  facet_grid(cols = vars(gene_symbol)) +
  theme_bw()


# check expression in time points
FeaturePlot(integrated_labeled,
            features = c("FOS", "NPAS4", "BDNF", "VGF"))
h0 <- subset(integrated_labeled, time.ident == "0hr")
h1 <- subset(integrated_labeled, time.ident == "1hr")
h6 <- subset(integrated_labeled, time.ident == "6hr")

FeaturePlot(h0, features = c("FOS", "NPAS4"))
FeaturePlot(h1, features = c("FOS", "NPAS4"))
FeaturePlot(h6, features = c("FOS", "NPAS4"))

FeaturePlot(h0, features = c("BDNF", "VGF"))
FeaturePlot(h1, features = c("BDNF", "VGF"))
FeaturePlot(h6, features = c("BDNF", "VGF"))


# test with wilcox, min.pct = 0.1
test1 <- find_time_DEG(obj = integrated_labeled,
                       metadata = integrated_labeled$cell.type, 
                       metadata_name ="cell.type",
                       test2use = "wilcox", 
                       min.pct = 0.1)



# assign case/control idents ####
MGS <- read_excel("MGS_iPSC lines_60samples_scRNA_ATAC-seq_bulkATAC-seq status_AK_HZ.xlsx")
MGS$`new iPSC line ID` <- str_replace(string = MGS$`new iPSC line ID`,
                                      pattern = "00000",
                                      replacement = "_")
MGS <- data.frame(line = MGS$`new iPSC line ID`,
                  aff = MGS$Aff)

integrated_labeled$case.ident <- NA
for (line in unique(integrated_labeled$cell.line.ident)){
  print(line)
  integrated_labeled$case.ident[integrated_labeled$cell.line.ident == line] <-
    MGS$aff[MGS$line == line]
}
unique(integrated_labeled$case.ident)
sum(integrated_labeled$case.ident == "control") #34379
sum(integrated_labeled$case.ident == "case") #33006


# differential expression comparing 3 times points in case and in control ####
case <- subset(integrated_labeled, case.ident == "case")
control <- subset(integrated_labeled, case.ident == "control")

case_dfs <- find_time_DEG(obj = case, 
                          metadata = case$cell.type, 
                          metadata_name = "cell.type", 
                          test2use = "MAST",
                          min.pct = 0)
control_dfs <- find_time_DEG(obj = control, 
                          metadata = control$cell.type, 
                          metadata_name = "cell.type", 
                          test2use = "MAST",
                          min.pct = 0)


# differential expression comparing case and control in 3 times points ####
timelist <- c("0hr", "1hr", "6hr")
objlist <- vector(mode = "list", length = 3L)
for (i in 1:3){
  objlist[[i]] <- subset(integrated_labeled, time.ident == timelist[i])
}
case_df_list <- vector(mode = "list", length = 3L)

for (i in 1:3){
  obj <- objlist[[i]]
  Idents(obj) <- "cell.type"
  j <- 0
  # calc hour 1 vs 0
  for (k in sort(unique(obj$cell.type))) {
    print(k)
    j = j + 1
    if (j == 1) {
      print("first cell type")
      df_temp <- FindMarkers(object = obj,
                             ident.1 = "case",
                             group.by = 'case.ident',
                             subset.ident = k,
                             ident.2 = "control",
                             min.pct = 0.0,
                             logfc.threshold = 0.0,
                             test.use = "MAST")
      df_temp$gene_symbol <- rownames(df_temp)
      df_temp$cell_type <- k
    } else {
      print("others")
      df_to_append <- FindMarkers(object = obj,
                                  ident.1 = "case",
                                  group.by = 'case.ident',
                                  subset.ident = k,
                                  ident.2 = "control",
                                  min.pct = 0.0,
                                  logfc.threshold = 0.0,
                                  test.use = "MAST")
      df_to_append$gene_symbol <- rownames(df_to_append)
      df_to_append$cell_type <- k
      df_temp <- rbind(df_temp,
                       df_to_append,
                       make.row.names = F,
                       stringsAsFactors = F)
    }
  }
  case_df_list[[i]] <- df_temp
} 


#save.image(file = "11Feb2022_find_many_markers.RData")

# output the gene list for GO term analysis
unique(case_df_list[[1]]$cell_type)
filtered_case_control_list_up <- vector(mode = "list", length = 15L)
filtered_case_control_list_do <- vector(mode = "list", length = 15L)
name_nums <- c("0", "1", "6")
j = 0
for (i in 1:3){
  for (t in unique(case_df_list[[i]]$cell_type)){
    j = j + 1
    up_temp <- case_df_list[[i]][case_df_list[[i]]$p_val_adj < 0.05 &
                                   case_df_list[[i]]$cell_type == t &
                                   case_df_list[[i]]$avg_log2FC > 0, ]
    do_temp <- case_df_list[[i]][case_df_list[[i]]$p_val_adj < 0.05 &
                                   case_df_list[[i]]$cell_type == t &
                                   case_df_list[[i]]$avg_log2FC < 0, ]
    if (length(up_temp$gene_symbol) < 100){
      
      up_temp <- case_df_list[[i]][case_df_list[[i]]$cell_type == t &
                                     case_df_list[[i]]$avg_log2FC > 0, ]
      if (length(up_temp$gene_symbol) < 100){
        filtered_case_control_list_up[[j]] <- up_temp
      } else {
        filtered_case_control_list_up[[j]] <- up_temp[1:100, ]
      }
      
      do_temp <- case_df_list[[i]][case_df_list[[i]]$cell_type == t &
                                     case_df_list[[i]]$avg_log2FC < 0, ]
      if (length(do_temp$gene_symbol) < 100){
        filtered_case_control_list_do[[j]] <- do_temp
      } else {
        filtered_case_control_list_do[[j]] <- do_temp[1:100, ]
      }
      
    } else {
      filtered_case_control_list_up[[j]] <- up_temp
      filtered_case_control_list_do[[j]] <- do_temp
    }
    write.table(filtered_case_control_list_up[[j]]$gene_symbol, 
                file = paste0("filtered_deg_expanded_", name_nums[i], "hr_", t, "_up.txt"), 
                quote = F, sep = "\t", row.names = F, col.names = F)
    write.table(filtered_case_control_list_do[[j]]$gene_symbol, 
                file = paste0("filtered_deg_expanded_", name_nums[i], "hr_", t, "_do.txt"), 
                quote = F, sep = "\t", row.names = F, col.names = F)
  }
  
}

# separate up and downregulated genes
for (i in 1:3){
  write.table(filtered_case_control_list[[i]][filtered_case_control_list[[i]]$avg_log2FC > 0, ]$gene_symbol,
              file = paste0("filtered_deg_list_up", name_nums[i], ".txt"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(filtered_case_control_list[[i]][filtered_case_control_list[[i]]$avg_log2FC < 0, ]$gene_symbol,
              file = paste0("filtered_deg_list_do", name_nums[i], ".txt"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
}

for (i in 1:3){
  write.table(case_df_list[[i]][case_df_list[[i]]$avg_log2FC > 0, ],
              file = paste0("unfiltered_deg_df_up", name_nums[i], ".txt"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(case_df_list[[i]][case_df_list[[i]]$avg_log2FC < 0, ],
              file = paste0("unfiltered_deg_list_do", name_nums[i], ".txt"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
}
# volcano plots
zero <- case_df_list[[1]]
one <- case_df_list[[2]]
six <- case_df_list[[3]]
make_volcano
zero$significance <- "nonsig"
zero$significance[zero$p_val_adj < 0.05 & zero$avg_log2FC > 0] <- "pos"
zero$significance[zero$p_val_adj < 0.05 & zero$avg_log2FC < 0] <- "neg"
unique(zero$significance)
zero$neg_log_pval <- (0 - log10(zero$p_val))
zero$labelling <- ""
for (i in c("BDNF", "FOS", "EGR1")){
  zero$labelling[zero$gene_symbol %in% i] <- i
}
# unique(zero$labelling)
for (t in unique(zero$cell_type)){
  print(t)
  p <- ggplot(data = zero[zero$cell_type == t, ],
         aes(x = avg_log2FC, 
             y = neg_log_pval, 
             color = significance,
             label = labelling)) + 
    geom_point(size = 0.2) +
    scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
    theme_minimal() +
    geom_text_repel(box.padding = unit(0.05, 'lines'),
                    min.segment.length = 0,
                    force = 2,
                    max.overlaps = 10000,
                    force_pull = 0.5) +
    ggtitle(paste0("DEG at 0hr for ", t, " cells"))
  
  print(p)
}

# read GO term results
setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/DEG_case_control_by_time")
analysis_0 <- 
  read.delim("analysis_0.txt")
analysis_1 <- 
  read.delim("analysis_1.txt")
analysis_6 <- 
  read.delim("analysis_6.txt")
col_names <- c("term", "REF", "#", "expected", "over_under", 
               "fold_enrichment", "raw_p_value", "FDR")
colnames(analysis_0) <- col_names
colnames(analysis_1) <- col_names
colnames(analysis_6) <- col_names
combined_0_up <- analysis_0[analysis_0$over_under == "+", ]
combined_0_do <- analysis_0[analysis_0$over_under == "-", ]
combined_1_up <- analysis_1[analysis_1$over_under == "+", ]
combined_1_do <- analysis_1[analysis_1$over_under == "-", ]
combined_6_up <- analysis_6[analysis_6$over_under == "+", ]
combined_6_do <- analysis_6[analysis_6$over_under == "-", ]

combined_up <- rbind(combined_0_up[1:15, ], combined_1_up[1:15, ], combined_6_up[1:15, ])
combined_up$time <- c(rep_len("0hr", 15), rep_len("1hr", 15), rep_len("6hr", 15))
combined_up$term <- gsub(pattern = '\\(GO:[0-9]+\\)', replacement = "", x = combined_up$term)

combined_do <- rbind(combined_0_do, combined_1_do, combined_6_do)
combined_do$time <- c(rep_len("0hr", length(combined_0_do$term)), 
                      rep_len("1hr", length(combined_1_do$term)), 
                      rep_len("6hr", length(combined_6_do$term)))
combined_do$term <- gsub(pattern = '\\(GO:[0-9]+\\)', replacement = "", x = combined_do$term)

# plot with dotplot
ggplot(combined_up,
       aes(x = as.factor(time),
           y = as.factor(term),
           size = (0 - log10(raw_p_value)),
           fill = as.numeric(fold_enrichment))) +
  geom_point(shape = 21) +
  labs(title = "Gene set enrichment (upregulated)\nin differentially expressed genes\n comparing case/control", 
       x = "", y = "", 
       fill = "fold\nenrichment \n", size = "p value") +
  scale_fill_gradientn(colors = brewer.pal(7, "YlGnBu")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))

ggplot(combined_do,
       aes(x = as.factor(time),
           y = as.factor(term),
           size = (0 - log10(raw_p_value)),
           fill = as.numeric(fold_enrichment))) +
  geom_point(shape = 21) +
  labs(title = "Gene set enrichment (downregulated)\nin differentially expressed genes\n comparing case/control", 
       x = "", y = "", 
       fill = "fold\nenrichment \n", size = "p value") +
  scale_fill_gradientn(colors = brewer.pal(7, "YlGnBu")) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))

# case/control specific deg analysis
up_list <- list.files(path = ".", pattern = "analysis_up", recursive = F)
for (i in 1:3){
  analysis <- read.delim(up_list[i])
  col_names <- c("term", "REF", "#", "expected", "over_under", 
               "fold_enrichment", "raw_p_value", "FDR")
  colnames(analysis) <- col_names
  analysis$term <- gsub(pattern = '\\(GO:[0-9]+\\)', replacement = "", x = analysis$term)
  analysis$sig <- "nonsig"
  analysis$sig[analysis$over_under == "+" & analysis$FDR < 0.05] <- "up"
  analysis$sig[analysis$over_under == "-" & analysis$FDR < 0.05] <- "down"
  p <- ggplot(analysis,
              aes(x = as.numeric(fold_enrichment),
                  y = (0 - log10(raw_p_value)),
                  fill = sig)) +
    geom_point(shape = 21, colour = "transparent") + 
    scale_fill_manual(values = c("grey", "blue", "red")) +
    theme(text = element_text(size = 10)) +
    theme_light()
  print(p)
}

# check pct.expressed and exp level ####
features <- c("FOS", "BDNF", "VGF", "EGR1", "RGS2")
DefaultAssay(integrated_labeled) <- "SCT"

integrated_labeled$cell.type.time <- NA
for (i in unique(integrated_labeled$cell.type)){
  for (j in unique(integrated_labeled$time.ident)){
    id <- paste(i, j, sep = "_")
    print(id)
    integrated_labeled$cell.type.time[integrated_labeled$cell.type == i &
                                        integrated_labeled$time.ident == j] <- id
  }
  
}
Idents(integrated_labeled)
testobj <- subset(integrated_labeled, cell.type %in% c("GABA", "NEFM_pos_glut", "NEFM_neg_glut"))
DPnew(integrated_labeled, 
      features = features,
      cols = c("darkred", "white"),
      group.by = "cell.type.time") + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        axis.title = element_text(size = 12))

# plot feature plot showing expression levels for several marker genes ####
for(k in unique(integrated_labeled$lib.ident)){
  print(k)
  if (k == "39_0" | k == "22_0" | k == "49_0" ){
    jpeg(file = paste0("feature_plot_", k, ".jpeg"))
    p <- FeaturePlot(subset(integrated_labeled, lib.ident == k), 
                     features = c("FOS", "BDNF", "VGF"), 
                     pt.size = 0.2) 
    print(p)
    dev.off()
  } else {
    jpeg(file = paste0("feature_plot_", k, ".jpeg"))
    p <- FeaturePlot(subset(integrated_labeled, lib.ident == k), 
                     features = c("FOS", "NPAS4", "BDNF", "VGF"), 
                     pt.size = 0.2) 
    print(p)
    dev.off()
  }
}

