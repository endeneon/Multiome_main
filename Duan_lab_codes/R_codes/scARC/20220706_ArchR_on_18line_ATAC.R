# Chuxuan Li 07/06/2022
# Use ArchR to generate gene activity score for heatmap

# init ####
library(ArchR)
library(Seurat)
library(Signac)
library(colorRamps)
library(stringr)

addArchRThreads(threads = 32)
addArchRGenome("hg38")

# input file ####
inputFiles <- "/data/FASTQ/Duan_Project_024/hybrid_output/hybrid_aggr_5groups_no2ndary/outs/atac_fragments.tsv.gz"
names(inputFiles) <- "aggr"
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 0, #Dont set this too high because you can always increase later
  filterFrags = 0, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "./",
  copyArrows = F
)
getAvailableMatrices(proj)
saveArchRProject(ArchRProj = proj, outputDirectory = "./", load = FALSE)

proj <- loadArchRProject(path = "./")
load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/multiomic_obj.RData")

# add metadata ####
# 1. add time point
times <- unique(multiomic_obj$time.ident)
proj$time.ident <- "NA"
for (t in times) {
  obj <- subset(multiomic_obj, time.ident == t)
  proj$time.ident[str_remove(proj$cellNames, "aggr#") %in% obj@assays$ATAC@counts@Dimnames[[2]]] <- t
}
unique(proj$time.ident)
# 2. add cell type
types <- unique(multiomic_obj$cell.type)
proj$cell.type <- "NA"
for (t in types) {
  obj <- subset(multiomic_obj, cell.type == t)
  proj$cell.type[str_remove(proj$cellNames, "aggr#") %in% obj@assays$ATAC@counts@Dimnames[[2]]] <- t
}
unique(proj$cell.type)
# 3. add cluster
atac.cl <- unique(multiomic_obj$seurat_clusters)
proj$ATAC.cluster.ident <- "NA"
for (c in atac.cl) {
  print(c)
  obj <- subset(multiomic_obj, seurat_clusters == c)
  proj$ATAC.cluster.ident[str_remove(proj$cellNames, "aggr#") %in% obj@assays$ATAC@counts@Dimnames[[2]]] <- c
}
unique(proj$ATAC.cluster.ident)
# 4. add group
groups <- unique(multiomic_obj$group.ident)
proj$group.ident <- "NA"
for (g in groups) {
  print(g)
  obj <- subset(multiomic_obj, group.ident == g)
  proj$group.ident[str_remove(proj$cellNames, "aggr#") %in% obj@assays$ATAC@counts@Dimnames[[2]]] <- g
}
unique(proj$group.ident)

# get gene matrix ####
idxType <- BiocGenerics::which(proj$cell.type != "NA")
cells <- proj$cellNames[idxType]
proj[cells]
proj_sub <- proj[cells]
gsse <- getMatrixFromProject(proj_sub)
gsmat <- gsse@assays@data@listData$GeneScoreMatrix
rownames(gsmat) <- gsse@elementMetadata@listData$name

# follow through ArchR pipeline ####
proj_sub <- addIterativeLSI(
  ArchRProj = proj_sub,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.5), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 30000, 
  dimsToUse = 1:30
)
proj_sub <- addHarmony(
  ArchRProj = proj_sub,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "group.ident"
)
proj_sub <- addClusters(
  input = proj_sub,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

table(proj_sub$Clusters)
proj_sub <- addUMAP(
  ArchRProj = proj_sub, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
plotEmbedding(
  ArchRProj = proj_sub, 
  colorBy = "GeneScoreMatrix", 
  name = genelist, 
  embedding = "UMAP", 
  bins = 100,
  continuousSet = "horizonExtra", 
  quantileCut = c(0, 0.999),
  log2norm = TRUE
)
markerHeatmap(
  seMarker = genelist, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = genelist,
  transpose = TRUE
)

# create gact object for heatmap ####
obj_sub <- multiomic_obj
obj_sub$in.archr <- "no"
obj_sub$in.archr[colnames(obj_sub) %in% str_remove(gsmat@Dimnames[[2]], "aggr#")] <- "yes"
sum(obj_sub$in.archr == "yes")
obj_sub <- subset(obj_sub, in.archr == "yes")
gsmat@Dimnames[[2]] <- str_remove(gsmat@Dimnames[[2]], "aggr#")
Idents(obj_sub)
obj_sub[["gact"]] <- CreateAssayObject(counts = gsmat)
obj_sub <- NormalizeData(object = obj_sub, assay = "gact")

# make genelist ####
DEG_paths <- list.files("../Analysis_RNAseq_pseudobulk/5line_DE_results/filtered_by_basemean/significant/", 
                        full.names = T, pattern = ".csv")
DEG_paths <- DEG_paths[str_detect(DEG_paths, "NPC", T)]
DEG_paths_up <- DEG_paths[str_detect(DEG_paths, "upregulated")]
DEG_paths_do <- DEG_paths[str_detect(DEG_paths, "downregulated")]
DEG_lists_up <- vector("list", length(DEG_paths_up))
DEG_lists_do <- vector("list", length(DEG_paths_do))

for (i in 1:length(DEG_lists_up)) {
  DEG_lists_up[[i]] <- read.csv(DEG_paths_up[i], row.names = 1)
  DEG_lists_do[[i]] <- read.csv(DEG_paths_do[i], row.names = 1)
  DEG_lists_up[[i]] <- arrange(DEG_lists_up[[i]], q_value)
  DEG_lists_do[[i]] <- arrange(DEG_lists_do[[i]], q_value)
}
for (i in 1:length(DEG_lists_up)) {
  if (i == 1) {
    genelist_up <- row.names(DEG_lists_up[[i]])
    genelist_do <- row.names(DEG_lists_do[[i]])
  } else {
    genelist_up <- c(genelist_up, row.names(DEG_lists_up[[i]]))
    genelist_do <- c(genelist_do, row.names(DEG_lists_do[[i]]))
  }
}
genelist_up <- unique(genelist_up)
genelist_do <- unique(genelist_do)
genelist_up <- genelist_up[1:1000]
genelist_do <- genelist_do[1:1000]
genelist <- c(genelist_up, genelist_do)
genelist <- genelist[genelist %in% gsmat@elementMetadata@listData$name]

# make heatmaps ####
obj_sub$clusterxtime.ident <- "others"
cls <- sort(unique(obj_sub$seurat_clusters))
for (k in 1:length(cls)) {
  cl <- cls[k]
  for (l in 1:length(times)) {
    t <- times[l]
    obj_sub$clusterxtime.ident[obj_sub$seurat_clusters == cl & 
                                 obj_sub$time.ident == t] <- paste(cl, t, sep = "_")
  }
}
cts <- sort(unique(obj_sub$clusterxtime.ident))

for (j in 1:length(cts)) {
  obj <- subset(obj_sub, clusterxtime.ident == cts[j])
  for (i in 1:length(genelist)) {
    g <- genelist[i]
    if (i == 1) {
      DefaultAssay(obj) <- "SCT"
      gex_col <- sum(exp(obj@assays$SCT@data[rownames(obj) == g, ])) /
        ncol(obj)
      DefaultAssay(obj) <- "gact"
      gact_col <- sum(exp(obj@assays$gact@data[rownames(obj) == g, ])) /
        ncol(obj)
    } else {
      DefaultAssay(obj) <- "SCT"
      gex_col_app <- sum(exp(obj@assays$SCT@data[rownames(obj) == g, ])) /
        ncol(obj)
      gex_col <- rbind(gex_col, gex_col_app)
      DefaultAssay(obj) <- "gact"
      gact_col_app <- sum(exp(obj@assays$gact@data[rownames(obj) == g, ])) /
        ncol(obj)
      gact_col <- rbind(gact_col, gact_col_app)
    }
  }
  rownames(gex_col) <- genelist
  rownames(gact_col) <- genelist
  if (j == 1) {
    print("entered if")
    gex_mat <- gex_col
    gact_mat <- gact_col
  } else {
    gex_mat <- cbind(gex_mat, gex_col)
    gact_mat <- cbind(gact_mat, gact_col)
  }
  print(j)
}
gact_df_nonzero <- gact_mat[rowSums(gact_mat) != 0, ]
gex_df_nonzero <- gex_mat[rowSums(gex_mat) != 0, ]
gex_clust <- hclust(dist(gex_df_nonzero, method = "euclidean"))
heatmap.2(x = as.matrix(gex_df_nonzero),
          scale = "row", 
          trace = "none", 
          #Rowv = T,
          Rowv = gex_clust$order, 
          Colv = T, 
          dendrogram = "column",
          col = "bluered",
          main = "gene expression")
p <- heatmap.2(x = as.matrix(gact_df_nonzero),
               scale = "row", 
               trace = "none", 
               Rowv = gex_clust$order, 
               Colv = T, 
               dendrogram = "none",
               col = colorRampPalette(c("white", "blue", "blue", "blue", "darkblue",
                                        "darkblue", "black", "black", "red", "darkred"))(100),
               symbreaks = F,
               main = "gene activity")
