# Chuxuan Li 11/12/2021
# Using ArchR to analyze group 2 and group 8 ATACseq data


# init
library(Seurat)
library(gplots)
library(ArchR)
library(future)
library(stringr)
library(pheatmap)

# set threads and parallelization
plan("multisession", workers = 2)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 11474836480)

addArchRThreads(threads = 10) 
addArchRGenome("hg38")

#inputFile <- "/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/atac_fragments.tsv.gz"

inputFile <- "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_0/outs/atac_fragments.tsv.gz"
names(inputFile) <- "scATAC_10x_aggregated"

inputFile

# create arrow files ####
ArrowFiles <- createArrowFiles(
  inputFile = inputFile,
  sampleNames = names(inputFile),
  minTSS = 3, #Dont set this too high because you can always increase later
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = F,
  #TileMatParams = list(tileSize = 5000),
  addGeneScoreMat = TRUE,
  force = TRUE
)
ArrowFiles

# unfiltered_ArrowFiles <- createArrowFiles(
#   inputFile = inputFile,
#   sampleNames = names(inputFile),
#   minTSS = 0, #Dont set this too high because you can always increase later
#   minFrags = 0,
#   addTileMat = TRUE,
#   addGeneScoreMat = TRUE,
#   force = T
# )

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)


projAggr <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ArchR_10x_aggr",
  copyArrows = TRUE #This is recommended so that if you modify the Arrow files you have an original copy for later usage.
)
projAggr

unfiltered_projAggr <- ArchRProject(
  ArrowFiles = unfiltered_ArrowFiles, 
  outputDirectory = "ArchR_10x_aggr",
  copyArrows = TRUE #This is recommended so that if you modify the Arrow files you have an original copy for later usage.
)


# QC plots ####
df <- getCellColData(projAggr, select = c("log10(nFrags)", "TSSEnrichment"))
df

ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 1, lty = "dashed") + geom_vline(xintercept = 0, lty = "dashed")

unfiltered_df <- getCellColData(unfiltered_projAggr, select = c("log10(nFrags)", "TSSEnrichment"))
unfiltered_df

ggPoint(
  x = unfiltered_df[,1], 
  y = unfiltered_df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(unfiltered_df[,1], probs = 0.99)),
  ylim = c(0, quantile(unfiltered_df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")


p1 <- plotGroups(
  ArchRProj = projAggr, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)

p2 <- plotGroups(
  ArchRProj = projAggr, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

p3 <- plotGroups(
  ArchRProj = projAggr, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)

p4 <- plotGroups(
  ArchRProj = projAggr, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotPDF(p1,p2,p3,p4, 
        name = "QC_stats_filtered.pdf", 
        ArchRProj = projAggr, 
        addDOC = FALSE, 
        width = 4, 
        height = 4)

p1 <- plotGroups(
  ArchRProj = unfiltered_projAggr, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)

p2 <- plotGroups(
  ArchRProj = unfiltered_projAggr, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

p3 <- plotGroups(
  ArchRProj = unfiltered_projAggr, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)

p4 <- plotGroups(
  ArchRProj = unfiltered_projAggr, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
plotPDF(p1,p2,p3,p4, 
        name = "QC_stats_unfiltered.pdf", 
        ArchRProj = unfiltered_projAggr, 
        addDOC = FALSE, 
        width = 4, 
        height = 4)

saveArchRProject(ArchRProj = projAggr, outputDirectory = "projAggr", load = FALSE)
saveArchRProject(ArchRProj = unfiltered_projAggr, outputDirectory = "unfiltered_projAggr", load = FALSE)

# Dimensional reduction ####
projAggr_test_0 <- addIterativeLSI(
  ArchRProj = projAggr,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 0, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = 0.8,
    maxClusters = NULL
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:50,
  totalFeatures = 1e+06,
  force = TRUE,
  verbose = T
)


# clustering
projAggr_test_0 <- addClusters(
  input = projAggr_test_0,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  force = TRUE,
  resolution = 0.5
)
# understand where each group is in which cluster
# cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
# cM

# run UMAP
projAggr_test_0 <- addUMAP(
  ArchRProj = projAggr_test_0, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)
# check clustering
plotEmbedding(ArchRProj = projAggr_test_0, 
              colorBy = "cellColData", 
              name = "Clusters", 
              embedding = "UMAP")

new_markers <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
                 "EBF1", "SERTAD4", # striatal
                 "FOXG1", # forebrain 
                 "FOXP2", "TLE4", # pallial glutamatergic
                 "ADCYAP1", "TCERG1L", "SEMA3E", # subcerebral
                 "POU3F2", "CUX1", "BCL11B",  # cortical
                 "LHX2", # general cortex
                 "SST", # inhibitory
                 "SATB2", "CUX2", "NEFM", # excitatory
                 "VIM", "SOX2", "NES", #NPC
                 "MAP2", "DCX"
)
p <- plotEmbedding(
  ArchRProj = projAggr_test_0, 
  colorBy = "GeneScoreMatrix", 
  name = new_markers, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

p$GAD1
p$SLC17A6
p$SOX2
p$VIM
p$NEFM

# use one iteration to try again
projAggr_test_1 <- addIterativeLSI(
  ArchRProj = projAggr,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 1, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = 0.8,
    maxClusters = NULL
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:50,
  totalFeatures = 1e+06,
  force = TRUE,
  verbose = T
)


# clustering
projAggr_test_1 <- addClusters(
  input = projAggr_test_1,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  force = TRUE,
  resolution = 0.5
)
# understand where each group is in which cluster
# cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
# cM

# run UMAP
projAggr_test_1 <- addUMAP(
  ArchRProj = projAggr_test_1, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)
# check clustering
plotEmbedding(ArchRProj = projAggr_test_1, 
              colorBy = "cellColData", 
              name = "Clusters", 
              embedding = "UMAP")


p <- plotEmbedding(
  ArchRProj = projAggr_test_1, 
  colorBy = "GeneScoreMatrix", 
  name = new_markers, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

p$GAD1
p$SLC17A6
p$SOX2
p$VIM


# use 5 iterations to try again
projAggr_test_5 <- addIterativeLSI(
  ArchRProj = projAggr,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 5, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = 0.8,
    maxClusters = NULL
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:50,
  totalFeatures = 1e+06,
  force = TRUE,
  verbose = T
)
# clustering
projAggr_test_5 <- addClusters(
  input = projAggr_test_5,
  reducedDims = "IterativeLSI",
  method = "Seurat", 
  name = "Clusters",
  force = TRUE,
  resolution = 0.8
)
# understand where each group is in which cluster
# cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
# cM

# run UMAP
projAggr_test_5 <- addUMAP(
  ArchRProj = projAggr_test_5, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)
# check clustering
plotEmbedding(ArchRProj = projAggr_test_5, 
              colorBy = "cellColData", 
              name = "Clusters", 
              embedding = "UMAP")


p <- plotEmbedding(
  ArchRProj = projAggr_test_5, 
  colorBy = "GeneScoreMatrix", 
  name = new_markers, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

p$GAD1
p$SLC17A6
p$SOX2
p$VIM



# use 10 iterations to try again
projAggr_test_10 <- addIterativeLSI(
  ArchRProj = projAggr,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 10, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = 0.8,
    maxClusters = NULL
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:50,
  totalFeatures = 1e+06,
  force = TRUE,
  verbose = T
)
# clustering
projAggr_test_10 <- addClusters(
  input = projAggr_test_10,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  force = TRUE,
  resolution = 0.5
)
# understand where each group is in which cluster
# cM <- confusionMatrix(paste0(projHeme2$Clusters), paste0(projHeme2$Sample))
# cM

# run UMAP
projAggr_test_10 <- addUMAP(
  ArchRProj = projAggr_test_10, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)
# check clustering
plotEmbedding(ArchRProj = projAggr_test_10, 
              colorBy = "cellColData", 
              name = "Clusters", 
              embedding = "UMAP")


p <- plotEmbedding(
  ArchRProj = projAggr_test_10, 
  colorBy = "GeneScoreMatrix", 
  name = new_markers, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

p$GAD1
p$SLC17A6
p$SOX2
p$VIM



markersGS <- getMarkerFeatures(
  ArchRProj = projAggr, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.2")

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.2", 
  labelMarkers = new_markers,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
