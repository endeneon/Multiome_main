# Chuxuan Li 10/25/2021
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
options(future.globals.maxSize = 21474836480)

addArchRThreads(threads = 16) 
addArchRGenome("hg38")

# create arrow files
inputFiles <- c("/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_0/outs/atac_fragments.tsv.gz",
                "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_1/outs/atac_fragments.tsv.gz",
                "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_2_6/outs/atac_fragments.tsv.gz",
                "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_0/outs/atac_fragments.tsv.gz",
                "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_1/outs/atac_fragments.tsv.gz",
                "/home/cli/NVME/scARC_Duan_018/GRCh38_mapped_only/libraries_8_6/outs/atac_fragments.tsv.gz")
names(inputFiles) <- c("g_2_0", "g_2_1", "g_2_6",
                       "g_8_0", "g_8_1", "g_8_6")
inputFiles

#
# "For each sample, this step will:
#   
# Read accessible fragments from the provided input files.
# Calculate quality control information for each cell (i.e. TSS enrichment scores and nucleosome info).
# Filter cells based on quality control parameters.
# Create a genome-wide TileMatrix using 500-bp bins.
# Create a GeneScoreMatrix using the custom geneAnnotation that was defined when we called addArchRGenome()."
#
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, 
  filterFrags = 3200, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
ArrowFiles
#"g_8_0.arrow" "g_2_1.arrow" "g_2_0.arrow" "g_8_6.arrow" "g_8_1.arrow" "g_2_6.arrow"

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

proj_g2_g8 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ArchR",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
proj_g2_g8
# class: ArchRProject 
# outputDirectory: /nvmefs/scARC_Duan_018/R_scARC_Duan_018/ArchR 
# samples(6): g_8_0 g_2_1 ... g_8_1 g_2_6
# sampleColData names(1): ArrowFiles
# cellColData names(15): Sample TSSEnrichment ... DoubletEnrichment BlacklistRatio
# numberOfCells(1): 86066
# medianTSS(1): 7.658
# medianFrags(1): 12149
getAvailableMatrices(proj_g2_g8)
# "GeneScoreMatrix" "TileMatrix"


proj_g2_g8 <- filterDoublets(proj_g2_g8)

# some QC plots
p1 <- plotGroups(
  ArchRProj = proj_g2_g8, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges",
  title = "TSS enrichment",
)
p1
p2 <- plotGroups(
  ArchRProj = proj_g2_g8, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  title = "TSS enrichment",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p2

p3 <- plotGroups(
  ArchRProj = proj_g2_g8, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges",
  title = "log10(unique nuclear fragments)"
)
p3
p4 <- plotGroups(
  ArchRProj = proj_g2_g8, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE,
  title = "log10(unique nuclear fragments)"
)
p4
plotPDF(p1,p2,p3,p4, name = "QC_Statistics.pdf", ArchRProj = proj_g2_g8, addDOC = FALSE, width = 10, height = 10)

p1 <- plotFragmentSizes(ArchRProj = proj_g2_g8)
p1
p2 <- plotTSSEnrichment(ArchRProj = proj_g2_g8)
p2
plotPDF(p1,p2, name = "QC_FragSizes_TSSProfile.pdf", ArchRProj = proj_g2_g8, addDOC = FALSE, width = 10, height = 10)

# dimensionality reduction
proj_g2_g8_LSI <- addIterativeLSI(
  ArchRProj = proj_g2_g8,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

proj_g2_g8_LSI <- addHarmony(
  ArchRProj = proj_g2_g8_LSI,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  theta = 0
)

proj_g2_g8_LSI <- addClusters(
  input = proj_g2_g8_LSI,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5
)

table(proj_g2_g8_LSI$Clusters)

# check which samples are in which clusters with confusion matrix
cM <- confusionMatrix(paste0(proj_g2_g8_LSI$Clusters), paste0(proj_g2_g8_LSI$Sample))
cM
cM_rearr <- cM[, order(colnames(cM))]
cM_rearr <- cM[order(as.numeric(str_sub(rownames(cM), start = 2L))), ]

cM_rearr <- cM_rearr / Matrix::rowSums(cM_rearr)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM_rearr), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)


# UMAP
proj_g2_g8_UMAP <- addUMAP(
  ArchRProj = proj_g2_g8_LSI, 
  reducedDims = "Harmony", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

# plot dimension plots
# color by sample
p1 <- plotEmbedding(ArchRProj = proj_g2_g8_UMAP, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
# color by cluster
p2 <- plotEmbedding(ArchRProj = proj_g2_g8_UMAP, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h") # horizontally align
plotPDF(p1,p2, name = "UMAP_by_Sample_and_Clusters.pdf", ArchRProj = proj_g2_g8_UMAP, addDOC = FALSE, width = 12, height = 5)


