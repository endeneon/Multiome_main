library(ArchR)
library(tidyverse)
library(data.table)
library(gtools)

addArchRThreads(16)

# Load and filter cell metadata
subset_cellMeta <- readRDS("/project/xinhe/zicheng/neuron_stimulation/data/cellMeta.rds") %>% 
  as.data.frame() %>% 
  filter(grepl("Batch_024", batch_wo_time))

# Load ArchR project and update sample names
ArchR_subset <- loadArchRProject("/scratch/midway3/zichengwang/ArchR_subset/")
ArchR_subset$batch_wo_time <- gsub("_.$", "", ArchR_subset$Sample)

# Subset ArchR project
ArchR_subset <- ArchR_subset[rownames(subset_cellMeta),]

# Save the subsetted ArchR project
saveArchRProject(ArchR_subset, dropCells = TRUE)

# Update project metadata
arrow_proj <- readRDS("/scratch/midway3/zichengwang/ArchR_subset/Save-ArchR-Project.rds")
arrow_proj@sampleColData <- arrow_proj@sampleColData[grepl("Batch_024", rownames(arrow_proj@sampleColData), fixed = TRUE), , drop = FALSE]
arrow_proj@sampleMetadata@listData <- arrow_proj@sampleMetadata@listData[grepl("Batch_024", names(arrow_proj@sampleMetadata@listData), fixed = TRUE)]
saveRDS(arrow_proj, "/scratch/midway3/zichengwang/ArchR_subset/Save-ArchR-Project.rds")

# Extract fragments from arrow files
ArchR_subset <- loadArchRProject("/scratch/midway3/zichengwang/ArchR_subset/")
ArchR_subset$pseudo_bulk_id <- paste0(ArchR_subset$cell_cluster, "__", ArchR_subset$sample_id)
cell_meta <- getCellColData(ArchR_subset)

generate_tagAlign <- function(dt) {
  dt[, mid := as.integer((V2 + V3) / 2)]  # Calculate mid-point
  
  # Create two new data.tables for each row
  dt_left <- dt[, .(V1, V2, mid, "N", 1000, "+")]
  dt_right <- dt[, .(V1, mid + 1L, V3, "N", 1000, "-")]
  
  # Combine left and right tables
  tagAlign <- rbind(dt_left, dt_right, use.names = FALSE)
  
  return(tagAlign)
}

unique_batch <- unique(cell_meta$Sample)
threads <- 20

lapply(unique_batch, function(batch) {
  cell_frags <- getFragmentsFromArrow(ArrowFile = paste0("/scratch/midway3/zichengwang/ArchR_subset/ArrowFiles/", batch, ".arrow"), verbose = FALSE)
  subset_cell_meta <- cell_meta[cell_meta$Sample == batch, ]
  subset_cellnames <- rownames(subset_cell_meta)
  
  lapply(c("GABA", "nmglut", "npglut"), function(cell_type) {
    context_frags <- cell_frags[cell_frags$RG %in% subset_cellnames[subset_cell_meta$cell_type == cell_type]] %>% 
      as.data.table() %>% 
      .[, 1:3] %>% 
      .[, start := start - 1L]
    
    setnames(context_frags, c("V1", "V2", "V3"))
    
    context_tagAlign <- generate_tagAlign(context_frags)
    
    output_file <- paste0("/scratch/midway3/zichengwang/tagAlign_files/", batch, "_", cell_type, "_tagAlign.gz")
    con <- pipe(paste0("LC_ALL=C sort -k 1,1V -k 2,2n -k3,3n -T /scratch/midway3/zichengwang/temp_dir -S 5G --parallel ", threads,
                       " | /scratch/midway2/zichengwang/software/htslib/bin/bgzip -@ ", threads, " -c > ", output_file), "wb")
    write_tsv(context_tagAlign, file = con, col_names = FALSE)
    close(con)
    
    system(paste0("/scratch/midway2/zichengwang/software/htslib/bin/tabix -@ ", threads, " -p bed ", output_file))
    
    cat("Finished", batch, cell_type, "\n")
    return("Done")
  })
})

# ABC snakemake config
config <- read_tsv("/scratch/midway3/zichengwang/Github_Repo/ABC-Enhancer-Gene-Prediction/config/config_biosamples_template.tsv")
config[1:9, ] <- NA

config$biosample <- paste0(rep(c("0hr", "1hr", "6hr"), 3), "_", c(rep("GABA", 3), rep("nmglut", 3), rep("npglut", 3)))

config$ATAC <- sapply(paste0(rep(c(0, 1, 6), 3), "_", c(rep("GABA", 3), rep("nmglut", 3), rep("npglut", 3))), function(x) {
  list.files("/scratch/midway3/zichengwang/tagAlign_files", pattern = paste0(".*", x, "_tagAlign.gz$"), full.names = TRUE) %>% paste0(., collapse = ",")
})
config$default_accessibility_feature <- "ATAC"
config$HiC_file <- paste0("/project2/xinhe/zicheng/archive/hicFiles/CD11-", rep(c(0, 1, 6), 3), ".hic")

config$HiC_type <- "hic"
config$HiC_resolution <- 5000

write_tsv(config, "/scratch/midway3/zichengwang/Github_Repo/ABC-Enhancer-Gene-Prediction/config/microc_config_biosamples.tsv")

# Add UMAP embedding
addArchRThreads(8)
ArchR_subset <- loadArchRProject("/scratch/midway3/zichengwang/ArchR_subset/")

# Add chromatin LSI embedding
ArchR_subset <- addIterativeLSI(
  ArchRProj = ArchR_subset,
  useMatrix = "TileMatrix", 
  name = "LSI_ATAC", 
  iterations = 2, 
  force = TRUE,
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  )
)

# Correct for batch effects
ArchR_subset <- addHarmony(
  ArchRProj = ArchR_subset,
  reducedDims = "LSI_ATAC",
  groupBy = "batch_wo_time",
  name = "Harmony_ATAC",
  force = TRUE
)

# Add expression LSI embedding
ArchR_subset <- addIterativeLSI(
  ArchRProj = ArchR_subset, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA",
  force = TRUE
)

# Correct for batch effects
ArchR_subset <- addHarmony(
  ArchRProj = ArchR_subset,
  reducedDims = "LSI_RNA",
  groupBy = "batch_wo_time",
  name = "Harmony_RNA",
  force = TRUE
)

# Combine the two embeddings
ArchR_subset <- addCombinedDims(ArchR_subset, reducedDims = c("Harmony_ATAC", "Harmony_RNA"), name = "Harmony_Combined")
ArchR_subset <- addUMAP(ArchR_subset, reducedDims = "Harmony_Combined", name = "UMAP_Harmony_Combined", minDist = 0.8, force = TRUE)

plotEmbedding(ArchR_subset, name = "cell_cluster", embedding = "UMAP_Harmony_Combined", size = 1.5, labelAsFactors = FALSE, labelMeans = FALSE, width = 5, height = 5)

saveArchRProject(ArchR_subset)

# Pseudo time trajectory
addArchRThreads(8)
ArchR_subset <- loadArchRProject("/scratch/midway3/zichengwang/ArchR_subset/")

time_points <- c("0hr", "1hr", "6hr")
cell_types <- c("GABA", "nmglut", "npglut")

# GABA trajectory
cell_type <- "GABA"
trajectory <- paste0(cell_type, "__", time_points)


###Infer trajectories and plot Fig S6D

ArchR_subset <- addTrajectory(
  ArchRProj = ArchR_subset,
  name = paste0(cell_type, "U"),
  groupBy = "cell_cluster",
  trajectory = trajectory,
  embedding = "UMAP_Harmony_Combined",
  spar = 1.5,
  force = TRUE
)
p <- plotTrajectory(ArchR_subset, trajectory = "GABAU", colorBy = "cellColData", name = "GABAU", embedding = "UMAP_Harmony_Combined")
p[[1]]

# nmglut trajectory
cell_type <- "nmglut"
trajectory <- paste0(cell_type, "__", time_points)

ArchR_subset <- addTrajectory(
  ArchRProj = ArchR_subset,
  name = paste0(cell_type, "U"),
  groupBy = "cell_cluster",
  trajectory = trajectory,
  embedding = "UMAP_Harmony_Combined",
  preFilterQuantile = 0.7,
  postFilterQuantile = 0.7,
  spar = 1.3,
  force = TRUE
)
p <- plotTrajectory(ArchR_subset, trajectory = "nmglutU", colorBy = "cellColData", name = "nmglutU", embedding = "UMAP_Harmony_Combined")
p[[1]]

# npglut trajectory
cell_type <- "npglut"
trajectory <- paste0(cell_type, "__", time_points)

ArchR_subset <- addTrajectory(
  ArchRProj = ArchR_subset,
  name = paste0(cell_type, "U"),
  groupBy = "cell_cluster",
  trajectory = trajectory,
  embedding = "UMAP_Harmony_Combined",
  spar = 1.3,
  force = TRUE
)
p <- plotTrajectory(ArchR_subset, trajectory = "npglutU", colorBy = "cellColData", name = "npglutU", embedding = "UMAP_Harmony_Combined")
p[[1]]

# Add ChromVar analysis
ArchR_subset <- addMotifAnnotations(
  ArchRProj = ArchR_subset, 
  motifSet = "cisbp", 
  name = "Motif", 
  force = TRUE
)

ArchR_subset <- addBgdPeaks(
  ArchRProj = ArchR_subset, 
  force = TRUE, 
  nIterations = 100, 
  seed = 88
)

ArchR_subset <- addDeviationsMatrix(
  ArchRProj = ArchR_subset, 
  peakAnnotation = "Motif", 
  force = TRUE, 
  threads = 8
)

saveArchRProject(ArchR_subset)
