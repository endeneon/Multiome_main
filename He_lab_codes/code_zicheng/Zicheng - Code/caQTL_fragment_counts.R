# Load necessary libraries
library(tidyverse)
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)

# Set working directory and seed for reproducibility
setwd("/project/xinhe/zicheng/neuron_stimulation/caQTL/script/")
set.seed(1)

# Set up ArchR environment
addArchRThreads(threads = 8) 
addArchRGenome("hg38")

# Load ArchR project
ArchR_AllCells <- loadArchRProject("/scratch/midway3/zichengwang/ArchR_AllCells/")

# Get unique cell clusters
unique_clusters <- unique(ArchR_AllCells$cell_cluster)

# Loop through each cluster
for (cluster in unique_clusters) {
  # Subset ArchR project by cluster
  subset_ArchR <- ArchR_AllCells[ArchR_AllCells$cell_cluster == cluster, ]
  
  # Generate pseudo-bulk count matrix
  pseudo_count_SE <- getGroupSE(subset_ArchR, useMatrix = "PeakMatrix", groupBy = "sample_id", divideN = FALSE)
  pseudo_count_df <- as.data.frame(assay(pseudo_count_SE))
  
  # Create peakset identifiers
  pseudo_count_rowdata <- rowData(pseudo_count_SE)
  peakset <- paste0(pseudo_count_rowdata$seqnames, "-", pseudo_count_rowdata$start, "-", pseudo_count_rowdata$end)
  rownames(pseudo_count_df) <- peakset
  
  # Save the count matrix to an RDS file
  saveRDS(pseudo_count_df, paste0("/project/xinhe/zicheng/neuron_stimulation/caQTL/data/count_matrix/archr/Raw/", cluster, ".rds"))
  
  # Print completion message for the cluster
  print(paste0("Finished ", cluster))
}
