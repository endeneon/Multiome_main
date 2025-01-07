library(tidyverse)
library(qvalue)

# Set working directory
setwd("/project/xinhe/zicheng/neuron_stimulation/script/")

# Function to aggregate results and calculate FDR and q-values
aggregate_results <- function(cell_type) {
  result_files <- list.files("/scratch/midway3/zichengwang/combined_time/log_mixed_res/", 
                             pattern = paste0(cell_type, ".*"), full.names = TRUE)
  agg_res <- lapply(result_files, readRDS) %>%
    bind_rows() %>%
    mutate(fdr = p.adjust(P, method = "fdr"), 
           qval = qvalue(P)$qvalues)
  return(agg_res)
}

# Aggregate results for each cell type
cell_types <- c("GABA", "nmglut", "npglut")
p2g_res_list <- lapply(cell_types, aggregate_results)
names(p2g_res_list) <- cell_types

# Save the aggregated results
saveRDS(p2g_res_list, "../output/scmixedGLM/log_Mito0.01_p2g_res_list.rds")
