library(tidyverse)
library(vroom)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(gridExtra)
# library(ComplexHeatmap)

# Filter cPeaks
time_vector <- rep(c("0hr", "1hr", "6hr"), 3)
cell_type_vector <- c(rep("GABA", 3), rep("nmglut", 3), rep("npglut", 3))
clusters <- paste0(time_vector, "__", cell_type_vector)
qtl_output_dir <- "/project/xinhe/zicheng/neuron_stimulation/caQTL/output/caQTL_mapping/"

load_cPeaks <- function(cluster) {
  cis_file <- list.files(qtl_output_dir, pattern = paste0("^Cis_", cluster, "__.*_25kb\\.cis_qtl\\.txt\\.gz$"))
  cPeaks <- vroom(paste0(qtl_output_dir, cis_file)) %>%
    filter(qval <= 0.05)
  return(cPeaks)
}

cPeak_list <- lapply(clusters, load_cPeaks)
names(cPeak_list) <- clusters
save(time_vector, cell_type_vector, clusters, cPeak_list, file = paste0(qtl_output_dir, "res_analysis/cPeak_list_25kb.RData"))
saveRDS(cPeak_list, "/project/xinhe/zicheng/share/cPeaks_25kb_fdr0.05.rds")

# Identify number of peaks and shared peaks
load("/project/xinhe/zicheng/neuron_stimulation/caQTL/output/caQTL_mapping/res_analysis/cPeak_list_25kb.RData")
cPeaks_per_cluster <- lapply(cPeak_list, nrow) %>% unlist()

barplot_df <- data.frame(time = time_vector, cell_type = cell_type_vector, num = cPeaks_per_cluster)

barplot_num_cPeaks <- ggplot(barplot_df, aes(x = cell_type, y = num, fill = time)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  xlab("Cell Type") + ylab("# cPeaks") + labs(fill = "Time") +
  scale_fill_brewer(palette = "Paired") + theme_minimal()

ggsave("/project/xinhe/zicheng/neuron_stimulation/caQTL/output/plot/barplot_num_cPeaks.pdf", plot = barplot_num_cPeaks, device = "pdf")

# Convert nominal parquet to RDS format
library(arrow)

parquet2rds <- function(time, cell_type) {
  parquet_file_list <- list.files("/project/xinhe/zicheng/neuron_stimulation/caQTL/output/caQTL_mapping", pattern = paste0("^Cis_Nominal_", time, "__", cell_type, ".*\\.parquet$"))
  stopifnot(length(parquet_file_list) == 22)
  file_prefix <- gsub(".cis_qtl_pairs.chr.*\\.parquet", "", parquet_file_list[1])
  parquet_file_list <- paste0(file_prefix, ".cis_qtl_pairs.chr", 1:22, ".parquet")
  output_name <- paste0(file_prefix, ".rds")
  
  nominal_by_chr <- lapply(parquet_file_list, function(x) read_parquet(paste0("/project/xinhe/zicheng/neuron_stimulation/caQTL/output/caQTL_mapping/", x))) %>% bind_rows()
  print(paste0(time, " ", cell_type, " parquet files were loaded!"))
  saveRDS(nominal_by_chr, paste0("/project/xinhe/zicheng/neuron_stimulation/caQTL/output/caQTL_mapping/", output_name))
  print(paste0(time, " ", cell_type, " RDS file was saved!"))
  return("Done!")
}

for (time in c("0hr", "1hr", "6hr")) {
  for (cell_type in c("GABA", "nmglut", "npglut")) {
    parquet2rds(time, cell_type)
  }
}

system("rm /project/xinhe/zicheng/neuron_stimulation/caQTL/output/caQTL_mapping/*.parquet")

