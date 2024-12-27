# Chuxuan Li 03/22/2023
# Generate aggregation csv for 10x arc aggr for 018-029 integration

library(stringr)

# library_id ####
library_id_018 <- unique(str_extract(list.dirs("/nvmefs/scARC_Duan_018/GRCh38_mapped_only",
                                               recursive = T), "[0-9]+_[0|1|6]"))
library_id <- c(library_id_018[2:7])
library_id_022 <- unique(str_extract(list.dirs("/nvmefs/scARC_Duan_022/cellranger_output",
                                    recursive = T), "[0-9]+-[0|1|6]"))
library_id <- c(library_id, library_id_022[2:16])
library_id_024 <- unique(str_extract(list.dirs("/nvmefs/scARC_Duan_024_raw_fastq/GRCh38_output",
                             recursive = T), "[0-9]+-[0|1|6]"))
library_id <- c(library_id, library_id_024[2:16])
library_id_025 <- unique(str_extract(list.dirs("/nvmefs/scARC_Duan_025_GRCh38",
                                        recursive = T), "[0-9][0-9]-[0|1|6]"))
library_id <- c(library_id, library_id_025[2:16])
library_id_029 <- str_extract(list.dirs("/data/FASTQ/Duan_Project_029/GRCh38_only", 
                                        recursive = F), "[0-9]+-[0|1|6]")
library_id <- c(library_id, library_id_029[2:16])
library_id

# atac_fragments ####
atac_fragments_018 <- list.files("/nvmefs/scARC_Duan_018/GRCh38_mapped_only", 
                                 "atac_fragments.tsv.gz$", full.names = T, 
                                 recursive = T, include.dirs = F)
atac_fragments <- c(atac_fragments_018)
atac_fragments_022 <- list.files("/nvmefs/scARC_Duan_022/cellranger_output", 
                                 "atac_fragments.tsv.gz$", full.names = T, 
                                 recursive = T, include.dirs = F)
atac_fragments <- c(atac_fragments, atac_fragments_022)
atac_fragments_024 <- list.files("/nvmefs/scARC_Duan_024_raw_fastq/GRCh38_output", 
                                 "atac_fragments.tsv.gz$", full.names = T, 
                                 recursive = T, include.dirs = F)
atac_fragments <- c(atac_fragments, atac_fragments_024)
atac_fragments_025 <- list.files("/nvmefs/scARC_Duan_025_GRCh38", 
                                 "atac_fragments.tsv.gz$", full.names = T, 
                                 recursive = T, include.dirs = F)
atac_fragments <- c(atac_fragments, atac_fragments_025[2:16])
atac_fragments_029 <- list.files("/data/FASTQ/Duan_Project_029/GRCh38_only", 
                             "atac_fragments.tsv.gz$", full.names = T, 
                             recursive = T, include.dirs = F)
atac_fragments <- c(atac_fragments, atac_fragments_029[1:15])
atac_fragments

# per_barcode_metrics ####
per_barcode_metrics_018 <- list.files("/nvmefs/scARC_Duan_018/GRCh38_mapped_only", 
                                 "per_barcode_metrics.csv", full.names = T, 
                                 recursive = T, include.dirs = F)
per_barcode_metrics <- c(per_barcode_metrics_018)
per_barcode_metrics_022 <- list.files("/nvmefs/scARC_Duan_022/cellranger_output", 
                                 "per_barcode_metrics.csv", full.names = T, 
                                 recursive = T, include.dirs = F)
per_barcode_metrics <- c(per_barcode_metrics, per_barcode_metrics_022)
per_barcode_metrics_024 <- list.files("/nvmefs/scARC_Duan_024_raw_fastq/GRCh38_output", 
                                 "per_barcode_metrics.csv", full.names = T, 
                                 recursive = T, include.dirs = F)
per_barcode_metrics <- c(per_barcode_metrics, per_barcode_metrics_024)
per_barcode_metrics_025 <- list.files("/nvmefs/scARC_Duan_025_GRCh38", 
                                 "per_barcode_metrics.csv", full.names = T, 
                                 recursive = T, include.dirs = F)
per_barcode_metrics <- c(per_barcode_metrics, per_barcode_metrics_025)
per_barcode_metrics_029 <- list.files("/data/FASTQ/Duan_Project_029/GRCh38_only", 
                                 "per_barcode_metrics.csv", full.names = T, 
                                 recursive = T, include.dirs = F)
per_barcode_metrics <- c(per_barcode_metrics, per_barcode_metrics_029)
per_barcode_metrics

# gex_molecule_info ####
gex_molecule_info_018 <- list.files("/nvmefs/scARC_Duan_018/GRCh38_mapped_only", 
                                      "gex_molecule_info.h5", full.names = T, 
                                      recursive = T, include.dirs = F)
gex_molecule_info <- c(gex_molecule_info_018)
gex_molecule_info_022 <- list.files("/nvmefs/scARC_Duan_022/cellranger_output", 
                                      "gex_molecule_info.h5", full.names = T, 
                                      recursive = T, include.dirs = F)
gex_molecule_info <- c(gex_molecule_info, gex_molecule_info_022)
gex_molecule_info_024 <- list.files("/nvmefs/scARC_Duan_024_raw_fastq/GRCh38_output", 
                                      "gex_molecule_info.h5", full.names = T, 
                                      recursive = T, include.dirs = F)
gex_molecule_info <- c(gex_molecule_info, gex_molecule_info_024)
gex_molecule_info_025 <- list.files("/nvmefs/scARC_Duan_025_GRCh38", 
                                      "gex_molecule_info.h5", full.names = T, 
                                      recursive = T, include.dirs = F)
gex_molecule_info <- c(gex_molecule_info, gex_molecule_info_025)
gex_molecule_info_029 <- list.files("/data/FASTQ/Duan_Project_029/GRCh38_only", 
                                      "gex_molecule_info.h5", full.names = T, 
                                      recursive = T, include.dirs = F)
gex_molecule_info <- c(gex_molecule_info, gex_molecule_info_029)
gex_molecule_info

# assemble ####
df <- data.frame(library_id = library_id,
                 atac_fragments = atac_fragments,
                 per_barcode_metrics = per_barcode_metrics,
                 gex_molecule_info = gex_molecule_info)
write.table(df, 
            file = "/nvmefs/scARC_Duan_018/018-029_combined_analysis/Duan_Project_018_to_029_samples.csv",
            quote = F, sep = ",", row.names = F, col.names = T)
