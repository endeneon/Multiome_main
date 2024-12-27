# Chuxuan Li 02/02/2023
# Generate aggregation csv for 10x arc aggr

library(stringr)
library_id <- str_extract(list.dirs("/data/FASTQ/Duan_Project_029/GRCh38_only", 
                                    recursive = F), "[0-9]+-[0|1|6]")
atac_fragments <- list.files("/data/FASTQ/Duan_Project_029/GRCh38_only", 
                             "atac_fragments.tsv.gz$", full.names = T, 
                             recursive = T, include.dirs = F)
per_barcode_metrics <- list.files("/data/FASTQ/Duan_Project_029/GRCh38_only", 
                                  "per_barcode_metrics.csv", full.names = T, 
                                  recursive = T, include.dirs = F)
gex_molecule_info <- list.files("/data/FASTQ/Duan_Project_029/GRCh38_only", 
                                  "gex_molecule_info.h5", full.names = T, 
                                  recursive = T, include.dirs = F)
df <- data.frame(library_id = library_id,
                 atac_fragments = atac_fragments,
                 per_barcode_metrics = per_barcode_metrics,
                 gex_molecule_info = gex_molecule_info)
write.table(df, 
          file = "/data/FASTQ/Duan_Project_029/GRCh38_only/Duan_Project_029_samples.csv",
          quote = F, sep = ",", row.names = F, col.names = T)
