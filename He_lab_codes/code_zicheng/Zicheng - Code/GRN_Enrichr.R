# Load necessary libraries
library(tidyverse)
library(enrichR)
library(openxlsx)
library(parallel)

# Set working directory
setwd("~/Documents/eGRN_Enrichr/")

# Load significant GRN data
sig_grn <- readRDS("log_Mito0.01_cellType_eGRN_corr_res_tf0.2_corr0.5.rds")

# Union targets
union_targets <- sig_grn %>% bind_rows()
union_targets <- split(union_targets$gene, union_targets$tf)

# Define databases for enrichment analysis
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", 
         "Reactome_2022", "GWAS_Catalog_2023", "DisGeNET", "MAGMA_Drugs_and_Diseases")

# Perform enrichment analysis
union_targets_enriched_list <- lapply(union_targets, function(genes) enrichr(unique(genes), dbs))
saveRDS(union_targets_enriched_list, "union_targets_enriched_list_log_Mito0.01.rds")

# Load necessary libraries for summarizing enriched terms
library(tidyverse)
library(openxlsx)

# Set working directory
setwd("~/Documents/eGRN_Enrichr/")

# Load enriched list
union_targets_enriched_list <- readRDS("union_targets_enriched_list_log_Mito0.01.rds")

# Define databases for enrichment analysis
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", 
         "Reactome_2022", "GWAS_Catalog_2023", "DisGeNET", "MAGMA_Drugs_and_Diseases")

# Aggregate enrichment results
enrichment_agg_list <- lapply(dbs, function(db) {
  db_agg <- lapply(names(union_targets_enriched_list), function(tf) {
    bind_cols(TF = tf, union_targets_enriched_list[[tf]][[db]])
  }) %>% bind_rows() %>% 
    mutate(fdr = p.adjust(P.value, method = "fdr")) %>%
    filter(fdr <= 0.1) %>% arrange(fdr)
  return(db_agg)
})
names(enrichment_agg_list) <- dbs

# Save aggregated enrichment results
write.xlsx(enrichment_agg_list, file = "Enrichr_enrichment_agg.xlsx")

# Check dimensions of enrichment results
lapply(enrichment_agg_list, dim)

# Summarize terms
term_summary <- lapply(enrichment_agg_list, function(df) {
  df %>% group_by(Term) %>% summarise(nTF = n(), min_fdr = min(fdr)) %>% arrange(desc(nTF))
})

# Save term summary
write.xlsx(term_summary, file = "Enrichr_term_summary.xlsx")