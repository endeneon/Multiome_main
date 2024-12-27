########GO enrichment of GRN targets
library(tidyverse)
library(enrichR)
library(openxlsx)
library(parallel)


setwd("~/Documents/eGRN_Enrichr/")

sig_grn = readRDS("log_Mito0.01_cellType_eGRN_corr_res_tf0.2_corr0.5.rds")

######### union_targets
union_targets = sig_grn %>% bind_rows(.)
union_targets = split(union_targets$gene,union_targets$tf)

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023","Reactome_2022","GWAS_Catalog_2023","DisGeNET","MAGMA_Drugs_and_Diseases")

union_targets_enriched_list <- lapply(union_targets, function(x) enrichr(unique(x), dbs))
saveRDS(union_targets_enriched_list,"union_targets_enriched_list_log_Mito0.01.rds")


############Summarize enriched terms
library(tidyverse)
library(openxlsx)

setwd("~/Documents/eGRN_Enrichr/")

union_targets_enriched_list = readRDS("union_targets_enriched_list_log_Mito0.01.rds")


dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023","Reactome_2022","GWAS_Catalog_2023","DisGeNET","MAGMA_Drugs_and_Diseases")

enrichment_agg_list = lapply(dbs, function(dat){
  dat_agg = lapply(names(union_targets_enriched_list), function(tf){
    bind_cols(TF = tf,union_targets_enriched_list[[tf]][[dat]])
  }) %>% bind_rows() %>% mutate(fdr = p.adjust(P.value, method = "fdr")) %>%
    filter(fdr <= 0.1) %>% arrange(fdr)
})
names(enrichment_agg_list) = dbs
write.xlsx(enrichment_agg_list, file = "Enrichr_enrichment_agg.xlsx")

lapply(enrichment_agg_list,dim)

term_summary = lapply(enrichment_agg_list, function(x){
  x %>% group_by(Term) %>% summarise(nTF = n(), min_fdr = min(fdr)) %>% arrange(desc(nTF))
})

## save term_summary using openxlsx
write.xlsx(term_summary, file = "Enrichr_term_summary.xlsx")