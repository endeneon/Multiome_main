library(tidyverse)
library(parallel)

conditions <- paste0(c("0hr", "1hr", "6hr"), "_", unlist(lapply(c("GABA", "nmglut", "npglut"), function(x) rep(x, 3))))

join_qtl_asoc <- function(condition) {
  asoc_file <- list.files("/project/xinhe/zicheng/neuron_stimulation/data/Siwe_ASOC/results_100_lines_22Dec2023", pattern = paste0(condition, ".*"), full.names = TRUE)[1]
  caQTL_file <- list.files("/project/xinhe/zicheng/neuron_stimulation/caQTL/output/caQTL_mapping", pattern = paste0("Cis_Nominal_", sub("_", "__", condition), ".*"), full.names = TRUE)[1]
  
  asoc_res <- read_tsv(asoc_file)
  caQTL_res <- readRDS(caQTL_file) %>%
    filter(variant_id %in% unique(asoc_res$ID), start_distance >= -251 & start_distance < 250)
  
  joined_res <- inner_join(asoc_res, caQTL_res, by = c("ID" = "variant_id"))
  
  return(joined_res)
}

join_qtl_asoc_list <- mclapply(conditions, join_qtl_asoc, mc.cores = 3)
names(join_qtl_asoc_list) <- conditions

join_qtl_asoc_agg <- lapply(names(join_qtl_asoc_list), function(x) {
  join_qtl_asoc_list[[x]] %>% mutate(condition = x)
}) %>% bind_rows()

### Filter ASoC with FDR <= 0.05
ascertain_asoc_joined <- join_qtl_asoc_agg %>% filter(FDR <= 0.05)

saveRDS(ascertain_asoc_joined, "/project/xinhe/zicheng/neuron_stimulation/caQTL/output/dynamic_QTL/asoc_filtered_caqtl.rds")

library(tidyverse)
library(vroom)

snp_pos <- vroom("/project2/xinhe/lifanl/neuron_stim/geno_info_100lines.txt", col_names = c("chr", "pos", "a2", "a1", "rsid")) %>%
  as.data.frame()
rownames(snp_pos) <- snp_pos$rsid

### Check alternative allele and reference allele
ascertain_asoc_joined <- readRDS("/project/xinhe/zicheng/neuron_stimulation/caQTL/output/dynamic_QTL/asoc_filtered_caqtl.rds") %>%
  mutate(caQTL_ALT = snp_pos[ID, "a1"], caQTL_REF = snp_pos[ID, "a2"]) %>%
  filter(caQTL_ALT == ALT, caQTL_REF == REF)

ascertain_asoc_joined <- ascertain_asoc_joined %>% mutate(ALT_ratio = 1 - REF_ratio)
lm(slope ~ ALT_ratio, ascertain_asoc_joined) %>% summary()

ascertain_asoc_joined %>% ggplot(aes(x = (1 - REF_ratio), y = slope)) +
  geom_point(size = 0.5) + 
  labs(x = "ASoC Allele Fraction", y = "caQTL Effect Size") + theme_bw() +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  geom_abline(intercept = -1.46489, slope = 2.77986, color = "red")

ggsave("/project/xinhe/zicheng/neuron_stimulation/output/final_figures/fig_7c.pdf", height = 4, width = 4)

cor(ascertain_asoc_joined$REF_ratio, ascertain_asoc_joined$slope) * -1



