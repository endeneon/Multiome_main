# recalc p values


# init ####
{
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(stringr)
  
  library(RColorBrewer)
  
  library(ggpubr)
  
  library(sva)
  library(harmony)
  
  library(reshape2)
  
  library(readr)
  library(readxl)
}

input_list <-
  list.files(path = "018-029_combat_limma_agesq_normalize_all_times_sep_DEG_results_csv",
             pattern = "*_all_genes.csv",
             full.names = T)

input_list <-
  list.files(path = "018-029_diff_response_filter_by_cpm_casevcontrol_lm_results/unfiltered_full_results/",
             pattern = "*_all_genes.csv",
             full.names = T)
sample_names <-
  str_split(string = input_list,
            pattern = "_sep_",
            simplify = T)[, 3]
sample_names <-
  str_split(string = sample_names,
            pattern = "_all_",
            simplify = T)[, 1]

full_data_set <-
  vector(mode = "list",
         length = length(sample_names))

full_data_set <-
  lapply(X = input_list, 
         FUN = function(x) {
           read_data <-
             read_csv(x)
           return(read_data)
         })
names(full_data_set) <-
  sample_names

# recalc Padj

full_data_set <-
  lapply(X = full_data_set, 
         FUN = function(x) {
           x$new_P_adj <-
             p.adjust(x$P.Value,
                      method = "fdr")
           return(x)
         })
