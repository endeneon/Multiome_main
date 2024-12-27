# compare ASoC datasets of 22 Dec 2023 to 4 eQTL sets
# Siwei 24 Mar 2024

# init ####
library(readr)
library(vcfR)
library(readxl)

library(ggplot2)
library(reshape2)

library(scales)

library(stringr)

library(RColorBrewer)

# load all 3x3 ASoC sets
ASoC_files <-
  dir(path = "results_100_lines_22Dec2023",
      pattern = ".*.tsv",
      full.names = T)

ASoC_full_set <-
  vector(mode = "list",
         length = length(ASoC_files))

for (i in 1:length(ASoC_files)) {
  print(i)
  ASoC_full_set[[i]] <-
    read.table(ASoC_files[i],
               header = T,
               sep = "\t")
}

ASoC_files
names(ASoC_full_set) <-
  str_split(string = ASoC_files,
            pattern = "\\/",
            simplify = T)[, 2]
names(ASoC_full_set) <-
  str_split(string = names(ASoC_full_set),
            pattern = "_SNP",
            simplify = T)[, 1]
names(ASoC_full_set)

ASoC_sig_list <-
  vector(mode = "list",
         length = length(ASoC_files))
ASoC_sig_list <-
  lapply(X = 1:9, 
         FUN = function(x) {
           return_val <- ASoC_full_set[[x]]
           return_val <-
             return_val[return_val$FDR < 0.05, ]
           # ! the Brain GTEx files do not have rsID, use chr+loc index, hg38
           return_val <-
             str_c(return_val$CHROM,
                   return_val$POS,
                   sep = "_")
           # return_val <- 
           #   return_val$ID
           return(return_val)
         })
# ASoC_sig_list_unique <-
#   sort(unique(unlist(ASoC_sig_list)))

ASoC_nonsig_list <-
  lapply(X = 1:9, 
         FUN = function(x) {
           return_val <- ASoC_full_set[[x]]
           return_val <-
             return_val[!(return_val$FDR < 0.05), ]
           return_val <-
             str_c(return_val$CHROM,
                   return_val$POS,
                   sep = "_")
           return_val <-
             sample(x = return_val,
                    size = length(ASoC_sig_list[[x]]),
                    replace = F)
           return(return_val)
         })
# ASoC_nonsig_list_unique <-
#   sort(unique(unlist(ASoC_nonsig_list)))

## load Brain Frontal Cortex dataset #####
df_brain_set <-
  read.table("Brain_Frontal_Cortex_BA9.v8.EUR.signif_pairs.txt",
             header = T, sep = "\t",
             quote = "")
df_brain_set$index <-
  str_c(str_split(df_brain_set$variant_id,
                  pattern = "_",
                  simplify = T)[, 1],
        str_split(df_brain_set$variant_id,
                  pattern = "_",
                  simplify = T)[, 2],
        sep = "_")

# count #####
## count intersections
enrichment_results <-
  vector(mode = "list",
         length = length(ASoC_sig_list))
names(enrichment_results) <-
  names(ASoC_full_set)
# c("sig",
#   "nonsig")

for (i in 1:length(ASoC_sig_list)) {
  print(i)
  sig_enriched_count <-
    sum(ASoC_sig_list[[i]] %in% df_brain_set$index)
  enriched_ratio <-
    sig_enriched_count / length(ASoC_sig_list[[i]])
  sig_list_length <-
    length(ASoC_sig_list[[i]])
  
  nonsig_enriched_count <-
    sum(ASoC_nonsig_list[[i]] %in% df_brain_set$index)
  non_enriched_ratio <-
    nonsig_enriched_count / length(ASoC_nonsig_list[[i]])
  
  chi2_test_result <-
    chisq.test(as.table(rbind(c(sig_enriched_count, nonsig_enriched_count),
                              c(sig_list_length - sig_enriched_count,
                                sig_list_length = nonsig_enriched_count))))$p.value
  
  enrichment_results[[i]] <-
    c(sig_enriched_count,
      enriched_ratio,
      sig_list_length,
      nonsig_enriched_count,
      non_enriched_ratio,
      enriched_ratio / non_enriched_ratio,
      chi2_test_result)
  names(enrichment_results[[i]]) <-
    c("sig_enriched_count",
      "enriched_ratio",
      "sig_list_length",
      "nonsig_enriched_count",
      "non_enriched_ratio",
      "enrichment_fold",
      "chi2_test_result")
}

df_sum_results <-
  as.data.frame(do.call(what = rbind,
                        args = enrichment_results))

write.table(df_sum_results,
            file = "Frontal_Cortex_ASoC_enrichment.txt",
            quote = F, sep = "\t",
            row.names = T, col.names = T)

## load Brain  Cortex dataset #####
df_brain_set <-
  read.table("Brain_Cortex.v8.EUR.signif_pairs.txt",
             header = T, sep = "\t",
             quote = "")
df_brain_set$index <-
  str_c(str_split(df_brain_set$variant_id,
                  pattern = "_",
                  simplify = T)[, 1],
        str_split(df_brain_set$variant_id,
                  pattern = "_",
                  simplify = T)[, 2],
        sep = "_")

# count #####
## count intersections
enrichment_results <-
  vector(mode = "list",
         length = length(ASoC_sig_list))
names(enrichment_results) <-
  names(ASoC_full_set)
# c("sig",
#   "nonsig")

for (i in 1:length(ASoC_sig_list)) {
  print(i)
  sig_enriched_count <-
    sum(ASoC_sig_list[[i]] %in% df_brain_set$index)
  enriched_ratio <-
    sig_enriched_count / length(ASoC_sig_list[[i]])
  sig_list_length <-
    length(ASoC_sig_list[[i]])
  
  nonsig_enriched_count <-
    sum(ASoC_nonsig_list[[i]] %in% df_brain_set$index)
  non_enriched_ratio <-
    nonsig_enriched_count / length(ASoC_nonsig_list[[i]])
  
  chi2_test_result <-
    chisq.test(as.table(rbind(c(sig_enriched_count, nonsig_enriched_count),
                              c(sig_list_length - sig_enriched_count,
                                sig_list_length = nonsig_enriched_count))))$p.value
  
  enrichment_results[[i]] <-
    c(sig_enriched_count,
      enriched_ratio,
      sig_list_length,
      nonsig_enriched_count,
      non_enriched_ratio,
      enriched_ratio / non_enriched_ratio,
      chi2_test_result)
  names(enrichment_results[[i]]) <-
    c("sig_enriched_count",
      "enriched_ratio",
      "sig_list_length",
      "nonsig_enriched_count",
      "non_enriched_ratio",
      "enrichment_fold",
      "chi2_test_result")
}

df_sum_results <-
  as.data.frame(do.call(what = rbind,
                        args = enrichment_results))
write.table(df_sum_results,
            file = "Cortex_ASoC_enrichment.txt",
            quote = F, sep = "\t",
            row.names = T, col.names = T)


## test GTEx V8 #####
## Use rsid

ASoC_sig_list <-
  vector(mode = "list",
         length = length(ASoC_files))
ASoC_sig_list <-
  lapply(X = 1:9, 
         FUN = function(x) {
           return_val <- ASoC_full_set[[x]]
           return_val <-
             return_val[return_val$FDR < 0.05, ]
           return_val <-
             return_val$ID
           # return_val <- 
           #   return_val$ID
           return(return_val)
         })
# ASoC_sig_list_unique <-
#   sort(unique(unlist(ASoC_sig_list)))

ASoC_nonsig_list <-
  lapply(X = 1:9, 
         FUN = function(x) {
           return_val <- ASoC_full_set[[x]]
           return_val <-
             return_val[!(return_val$FDR < 0.05), ]
           return_val <-
             return_val$ID
           # return_val <- 
           #   return_val$ID
           return(return_val)
         })
# ASoC_nonsig_list_unique <-
#   sort(unique(unlist(ASoC_nonsig_list)))

# count #####
## count intersections
enrichment_results <-
  vector(mode = "list",
         length = length(ASoC_sig_list))
names(enrichment_results) <-
  names(ASoC_full_set)
# c("sig",
#   "nonsig")


## load eQTL lists ####
eQTL_list_raw <-
  vector(mode = "list",
         length = 4L)
for (i in 1:length(eQTL_list_raw)) {
  eQTL_list_raw[[i]] <-
    read_excel(path = "eQTL_4_sets.xlsx",
               sheet = i)
}
names(eQTL_list_raw) <-
  excel_sheets("eQTL_4_sets.xlsx")

eQTL_list_unlisted <-
  vector(mode = "list",
         length = 4L)
for (i in 1:length(eQTL_list_raw)) {
  print(i)
  extracted_unlist <-
    unlist(eQTL_list_raw[[i]])
  extracted_unlist <-
    extracted_unlist[!is.na(extracted_unlist)]
  extracted_unlist <-
    str_split(string = extracted_unlist,
              pattern = ",",
              simplify = T)
  extracted_unlist <-
    sort(unique(unlist(extracted_unlist)))
  extracted_unlist <-
    extracted_unlist[!duplicated(extracted_unlist)]
  extracted_unlist <-
    extracted_unlist[!is.na(extracted_unlist)]
  extracted_unlist <-
    extracted_unlist[!(extracted_unlist == 'NA')]
  extracted_unlist <-
    extracted_unlist[str_length(extracted_unlist) > 0]
  eQTL_list_unlisted[[i]] <-
    extracted_unlist
}
names(eQTL_list_unlisted) <-
  excel_sheets("eQTL_4_sets.xlsx")

names(eQTL_list_unlisted)

df_brain_set <-
  eQTL_list_unlisted[['GTEx_V8']]

for (i in 1:length(ASoC_sig_list)) {
  print(i)
  sig_enriched_count <-
    sum(ASoC_sig_list[[i]] %in% df_brain_set)
  enriched_ratio <-
    sig_enriched_count / length(ASoC_sig_list[[i]])
  sig_list_length <-
    length(ASoC_sig_list[[i]])
  
  # nonsig_enriched_count <-
  #   sum(sample(x = dbSNP_index,
  #              size = length(ASoC_sig_list[[i]]),
  #              replace = F) %in%
  #         df_brain_set)
  nonsig_enriched_count <-
    sum(sample(x = ASoC_nonsig_list[[i]],
               size = length(ASoC_sig_list[[i]]),
               replace = F) %in% 
          df_brain_set)
  non_enriched_ratio <-
    nonsig_enriched_count / length(ASoC_sig_list[[i]])
  
  chi2_test_result <-
    chisq.test(as.table(rbind(c(sig_enriched_count, nonsig_enriched_count),
                              c(sig_list_length - sig_enriched_count,
                                sig_list_length = nonsig_enriched_count))))$p.value
  
  enrichment_results[[i]] <-
    c(sig_enriched_count,
      enriched_ratio,
      sig_list_length,
      nonsig_enriched_count,
      non_enriched_ratio,
      enriched_ratio / non_enriched_ratio,
      chi2_test_result)
  names(enrichment_results[[i]]) <-
    c("sig_enriched_count",
      "enriched_ratio",
      "sig_list_length",
      "nonsig_enriched_count",
      "non_enriched_ratio",
      "enrichment_fold",
      "chi2_test_result")
}

df_sum_results <-
  as.data.frame(do.call(what = rbind,
                        args = enrichment_results))

write.table(df_sum_results,
            file = "GTEx_v8_ASoC_enrichment.txt",
            quote = F, sep = "\t",
            row.names = T, col.names = T)


## test commonmind
names(eQTL_list_unlisted)
df_brain_set <-
  eQTL_list_unlisted[['CommonMind']]

for (i in 1:length(ASoC_sig_list)) {
  print(i)
  sig_enriched_count <-
    sum(ASoC_sig_list[[i]] %in% df_brain_set)
  enriched_ratio <-
    sig_enriched_count / length(ASoC_sig_list[[i]])
  sig_list_length <-
    length(ASoC_sig_list[[i]])
  
  # nonsig_enriched_count <-
  #   sum(sample(x = dbSNP_index,
  #              size = length(ASoC_sig_list[[i]]),
  #              replace = F) %in%
  #         df_brain_set)
  nonsig_enriched_count <-
    sum(sample(x = ASoC_nonsig_list[[i]],
               size = length(ASoC_sig_list[[i]]),
               replace = F) %in% 
          df_brain_set)
  non_enriched_ratio <-
    nonsig_enriched_count / length(ASoC_sig_list[[i]])
  
  chi2_test_result <-
    chisq.test(as.table(rbind(c(sig_enriched_count, nonsig_enriched_count),
                              c(sig_list_length - sig_enriched_count,
                                sig_list_length = nonsig_enriched_count))))$p.value
  
  enrichment_results[[i]] <-
    c(sig_enriched_count,
      enriched_ratio,
      sig_list_length,
      nonsig_enriched_count,
      non_enriched_ratio,
      enriched_ratio / non_enriched_ratio,
      chi2_test_result)
  names(enrichment_results[[i]]) <-
    c("sig_enriched_count",
      "enriched_ratio",
      "sig_list_length",
      "nonsig_enriched_count",
      "non_enriched_ratio",
      "enrichment_fold",
      "chi2_test_result")
}

df_sum_results <-
  as.data.frame(do.call(what = rbind,
                        args = enrichment_results))

write.table(df_sum_results,
            file = "CommonMind_ASoC_enrichment.txt",
            quote = F, sep = "\t",
            row.names = T, col.names = T)
