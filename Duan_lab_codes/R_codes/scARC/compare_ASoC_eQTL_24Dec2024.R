# compare ASoC datasets of 22 Dec 2023 to 4 eQTL sets
# Siwei 24 Mar 2024

# init ####
library(readr)
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
           return_val <-
             return_val$ID
           # return_val <- 
           #   return_val$ID
           return(return_val)
         })
ASoC_sig_list_unique <-
  sort(unique(unlist(ASoC_sig_list)))

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
ASoC_nonsig_list_unique <-
  sort(unique(unlist(ASoC_nonsig_list)))

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

## count intersections
enrichment_results <-
  vector(mode = "list",
         length = 4L)
names(enrichment_results) <-
  names(eQTL_list_unlisted)
  # c("sig",
  #   "nonsig")

for (i in 1:length(eQTL_list_unlisted)) {
  print(i)
  sig_enriched_count <-
    sum(ASoC_sig_list_unique %in% eQTL_list_unlisted[[i]])
  enriched_ratio <-
    sig_enriched_count / length(ASoC_sig_list_unique)
  sig_list_length <-
    length(ASoC_sig_list_unique)
  
  nonsig_enriched_count <-
    sum(sample(x = ASoC_nonsig_list_unique,
               size = length(ASoC_sig_list_unique),
               replace = F) %in% 
          eQTL_list_unlisted[[i]])
  non_enriched_ratio <-
    nonsig_enriched_count / length(ASoC_sig_list_unique)
  
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

for (i in 1:length(enrichment_results)) {
  df_2_plot <-
    as.data.frame(t(enrichment_results[[i]]))
  df_2_plot <-
    df_2_plot[, c(2, 5)]
  category_names <-
    colnames(df_2_plot)
  df_2_plot <-
    as.data.frame(t(df_2_plot))
  # df_2_plot$Category <-
  #   unlist(str_split(string = category_names,
  #                    pattern = "_ra", 
  #                    simplify = T)[, 1])
  df_2_plot$Category <-
    c("ASoC", "non-ASoC")
  # df_2_plot$
  plot_output <-
    ggplot(df_2_plot,
           aes(x = Category,
               y = V1,
               fill = Category)) +
    geom_bar(stat = "identity",
             colour = "black") +
    scale_fill_manual(values = brewer.pal(n = 4,
                                          name = "Dark2"),
                      guide = "none") +
    scale_y_continuous(expand = c(0, 0),
                       labels = scales::percent) +
    labs(y = "Percentage of enrichment") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black",
                                   size = 10)) +
    ggtitle(names(enrichment_results)[i])
  ggsave(filename  = paste0(names(enrichment_results)[i],
                            "_",
                            "enrichment",
                            '.pdf'),
         width = 2.10,
         height = 2.50,
         units = "in")
  print(plot_output)
}
# for (i in 1:length(enrichment_results)) {
#   for (j in 1:length(eQTL_list_unlisted)) {
#     print(paste(i, j))
#     
#   }
# }

# plot caQTLs ####
caQTL_files <-
  dir(path = "All_caQTLs_25kb_fdr0.05",
      pattern = ".*.csv",
      full.names = T)

caQTL_full_set <-
  vector(mode = "list",
         length = length(ASoC_files))

for (i in 1:length(caQTL_files)) {
  print(i)
  caQTL_full_set[[i]] <-
    read.table(caQTL_files[i],
               header = T,
               sep = ",")
}
names(caQTL_full_set) <-
  str_split(caQTL_files,
            pattern = "\\/",
            simplify = T)[, 2]
names(caQTL_full_set) <-
  str_split(names(caQTL_full_set),
            pattern = "\\.",
            simplify = T)[, 1]
names(caQTL_full_set)
names(ASoC_full_set)
names(caQTL_full_set) <-
  str_replace_all(string = names(caQTL_full_set),
                  pattern = "__",
                  replacement = "_")
all(names(caQTL_full_set) == names(ASoC_full_set))

for (i in 1:length(ASoC_full_set)) {
  cell_type_time <-
    names(ASoC_full_set)[i]
  print(paste(i,
              cell_type_time,
              sep = ","))
  ASoC_SNPs <-
    ASoC_full_set[[cell_type_time]]
  ASoC_SNPs <-
    ASoC_SNPs[ASoC_SNPs$FDR < 0.05, ]$ID
  # print(length(ASoC_SNPs))
  enriched_ratio <-
    sum(ASoC_SNPs %in% caQTL_full_set[[cell_type_time]]$variant_id) / length(ASoC_SNPs)
  
  non_ASoC_SNPs <-
    ASoC_full_set[[cell_type_time]]
  non_ASoC_SNPs <-
    non_ASoC_SNPs[!(non_ASoC_SNPs$FDR < 0.05), ]$ID
  non_ASoC_SNPs <-
    sample(x = non_ASoC_SNPs,
           size = length(ASoC_SNPs),
           replace = F)
  non_enriched_ratio <-
    sum(non_ASoC_SNPs %in% caQTL_full_set[[cell_type_time]]$variant_id) / length(ASoC_SNPs)
  
  df_each <-
    data.frame(V1 = c(enriched_ratio,
                      non_enriched_ratio),
               Category = c("ASoC",
                            "non-ASoC"), 
               Time = str_split(cell_type_time,
                                pattern = "_",
                                simplify = T)[, 1],
               Cell_type = str_split(cell_type_time,
                                     pattern = "_",
                                     simplify = T)[, 2])
  if (i == 1) {
    df_sum <- df_each
  } else {
    df_sum <-
      rbind(df_sum,
            df_each)
  }
  
  # plot_output <-
  #   ggplot(df_2_plot,
  #          aes(x = Category,
  #              y = V1,
  #              fill = Category)) +
  #   geom_bar(stat = "identity",
  #            colour = "black") +
  #   scale_fill_manual(values = brewer.pal(n = 4,
  #                                         name = "Dark2"),
  #                     guide = "none") +
  #   scale_y_continuous(expand = c(0, 0),
  #                      labels = scales::percent) +
  #   labs(y = "Percentage of enrichment") +
  #   theme_classic() +
  #   theme(axis.text = element_text(colour = "black",
  #                                  size = 10)) +
  #   ggtitle(paste("caQTL of",
  #                 cell_type_time))
  # ggsave(filename  = paste0(names(enrichment_results)[i],
  #                           "_",
  #                           "enrichment",
  #                           '.pdf'),
  #        width = 2.10,
  #        height = 2.50,
  #        units = "in")
  # print(plot_output)
}

df_sum$Category <-
  factor(df_sum$Category)
df_sum$Time <-
  factor(df_sum$Time)
df_sum$Cell_type <-
  factor(df_sum$Cell_type)


df_2_plot <-
  df_sum

plot_output <-
  ggplot(df_2_plot,
         aes(x = Category,
             y = V1,
             fill = Category)) +
  geom_bar(stat = "identity",
           colour = "black") +
  scale_fill_manual(values = brewer.pal(n = 4,
                                        name = "Dark2"),
                    guide = "none") +
  scale_y_continuous(expand = c(0, 0),
                     labels = scales::percent) +
  labs(y = "Percentage of enrichment") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black",
                                 size = 10)) +
  facet_grid(Cell_type ~ Time) +
  ggtitle("caQTL time x cell type")
print(plot_output)

df_2_fold <-
  data.frame(Fold = df_2_plot$V1[df_2_plot$Category == "ASoC"] /
               df_2_plot$V1[df_2_plot$Category == "non-ASoC"],
             Time = df_2_plot$Time[df_2_plot$Category == "ASoC"],
             Cell_type = df_2_plot$Cell_type[df_2_plot$Category == "ASoC"],
             stringsAsFactors = T)
ggplot(df_2_fold,
       aes(x = Cell_type,
           y = Fold,
           group = Time,
           fill = Time)) +
  geom_bar(stat = "identity",
           position = "dodge",
           colour = "black") +
  scale_fill_manual(values = brewer.pal(n = 4,
                                        name = "Set2")) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 6)) +
  labs(y = "Fold of enrichment") +
  theme_classic() +
  theme(axis.text = element_text(colour = "black",
                                 size = 10)) +
  ggtitle("caQTL time x cell type")
