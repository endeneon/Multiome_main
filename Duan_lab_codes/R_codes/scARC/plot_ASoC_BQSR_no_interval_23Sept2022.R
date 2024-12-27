# Siwei 23 Sept 2022
# use No intervals
# all other conditions same to 15 Aug 
# calculate ASoC ratio from the aggregation
# of down sampled results

# init
library(ggplot2)
library(readr)
library(vcfR)

library(parallel)
library(stringr)

library(MASS)

library(RColorBrewer)
library(grDevices)

# input_file_list <-
#   list.files(path = "BQSR_downsampled_w_global_interval_26Sept2022",
#              pattern = "\\.txt$",
#              full.names = T)
# output_path <- "BQSR_downsampled_w_global_interval_26Sept2022/output"

# input_file_list <-
#   list.files(path = "BQSR_downsampled_no_interval_25Sept2022",
#              pattern = "\\.txt$",
#              full.names = T)
# output_path <- "BQSR_downsampled_no_interval_25Sept2022/output"

input_file_list <-
  list.files(path = "known_genotype_only",
             pattern = "all_mgs_interval_BQSR_28Dec2022_4_R\\.txt$",
             full.names = T)
output_path <- "known_genotype_only/output"

input_file_list <-
  list.files(path = "known_genotype_only",
             pattern = "16Dec2022_4_R\\.txt$",
             full.names = T)
output_path <- "known_genotype_only/output"

input_file_list <-
  list.files(path = "NotDuplicateReadFilter_15Aug2022",
             pattern = "NEFM_neg_0hr.*\\.txt$",
             full.names = T)
output_path <- "NotDuplicateReadFilter_15Aug2022/output"


# input_file_list <-
#   list.files(path = "same_size_w_upsampling", 
#              pattern = "12Sept2022\\.txt$",
#              full.names = T)
# output_path <- "same_size_w_upsampling/output"

# date_prefix <-
#   format(Sys.Date(), "%b_%d_%Y_")
# date_prefix <- "Sep_25_2022_"

input_file_list <-
  list.files(path = "BQSR_no_interval_14Oct2022",
             pattern = "\\.txt$",
             full.names = T)
output_path <- "BQSR_no_interval_14Oct2022/output/"


input_file_list <-
  list.files(path = "BQSR_no_interval_23Sept2022",
             pattern = "NEFM_neg_0hr.*\\.txt$",
             full.names = T)
output_path <- "BQSR_no_interval_23Sept2022/output"

input_file_list <-
  list.files(path = "known_genotype_only",
             pattern = "\\.txt$",
             full.names = T)
output_path <- "known_genotype_only/output"

input_file_list <-
  list.files(path = "NotDuplicateReadFilter_15Aug2022",
             pattern = "NEFM_neg_.*\\.txt$",
             full.names = T)
output_path <- "NotDuplicateReadFilter_15Aug2022/output"

input_file_list <-
  list.files(path = "BQSR_no_interval_23Sept2022",
             pattern = "NEFM_neg_.*\\.txt$",
             full.names = T)
output_path <- "BQSR_no_interval_23Sept2022/output"

input_file_list <-
  list.files(path = "masked_genome_05Dec2022",
             pattern = "NEFM_neg_.*\\.txt$",
             full.names = T)
output_path <- "masked_genome_05Dec2022/output"
dir.create(output_path)

i <- 1L
for (i in 1:length(input_file_list)) {
  df_raw <- read_delim(input_file_list[i],
                       delim = "\t", escape_double = FALSE,
                       trim_ws = TRUE)
  
  output_file_full_name <-
    str_split(string = input_file_list[i],
              pattern = "\\/",
              simplify = T)[2]
  
  # output_file_full_name <-
  #   str_remove_all(string = output_file_full_name,
  #                  pattern = '_scATAC_Dedup')

  output_file_full_name <-
    str_remove_all(string = output_file_full_name,
                   pattern = 'all_mgs_interval_BQSR_28Dec2022_4_R')
  # output_file_full_name <-
  #   str_remove_all(string = output_file_full_name,
  #                  pattern = '_scATAC_Dedup_08Nov2022_4_R')
  
  # output_file_full_name <-
  #   str_remove_all(string = output_file_full_name,
  #                  pattern = '_scATAC_BQSR_16Dec2022_4_R')
  output_file_full_name_short <- 
    str_replace_all(string = output_file_full_name,
                    pattern = "\\.txt",
                    replacement = "")
    # output_file_full_name
  # output_file_full_name <-
  #   str_remove_all(string = output_file_full_name,
  #                  pattern = '_scATAC_NotDuplicatedReadFilter_15Aug2022_4_R')
  # output_file_full_name <-
  #   str_replace_all(string = output_file_full_name,
  #                   pattern = "\\.txt",
  #                   replacement = "_15Aug2022_NCALLED_less_than_14\\.png")
  
  # output_file_full_name <-
  #   str_remove_all(string = output_file_full_name,
  #                  pattern = '_scATAC_Dedup_23Sept2022_4_R')
  # output_file_full_name_short <- output_file_full_name
  # output_file_full_name <-
  #   str_replace_all(string = output_file_full_name,
  #                   pattern = "\\.txt",
  #                   replacement = "_23Sept2022_NCALLED_less_than_14\\.png")
  
  # pattern = '_downsampled_median_downsampled_26Sept2022'
  
  
  # output_file_full_name <-
  #   str_replace_all(string = output_file_full_name,
  #                   pattern = "\\.txt",
  #                   replacement = "\\.png")
  # 
  # output_file_full_name <-
  #   paste(output_path,
  #         output_file_full_name,
  #         sep = "/")
  
  
  print(output_file_full_name)
  # print(str_replace_all(string = output_file_full_name,
  #                       pattern = "\\.txt",
  #                       replacement = "_15Aug2022_NCALLED_less_than_14\\.png"))
  
  df_to_plot <- df_raw
  df_to_plot <- df_to_plot[df_to_plot$ID != ".", ]
  
  # remove all SNPs with NCALLED >=14
  # df_to_plot <-
  #   df_to_plot[df_to_plot$NCALLED < 14, ]
  
  # make a 8 core cluster
  clust_df_plot_ref <- makeCluster(type = "FORK", 4)
  clusterExport(clust_df_plot_ref, "df_to_plot")
  
  df_to_plot$REF_C <- 
    parApply(cl = clust_df_plot_ref, X = df_to_plot, MARGIN = 1, 
             FUN = function(x)(sum(as.numeric(x[seq(from = 7,
                                                    to = ncol(df_to_plot) - 1,
                                                    by = 2)]))))
  stopCluster(clust_df_plot_ref)
  rm(clust_df_plot_ref)
  print("-")
  
  clust_df_plot_alt <- makeCluster(type = "FORK", 4)
  clusterExport(clust_df_plot_alt, "df_to_plot")
  
  df_to_plot$ALT_C <- 
    parApply(cl = clust_df_plot_alt, X = df_to_plot, MARGIN = 1,
             FUN = function(x)(sum(as.numeric(x[seq(from = 8,
                                                    to = ncol(df_to_plot),
                                                    by = 2)]))))
  stopCluster(clust_df_plot_alt)
  rm(clust_df_plot_alt)
  print("-")
  
  df_to_plot <- df_to_plot[df_to_plot$REF_C > df_to_plot$NCALLED, ]
  df_to_plot <- df_to_plot[df_to_plot$ALT_C > df_to_plot$NCALLED, ]
  df_to_plot$DP <- df_to_plot$REF_C + df_to_plot$ALT_C
  print(sum(df_to_plot$DP))
  
  # sum(df_to_plot$DP > 19) 
  
  df_to_plot <- df_to_plot[df_to_plot$DP > 19, ]
  # df_to_plot <- df_to_plot[df_to_plot$DP < 1000, ]
  
  
  
  
  clust_df_plot_use <- makeCluster(type = "FORK", 4)
  clusterExport(clust_df_plot_use, "df_to_plot")
  
  df_to_plot$pVal <- 
    parApply(cl = clust_df_plot_use, X = df_to_plot, MARGIN = 1,
             FUN = function(x)(binom.test(x = as.numeric(x[ncol(df_to_plot) - 1]),
                                          n = sum(as.numeric(x[ncol(df_to_plot) - 1]), 
                                                  as.numeric(x[ncol(df_to_plot) - 2])),
                                          p = 0.5,
                                          alternative = "t")$p.value))
  stopCluster(clust_df_plot_use)
  rm(clust_df_plot_use)
  print("-")
  
  
  df_to_plot$FDR <- p.adjust(p = df_to_plot$pVal,
                             method = "fdr")
  df_to_plot_output <- df_to_plot[, c(1:6, (ncol(df_to_plot) - 4):ncol(df_to_plot))]
  df_to_plot_output <- df_to_plot_output[order(df_to_plot_output$pVal), ]
  ### cut only FDR < 0.05 for annotation
  # df_to_plot_output <- df_to_plot_output[df_to_plot_output$FDR < 0.05, ]
  
  # sum(df_to_plot$FDR < 0.05)
  
  ### make PNG plot and send out
  ## !! assign time and cell type !!
  png(filename = paste(output_path,
                       str_replace_all(string = output_file_full_name,
                                       pattern = "\\.txt",
                                       replacement = "\\.png"),
                       sep = "/"),
      width = 400, height = 500, 
      bg = "white")
  
  fig_to_plot <-
    ggplot(df_to_plot,
           aes(x = (REF_C/DP),
               y = (0 - log10(pVal)))) +
    geom_point(size = 0.5,
               aes(colour = factor(ifelse(FDR < 0.05, 
                                          "FDR < 0.05", 
                                          "FDR > 0.05")))) +
    scale_color_manual(values = c("red", "black")) +
    labs(colour = "FDR value") +
    ggtitle(paste(#"",
                  output_file_full_name_short,
                  " DP >= 20, 18 lines\n",
                  # "SNP lists from 23 Sept results, removed all NCALLED >= 14\n",
                  nrow(df_to_plot),
                  "SNPs, of",
                  sum(df_to_plot$FDR < 0.05),
                  "FDR < 0.05\n",
                  "FDR < 0.05 & REF/DP < 0.5 = ", 
                  sum((df_to_plot$REF_C/df_to_plot$DP < 0.5) & (df_to_plot$FDR < 0.05)),
                  ";\n FDR < 0.05 & REF/DP > 0.5 = ", "",
                  sum((df_to_plot$REF_C/df_to_plot$DP > 0.5) & (df_to_plot$FDR < 0.05)), "\n",
                  "rs2027349: REF_C =",
                  df_to_plot[df_to_plot$ID %in% "rs2027349", ]$REF_C,
                  ", DP = ",
                  df_to_plot[df_to_plot$ID %in% "rs2027349", ]$DP, 
                  ",\n Pval =",
                  format(signif(df_to_plot[df_to_plot$ID %in% "rs2027349", ]$pVal, digits = 2),
                         nsmall = 2), 
                  ", FDR =",
                  format(signif(df_to_plot[df_to_plot$ID %in% "rs2027349", ]$FDR, digits = 2),
                         nsmall = 2))) +
    ylab("-log10P") +
    ylim(0, 50) +
    theme_classic()
  # fig_to_plot <-
  #   ggplot(df_to_plot,
  #          aes(x = (REF_C/DP),
  #              y = (0 - log10(pVal)))) +
  #   geom_point(size = 0.5,
  #              aes(colour = factor(ifelse(FDR < 0.05, 
  #                                         "FDR < 0.05", 
  #                                         "FDR > 0.05")))) +
  #   scale_color_manual(values = c("red", "black")) +
  #   labs(colour = "FDR value") +
  #   ggtitle(paste(output_file_full_name, "\tDP >= 20, 18 lines\n",
  #                 nrow(df_to_plot),
  #                 "SNPs, of",
  #                 sum(df_to_plot$FDR < 0.05),
  #                 "FDR < 0.05\n",
  #                 "FDR < 0.05 & REF_C/DP < 0.5 = ", 
  #                 sum((df_to_plot$REF_C/df_to_plot$DP < 0.5) & (df_to_plot$FDR < 0.05)), ";\n",
  #                 "FDR < 0.05 & REF_C/DP > 0.5 = ",
  #                 sum((df_to_plot$REF_C/df_to_plot$DP > 0.5) & (df_to_plot$FDR < 0.05)))) +
  #   ylab("-log10P") +
  #   theme_classic() +
  #   theme(legend.position = "none")
  
  # fig_to_plot <-
  #   ggplot(df_to_plot,
  #          aes(x = (REF_C/DP),
  #              y = (0 - log10(pVal)))) +
  #   geom_point(size = 0.5,
  #              aes(colour = factor(ifelse(FDR < 0.05, 
  #                                         "FDR < 0.05", 
  #                                         "FDR > 0.05")))) +
  #   scale_color_manual(values = c("red", "black")) +
  #   labs(colour = "FDR value") +
  #   ggtitle(paste(output_file_full_name,"\n", "DP >= 20, 18 lines\n",
  #                 # "Duplicates from picard MarkDuplicates\n",
  #                 # "Duplicates marked by 10x Cell Ranger\n",
  #                 nrow(df_to_plot),
  #                 "SNPs, of",
  #                 sum(df_to_plot$FDR < 0.05),
  #                 "FDR < 0.05",
  #                 )) +
  #   ylab("-log10P") +
  #   theme_classic()
  
  print(fig_to_plot)
  
  dev.off()
  ###
  
  print(df_to_plot_output[df_to_plot_output$ID %in% "rs2027349", ])
  

  
  df_to_annot <- df_to_plot_output[, c(1, 2, 2, 4, 5, 3, 6:11)]
  
  write.table(df_to_annot,
              file = 
                str_replace_all(string = output_file_full_name,
                                pattern = "\\.txt",
                                replacement = "\\.avinput"),
              quote = F, sep = "\t",
              row.names = F, col.names = F)
  
}

sum((df_to_plot_output$REF_C/df_to_plot_output$DP) < 0.5)
sum((df_to_plot_output$REF_C/df_to_plot_output$DP) > 0.5)

date_prefix <-
  format(Sys.Date(), "%b_%d_%Y_")


# View(df_to_plot[df_to_plot$ID %in% "rs2027349", ])

## !! assign time and cell type !!
ggplot(df_to_plot,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle(paste(output_file_full_name, "DP >= 20, 18 lines\n",
                nrow(df_to_plot),
                "SNPs, of",
                sum(df_to_plot$FDR < 0.05),
                "FDR < 0.05\n",
                "FDR < 0.05 & REF/DP < 0.5 = ", 
                sum((df_to_plot$REF_C/df_to_plot$DP < 0.5) & (df_to_plot$FDR < 0.05)),
                ";\n FDR < 0.05 & REF/DP > 0.5 = ", "",
                sum((df_to_plot$REF_C/df_to_plot$DP > 0.5) & (df_to_plot$FDR < 0.05)), "\n",
                "rs2027349: REF_C =",
                df_to_plot[df_to_plot$ID %in% "rs2027349", ]$REF_C,
                ", DP = ",
                df_to_plot[df_to_plot$ID %in% "rs2027349", ]$DP, 
                ",\n Pval =",
                format(signif(df_to_plot[df_to_plot$ID %in% "rs2027349", ]$pVal, digits = 2),
                       nsmall = 2), 
                ", FDR =",
                format(signif(df_to_plot[df_to_plot$ID %in% "rs2027349", ]$FDR, digits = 2),
                       nsmall = 2))) +
  ylab("-log10P") +
  ylim(0, 40) +
  theme_classic()

ggplot(df_to_plot,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle(paste("",
                # output_file_full_name, 
                "DP >= 20, 18 lines, 08 Oct 2022 model\n",
                nrow(df_to_plot),
                "SNPs, of",
                sum(df_to_plot$FDR < 0.05),
                "FDR < 0.05\n",
                "FDR < 0.05 & REF_C/DP < 0.5 = ", 
                sum((df_to_plot$REF_C/df_to_plot$DP < 0.5) & (df_to_plot$FDR < 0.05)), ";\n",
                "FDR < 0.05 & REF_C/DP > 0.5 = ",
                sum((df_to_plot$REF_C/df_to_plot$DP > 0.5) & (df_to_plot$FDR < 0.05)))) +
  ylab("-log10P") +
  theme_classic() +
  theme(legend.position = "none")



# ggplot(df_to_plot,
#        aes(x = (REF_C/DP),
#            y = (0 - log10(pVal)))) +
#   geom_point(size = 0.5,
#              aes(colour = factor(ifelse(FDR < 0.05, 
#                                         "FDR < 0.05", 
#                                         "FDR > 0.05")))) +
#   scale_color_manual(values = c("red", "black")) +
#   labs(colour = "FDR value") +
#   ggtitle(paste("NEFM pos glut, 0hr, DP >= 20, 18 lines\n",
#                 # "Duplicates from picard MarkDuplicates\n",
#                 "Duplicates marked by 10x Cell Ranger\n",
#                 nrow(df_to_plot),
#                 "SNPs, of",
#                 sum(df_to_plot$FDR < 0.05),
#                 "FDR < 0.05")) +
#   ylab("-log10P") +
#   theme_classic()



sum(df_to_plot$REF_C/df_to_plot$DP < 0.5) #
sum(df_to_plot$REF_C/df_to_plot$DP > 0.5) #

sum((df_to_plot$REF_C/df_to_plot$DP < 0.5) & (df_to_plot$FDR < 0.05)) #29685

sum((df_to_plot$REF_C/df_to_plot$DP > 0.5) & (df_to_plot$FDR < 0.05)) #29685


sum(df_to_plot$REF_C/df_to_plot$DP > 0.5) #26342

df_to_plot[df_to_plot$ID %in% "rs2027349", 42:47]


df_test <- df_to_plot[, c(1,2,3)]
df_test <- df_test[, c(1,2,3,3)]


# 22484 of 94960

df_no_upsampling <- df_to_plot_output



##### use vcfR to assemble table and remove all homozygous with DP < 2

vcf_file_list <-
  list.files(path = "~/subsetted_bams/downsample_read_count_vcf/results_20Sept2022", 
             pattern = ".*vcf$",
             full.names = T)

make_violin_plot <-
  function(input_vcf) {
    vcf_test <-
      read.vcfR(file = input_vcf,
                limit = 1e+8)
    cell_type_time <-
      basename(input_vcf)
    cell_type_time <-
      str_remove_all(string = cell_type_time,
                     pattern = "_downsampled_median\\.vcf")
    
    
    # assemble table, get GT
    df_gt <-
      as.data.frame(extract.gt(vcf_test, 
                               element = c("GT"), 
                               return.alleles = F))
    # get DP
    df_DP <-
      as.data.frame(extract.gt(vcf_test, 
                               element = c("DP"), 
                               return.alleles = F))
    # get AD
    df_AD <-
      as.data.frame(extract.gt(vcf_test, 
                               element = c("AD"), 
                               return.alleles = F))
    df_assembled_list <-
      vector(mode = "list",
             length = 6L)
    names(df_assembled_list) <-
      c("GT", "DP", "AD", "REF_N", "ALT_N", "return_sum")
    
    df_assembled_list[["GT"]] <- df_gt
    df_assembled_list[["DP"]] <- df_DP
    df_assembled_list[["AD"]] <- df_AD
    
    df_assembled_list[["AD"]][(df_assembled_list[["DP"]] < 3) & 
                                ((df_assembled_list[["GT"]] == "1/1") | 
                                   (df_assembled_list[["GT"]] == "0/0"))] <-
      "0,0"
    
    # make a 16 core cluster
    clust_df_plot_ref <- makeCluster(type = "FORK", 16)
    df_count_temp <- df_assembled_list[["AD"]]
    clusterExport(clust_df_plot_ref, "df_count_temp")
    
    df_assembled_list[["REF_N"]] <-
      parApply(cl = clust_df_plot_ref, 
               X = df_count_temp, MARGIN = 1, 
               FUN = function(x)(str_split(x,
                                           pattern = ",",
                                           simplify = T)[, 1]))
    stopCluster(clust_df_plot_ref)
    rm(clust_df_plot_ref)
    df_assembled_list[["REF_N"]] <-
      as.numeric(df_assembled_list[["REF_N"]])
    
    df_assembled_list[["REF_N"]] <-
      matrix(unlist(df_assembled_list[["REF_N"]]),
             nrow = nrow(df_DP),
             # ncol = 18,
             byrow = T)
    dimnames(df_assembled_list[["REF_N"]]) <-
      dimnames(df_assembled_list[["GT"]])
    
    
    # make a 16 core cluster
    clust_df_plot_ref <- makeCluster(type = "FORK", 16)
    df_count_temp <- df_assembled_list[["AD"]]
    clusterExport(clust_df_plot_ref, "df_count_temp")
    
    df_assembled_list[["ALT_N"]] <-
      parApply(cl = clust_df_plot_ref, 
               X = df_count_temp, MARGIN = 1, 
               FUN = function(x)(str_split(x,
                                           pattern = ",",
                                           simplify = T)[, 2]))
    stopCluster(clust_df_plot_ref)
    rm(clust_df_plot_ref)
    df_assembled_list[["ALT_N"]] <-
      as.numeric(df_assembled_list[["ALT_N"]])
    df_assembled_list[["ALT_N"]] <-
      matrix(unlist(df_assembled_list[["ALT_N"]]),
             nrow = nrow(df_DP),
             # ncol = 18,
             byrow = T)
    dimnames(df_assembled_list[["ALT_N"]]) <-
      dimnames(df_assembled_list[["GT"]])
    
    
    df_assembled_list[["return_sum"]]$REF_N <-
      rowSums(df_assembled_list[["REF_N"]])
    df_assembled_list[["return_sum"]]$ALT_N <-
      rowSums(df_assembled_list[["ALT_N"]])
    df_assembled_list[["return_sum"]]$rsID <-
      rownames(df_assembled_list[["GT"]])
    df_assembled_list[["return_sum"]]$DP <-
      df_assembled_list[["return_sum"]]$REF_N +
      df_assembled_list[["return_sum"]]$ALT_N 
    
    df_assembled_list[["return_sum"]] <-
      as.data.frame(df_assembled_list[["return_sum"]],
                    stringsAsFactors = F)
    
    # minDP = 30
    df_assembled_list[["return_sum"]] <-
      df_assembled_list[["return_sum"]][df_assembled_list[["return_sum"]]$DP > 29, ]
    df_assembled_list[["return_sum"]] <-
      df_assembled_list[["return_sum"]][df_assembled_list[["return_sum"]]$DP < 1500, ]
    df_assembled_list[["return_sum"]] <-
      df_assembled_list[["return_sum"]][df_assembled_list[["return_sum"]]$REF_N > 19, ]
    df_assembled_list[["return_sum"]] <-
      df_assembled_list[["return_sum"]][df_assembled_list[["return_sum"]]$ALT_N > 19, ]
    
    clust_df_plot_use <- makeCluster(type = "FORK", 16)
    df_count_temp <- df_assembled_list[["return_sum"]] 
    clusterExport(clust_df_plot_use, "df_count_temp")
    
    df_assembled_list[["return_sum"]]$pVal <- 
      parApply(cl = clust_df_plot_use, 
               X = df_count_temp, MARGIN = 1,
               FUN = function(x)(binom.test(x = as.numeric(x[1]),
                                            n = as.numeric(x[4]),
                                            p = 0.5,
                                            alternative = "t")$p.value))
    stopCluster(clust_df_plot_use)
    rm(clust_df_plot_use)
    
    df_assembled_list[["return_sum"]]$FDR <-
      p.adjust(p = df_assembled_list[["return_sum"]]$pVal,
               method = "fdr")
    
    
    df_to_plot <-
      df_assembled_list[["return_sum"]]
    
    png(filename = paste(output_path,
                         cell_type_time,
                         '.png',
                         sep = ""),
        width = 400, height = 400, 
        bg = "white")
    
    fig_to_plot <-
      ggplot(df_to_plot,
             aes(x = (REF_N/DP),
                 y = (0 - log10(pVal)))) +
      geom_point(size = 0.5,
                 aes(colour = factor(ifelse(FDR < 0.05, 
                                            "FDR < 0.05", 
                                            "FDR > 0.05")))) +
      scale_color_manual(values = c("red", "black")) +
      labs(colour = "FDR value") +
      ggtitle(paste(cell_type_time,"\n", "DP = [30, 1500], 18 lines\n",
                    # "Duplicates from picard MarkDuplicates\n",
                    # "Duplicates marked by 10x Cell Ranger\n",
                    nrow(df_to_plot),
                    "SNPs, of",
                    sum(df_to_plot$FDR < 0.05),
                    "FDR < 0.05")) +
      ylab("-log10P") +
      theme_classic()
    
    print(fig_to_plot)
    
    dev.off()
    
    return(df_to_plot)
  }

df_results_type_list <-
  vector(mode = "list", length = 9L)
for (i in 1:9) {
  df_results_type_list[[i]] <-
    make_violin_plot(input_vcf = vcf_file_list[i])
}





ggplot(df_to_plot,
       aes(x = (REF_N/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle(paste(cell_type_time, ", DP = [30, 1500] 18 lines\n",
                # "Duplicates from picard MarkDuplicates\n",
                # "Duplicates marked by 10x Cell Ranger\n",
                nrow(df_to_plot),
                "SNPs, of",
                sum(df_to_plot$FDR < 0.05),
                "FDR < 0.05")) +
  ylab("-log10P") +
  theme_classic()

#######
sum(df_to_plot$REF_N/df_to_plot$DP < 0.5) #29685
sum(df_to_plot$REF_N/df_to_plot$DP > 0.5) #26342

sum((df_to_plot$REF_N/df_to_plot$DP < 0.5) & (df_to_plot$FDR < 0.05)) #29685
sum((df_to_plot$REF_N/df_to_plot$DP > 0.5) & (df_to_plot$FDR < 0.05)) #29685


for (i in 1:9) {
  df_to_plot <- df_results_type_list[[i]]
  cat(paste0('\ni =  ', i, "\n",
             'REF/DP < 0.5 = ',
             sum(df_to_plot$REF_N/df_to_plot$DP < 0.5), "\n",
             "REF/DP > 0.5 = ",
             sum(df_to_plot$REF_N/df_to_plot$DP > 0.5), "\n",
             "(REF/DP < 0.5) & FDR < 0.05 = ",
             sum((df_to_plot$REF_N/df_to_plot$DP < 0.5) & (df_to_plot$FDR < 0.05)), "\n",
             "(REF/DP > 0.5) & FDR < 0.05 = ",
             sum((df_to_plot$REF_N/df_to_plot$DP > 0.5) & (df_to_plot$FDR < 0.05)), "\n",
             "rs2027349: REF_C =",
             df_to_plot[df_to_plot$ID %in% "rs2027349", ]$REF_C,
             ", DP = ",
             df_to_plot[df_to_plot$ID %in% "rs2027349", ]$DP, 
             ",\n Pval =",
             format(signif(df_to_plot[df_to_plot$ID %in% "rs2027349", ]$pVal, digits = 2),
                    nsmall = 2), 
             ", FDR =",
             format(signif(df_to_plot[df_to_plot$ID %in% "rs2027349", ]$FDR, digits = 2),
                    nsmall = 2),))
}

for (i in 1:9) {
  df_to_plot <- df_results_type_list[[i]]
  print(df_to_plot[df_to_plot$rsID %in% "rs2027349", ])
}




df_assembled_list[["REF_N"]] <-
  str_split(as.character(df_assembled_list[["AD"]]),
            pattern = ",",
            simplify = T)[, 1]

# filter out all homozygous and DP <= 1
sum((df_assembled_list[["DP"]] < 2) & 
      ((df_assembled_list[["GT"]] == "1/1") | (df_assembled_list[["GT"]] == "0/0")), 
    na.rm = T)

length((df_assembled_list[["DP"]] < 2) & 
         ((df_assembled_list[["GT"]] == "1/1") | (df_assembled_list[["GT"]] == "0/0")))
