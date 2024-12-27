# Siwei 05 Aug 2022
# calculate ASoC ratio from the aggregation

# init
library(ggplot2)
library(readr)

library(parallel)
library(stringr)

library(MASS)
library(fitdistrplus)

library(RColorBrewer)
library(grDevices)


##########
input_file_list <-
  list.files(path = "SNP_18_lines_bwa_11Aug2022", 
             pattern = "_4_R\\.txt$",
             full.names = T)
output_path <- "SNP_18_lines_2_annot/"


####### Calc normalised REF/ALT functions, use round()
REF_C_normalised <- function(x) {
  
  ## assemble raw read matrix
  raw_read_matrix <-
    data.frame(REF = as.numeric(unlist(x[seq(from = 7,
                                             to = length(x) - 1,
                                             by = 2)])),
               ALT = as.numeric(unlist(x[seq(from = 8,
                                             to = length(x),
                                             by = 2)])))
  raw_read_matrix$DP <-
    rowSums(raw_read_matrix)
  # print(raw_read_matrix)
  
  ## find the lowest DP count > 1 !!!
  # minDP <- raw_read_matrix$DP - 1
  minDP <- min(raw_read_matrix$DP[raw_read_matrix$DP > 1],
               na.rm = T)
  
  # normalise at individual level, round
  output_matrix <-
    data.frame(REF_N = apply(X = raw_read_matrix,
                             MARGIN = 1,
                             FUN = function(x)(round(x[1] * minDP / x[3]))),
               ALT_N = apply(X = raw_read_matrix,
                             MARGIN = 1,
                             FUN = function(x)(round(x[2] * minDP / x[3]))))
  output_matrix$DP_N <-
    rowSums(output_matrix, 
            na.rm = T)
  
  # take the sum of all values and return
  normalised_count <-
    colSums(output_matrix, 
            na.rm = T)
  
  # return(output_matrix)
  # return(raw_read_matrix)
  return(c(normalised_count))
  
}

######  parApply test ####
i <- 4L

df_raw <- read_delim(input_file_list[i],
                     delim = "\t", escape_double = FALSE,
                     trim_ws = TRUE)

output_file_full_name <-
  str_split(string = input_file_list[i],
            pattern = "\\/",
            simplify = T)[2]
output_file_full_name <-
  str_remove_all(string = output_file_full_name,
                 pattern = '_scATAC_no_VQSR_MAPQ_30_10Aug2022_4_R')

print(output_file_full_name)
print(str_replace_all(string = output_file_full_name,
                      pattern = "\\.txt",
                      replacement = "\\.png"))

df_to_plot <- df_raw
df_to_plot <- df_to_plot[df_to_plot$ID != ".", ]

df_to_plot <- df_to_plot[1:10, ]

# make a 16 core cluster
clust_df_plot_ref <- makeCluster(type = "FORK", 8)
clusterExport(clust_df_plot_ref, "df_to_plot")

# c(df_to_plot$REF_N, df_to_plot$ALT_N, df_to_plot$DP_N) = 
#   parApply(cl = clust_df_plot_ref, X = df_to_plot, MARGIN = 1,
#            FUN = REF_C_normalised)
# df_to_plot[, 43:45] <-
  parApply(cl = clust_df_plot_ref, X = df_to_plot, MARGIN = 1,
           FUN = REF_C_normalised)

stopCluster(clust_df_plot_ref)
rm(clust_df_plot_ref)

#######
i <- 4L
for (i in 1:length(input_file_list)) {
  
  df_raw <- read_delim(input_file_list[i],
                       delim = "\t", escape_double = FALSE,
                       trim_ws = TRUE)
  
  output_file_full_name <-
    str_split(string = input_file_list[i],
              pattern = "\\/",
              simplify = T)[2]
  output_file_full_name <-
    str_remove_all(string = output_file_full_name,
                   pattern = '_scATAC_no_VQSR_MAPQ_30_10Aug2022_4_R')q

  print(output_file_full_name)
  print(str_replace_all(string = output_file_full_name,
                        pattern = "\\.txt",
                        replacement = "\\.png"))
  
  df_to_plot <- df_raw
  df_to_plot <- df_to_plot[df_to_plot$ID != ".", ]
  
  # # make a 16 core cluster #####
  # clust_df_plot_ref <- makeCluster(type = "FORK", 8)
  # clusterExport(clust_df_plot_ref, "df_to_plot")
  # 
  # df_to_plot$REF_C <-
  #   parApply(cl = clust_df_plot_ref, X = df_to_plot, MARGIN = 1,
  #            FUN = function(x)(sum(as.numeric(x[seq(from = 7,
  #                                                   to = ncol(df_to_plot) - 1,
  #                                                   by = 2)]))))
  # 
  # clust_df_plot_alt <- makeCluster(type = "FORK", 8)
  # clusterExport(clust_df_plot_alt, "df_to_plot")
  # 
  # df_to_plot$ALT_C <- 
  #   parApply(cl = clust_df_plot_alt, X = df_to_plot, MARGIN = 1,
  #            FUN = function(x)(sum(as.numeric(x[seq(from = 8,
  #                                                   to = ncol(df_to_plot),
  #                                                   by = 2)]))))
  # stopCluster(clust_df_plot_alt)
  # rm(clust_df_plot_alt)
  # 
  # df_to_plot <- df_to_plot[df_to_plot$REF_C > 1, ]
  # df_to_plot <- df_to_plot[df_to_plot$ALT_C > 1, ]
  # df_to_plot$DP <- df_to_plot$REF_C + df_to_plot$ALT_C
  # 
  # # sum(df_to_plot$DP > 19) 
  # 
  # df_to_plot <- df_to_plot[df_to_plot$DP > 6, ]
  ##### 
  ## !!!!! insert down-sampled, normalised counts here !!!
  ## !!!!! remember the output df is transposed from the original one !!!
  ## !!!!! hence the order of ROWs are REF_N, ALT_N, DP_N
  ## !!!!! attach the transposed ROWs to df_to_plot
  ## !!!!! pVal = binom.test(x=REF_N, n=sum(REF_N, ALT_N)) !!!!!
  
  clust_df_plot_ref <- makeCluster(type = "FORK", 8)
  clusterExport(clust_df_plot_ref, "df_to_plot")
  
  # c(df_to_plot$REF_N, df_to_plot$ALT_N, df_to_plot$DP_N) = 
  #   parApply(cl = clust_df_plot_ref, X = df_to_plot, MARGIN = 1,
  #            FUN = REF_C_normalised)
  # df_to_plot[, 43:45] <-
  returned_df <-
    parApply(cl = clust_df_plot_ref, X = df_to_plot, MARGIN = 1,
             FUN = REF_C_normalised)
  
  stopCluster(clust_df_plot_ref)
  rm(clust_df_plot_ref)
  
  df_to_plot$REF_N <- returned_df[1, ]
  df_to_plot$ALT_N <- returned_df[2, ]
  df_to_plot$DP_N <- returned_df[3, ]
  # df_to_plot$DP_N <- 
  #   df_to_plot$REF_N +
  #   df_to_plot$ALT_N
  
  # df_to_plot_bkup <- df_to_plot
  # df_to_plot <- df_to_plot_bkup
  df_to_plot <-
    df_to_plot[df_to_plot$DP_N > 9, ]
  
  #### calc p and FDR
  
  clust_df_plot_use <- makeCluster(type = "FORK", 8)
  clusterExport(clust_df_plot_use, "df_to_plot")
  
  df_to_plot$pVal <- 
    # parApply(cl = clust_df_plot_use, X = df_to_plot, MARGIN = 1,
    #          FUN = function(x)(binom.test(x = as.numeric(x[ncol(df_to_plot) - 1]),
    #                                       n = sum(as.numeric(x[ncol(df_to_plot) - 1]), 
    #                                               as.numeric(x[ncol(df_to_plot) - 2])),
    #                                       p = 0.5,
    #                                       alternative = "t")$p.value))
  parApply(cl = clust_df_plot_use, X = df_to_plot, MARGIN = 1,
           FUN = function(x)(binom.test(x = as.numeric(x[ncol(df_to_plot) - 2]),
                                        n = as.numeric(x[ncol(df_to_plot)]),
                                        p = 0.5,
                                        alternative = "t")$p.value))
  stopCluster(clust_df_plot_use)
  rm(clust_df_plot_use)
  
  
  df_to_plot$FDR <- p.adjust(p = df_to_plot$pVal,
                             method = "fdr")
  df_to_plot_output <- df_to_plot[, c(1:6, (ncol(df_to_plot) - 4):ncol(df_to_plot))]
  df_to_plot_output <- df_to_plot_output[order(df_to_plot_output$pVal), ]
  ### cut only FDR < 0.05 for annotation
  # df_to_plot_output <- df_to_plot_output[df_to_plot_output$FDR < 0.05, ]
  
  # sum(df_to_plot$FDR < 0.05)
  
  ### make PNG plot and send out
  ## !! assign time and cell type !!
  # png(filename = paste(output_path,
  #                      str_replace_all(string = output_file_full_name,
  #                                      pattern = "\\.txt",
  #                                      replacement = "\\.png")),
  #     width = 640, height = 480, 
  #     bg = "white")
  
  # fig_to_plot <- 
    ggplot(df_to_plot,
         aes(x = (REF_N/DP_N),
             y = (0 - log10(pVal)))) +
    geom_point(size = 0.5,
               aes(colour = factor(ifelse(FDR < 0.05, 
                                          "FDR < 0.05", 
                                          "FDR > 0.05")))) +
    scale_color_manual(values = c("red", "black")) +
    labs(colour = "FDR value") +
    ggtitle(paste(output_file_full_name, ", DP >= 10, 5 lines\n",
                  # "Duplicates from picard MarkDuplicates\n",
                  # "Duplicates marked by 10x Cell Ranger\n",
                  nrow(df_to_plot),
                  "SNPs, of",
                  sum(df_to_plot$FDR < 0.05),
                  "FDR < 0.05")) +
    ylab("-log10P") +
    theme_classic()
  
  # print(fig_to_plot)

  # dev.off()
  ###
    
  
  print(df_to_plot_output[df_to_plot_output$ID %in% "rs2027349", ])
  
  sum((df_to_plot$REF_N / df_to_plot$DP_N > 0.5) & df_to_plot$FDR < 0.05)
  sum((df_to_plot$REF_N / df_to_plot$DP_N < 0.5) & df_to_plot$FDR < 0.05)
  
  # use either REF/ALT, 325:59
  
  
  # output_file_full_name <-
  #   paste(output_path,
  #         output_file_full_name,
  #         sep = "")
  # 
  # df_to_annot <- df_to_plot_output[, c(1, 2, 2, 4, 5, 3, 6:10)]
  # 
  # write.table(df_to_annot,
  #             file = str_replace_all(string = output_file_full_name,
  #                                    pattern = "\\.txt",
  #                                    replacement = "\\.avinput"),
  #             quote = F, sep = "\t",
  #             row.names = F, col.names = F)
  
}





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
  ggtitle(paste("NEFM pos glut, 0hr, DP >= 20, 18 lines\n",
                # "Duplicates from picard MarkDuplicates\n",
                "Duplicates marked by 10x Cell Ranger\n",
                nrow(df_to_plot),
                "SNPs, of",
                sum(df_to_plot$FDR < 0.05),
                "FDR < 0.05")) +
  ylab("-log10P") +
  theme_classic()



sum(df_to_plot$REF_C/df_to_plot$DP < 0.5) #29685
sum(df_to_plot$REF_C/df_to_plot$DP > 0.5) #26342

df_test <- df_to_plot[, c(1,2,3)]
df_test <- df_test[, c(1,2,3,3)]
