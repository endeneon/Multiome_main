# Siwei 13 Sept 2022
# calculate ASoC ratio from the aggregation
# of up/down sampled results

# init
library(ggplot2)
library(readr)

library(parallel)
library(stringr)

library(MASS)

library(RColorBrewer)
library(grDevices)

input_file_list <-
  list.files(path = "downsize_no_upsampling", 
             pattern = "14Sept2022\\.txt$",
             full.names = T)
output_path <- "downsize_no_upsampling/output/"

input_file_list <-
  list.files(path = "same_size_w_upsampling", 
             pattern = "12Sept2022\\.txt$",
             full.names = T)
output_path <- "same_size_w_upsampling/output"




i <- 1L
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
                   pattern = '_downsampled_median_downsampled_no_upsampling_14Sept2022')
  
  print(output_file_full_name)
  print(str_replace_all(string = output_file_full_name,
                        pattern = "\\.txt",
                        replacement = "\\.png"))
  
  df_to_plot <- df_raw
  df_to_plot <- df_to_plot[df_to_plot$ID != ".", ]
  
  # make a 16 core cluster
  clust_df_plot_ref <- makeCluster(type = "FORK", 8)
  clusterExport(clust_df_plot_ref, "df_to_plot")
  
  df_to_plot$REF_C <- 
    parApply(cl = clust_df_plot_ref, X = df_to_plot, MARGIN = 1, 
             FUN = function(x)(sum(as.numeric(x[seq(from = 7,
                                                    to = ncol(df_to_plot) - 1,
                                                    by = 2)]))))
  stopCluster(clust_df_plot_ref)
  rm(clust_df_plot_ref)
  
  clust_df_plot_alt <- makeCluster(type = "FORK", 8)
  clusterExport(clust_df_plot_alt, "df_to_plot")
  
  df_to_plot$ALT_C <- 
    parApply(cl = clust_df_plot_alt, X = df_to_plot, MARGIN = 1,
             FUN = function(x)(sum(as.numeric(x[seq(from = 8,
                                                    to = ncol(df_to_plot),
                                                    by = 2)]))))
  stopCluster(clust_df_plot_alt)
  rm(clust_df_plot_alt)
  
  df_to_plot <- df_to_plot[df_to_plot$REF_C > df_to_plot$NCALLED * 2, ]
  df_to_plot <- df_to_plot[df_to_plot$ALT_C > df_to_plot$NCALLED * 2, ]
  df_to_plot$DP <- df_to_plot$REF_C + df_to_plot$ALT_C
  
  # sum(df_to_plot$DP > 19) 
  
  df_to_plot <- df_to_plot[df_to_plot$DP > 29, ]
  df_to_plot <- df_to_plot[df_to_plot$DP < 1000, ]
  
  
  clust_df_plot_use <- makeCluster(type = "FORK", 8)
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
                       sep = ""),
      width = 400, height = 400, 
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
    ggtitle(paste(output_file_full_name,"\n", "DP >= 30, 18 lines\n",
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
  ###
  
  print(df_to_plot_output[df_to_plot_output$ID %in% "rs2027349", ])
  
  output_file_full_name <-
    paste(output_path,
          output_file_full_name,
          sep = "")
  
  df_to_annot <- df_to_plot_output[, c(1, 2, 2, 4, 5, 3, 6:11)]
  
  write.table(df_to_annot,
              file = str_replace_all(string = output_file_full_name,
                                     pattern = "\\.txt",
                                     replacement = "\\.avinput"),
              quote = F, sep = "\t",
              row.names = F, col.names = F)
  
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

sum((df_to_plot$REF_C/df_to_plot$DP < 0.5) & (df_to_plot$FDR < 0.05)) #29685

sum((df_to_plot$REF_C/df_to_plot$DP > 0.5) & (df_to_plot$FDR < 0.05)) #29685


sum(df_to_plot$REF_C/df_to_plot$DP > 0.5) #26342


df_test <- df_to_plot[, c(1,2,3)]
df_test <- df_test[, c(1,2,3,3)]


# 22484 of 94960

df_no_upsampling <- df_to_plot_output
