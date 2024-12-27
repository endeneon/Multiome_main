# Chuxuan Li 04/28/2022
# Using the output of 26Apr2022_call_DA_peaks_on_ATAC_new_obj.R, 1) plot volcano
#plots for all peaks, 2) plot dotplots showing peaks near response genes 
# 08/01/2022 repurposed for processing da peaks minpct = 0.01 results

# init ####
library(readr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(graphics)
library(stringr)

set.seed(2022)

names_1v0 <- c("GABA_1v0hr", "nmglut_1v0hr", "npglut_1v0hr")
names_6v0 <- c("GABA_6v0hr", "nmglut_6v0hr", "npglut_6v0hr")
colors <- c("steelblue3", "grey", "indianred3")

labels <- c("chr11-27770213-27771025", "chr14-75236946-75238215",
            "chr11-66416626-66417402", "chr5-138459042-138459607")
names(labels) <- c("BDNF", "FOS", "NPAS4", "EGR1")

# da_peaks_by_time_list_1v0 <- list(GABA_1v0, NEFM_neg_glut_1v0, NEFM_pos_glut_1v0)
# da_peaks_by_time_list_6v0 <- list(GABA_6v0, NEFM_neg_glut_6v0, NEFM_pos_glut_6v0)
da_peaks_by_time_list_1v0 <- vector("list", length(fit_all))
da_peaks_by_time_list_6v0 <- vector("list", length(fit_all))
for (i in 1:length(fit_all)){
  da_peaks_by_time_list_1v0[[i]] <- topTable(fit_all[[i]], coef = "onevszero", p.value = Inf, 
                                             sort.by = "P", number = Inf)
  da_peaks_by_time_list_6v0[[i]] <- topTable(fit_all[[i]], coef = "sixvszero", p.value = Inf, 
                                             sort.by = "P", number = Inf)
}

# volcano plots ####
for (i in 1:3){
  name1 <- names_1v0[i]
  print(name1)
  da_peaks_by_time_list_1v0[[i]]$peak <- rownames(da_peaks_by_time_list_1v0[[i]])
  da_peaks_by_time_list_1v0[[i]]$genes <- NULL
  df <- da_peaks_by_time_list_1v0[[i]][, c("peak", "AveExpr", "logFC", "P.Value", "adj.P.Val")]
  df$significance <- "nonsig"
  df$significance[df$adj.P.Val < 0.05 & df$logFC > 0] <- "pos"
  df$significance[df$adj.P.Val < 0.05 & df$logFC < 0] <- "neg"
  df$neg_log_p_val <- (0 - log10(df$P.Value))
  #df$labelling <- ""
  # for (j in 1:length(labels)){
  #   df$labelling[df$peak %in% labels[j]] <- names(labels[j])
  # }
  # print(unique(df$labelling))
  
  # plot
#  file_name <- paste0("./DApeaks_plots/volcano_plots/", name1, ".pdf")
  file_name <- paste0("./pseudobulk_da_volcano_plots_updated/", name1, ".pdf")
  pdf(file = file_name)
  p <- ggplot(data = df,
              aes(x = logFC, 
                  y = neg_log_p_val, 
                  #label = labelling,
                  color = significance)) + 
    geom_point() +
    scale_color_manual(values = colors) +
    theme_minimal() +
    ggtitle(str_replace_all(name1, "_", " ")) +
    # geom_text_repel(box.padding = unit(0.05, 'lines'),
    #                 min.segment.length = 0,
    #                 force = 0,
    #                 force_pull = 10,
    #                 max.overlaps = 100000) +
    xlim(-5, 5) 
  print(p)
  dev.off()
  
  name6 <- names_6v0[i]
  print(name6)
  da_peaks_by_time_list_6v0[[i]]$peak <- rownames(da_peaks_by_time_list_6v0[[i]])
  da_peaks_by_time_list_6v0[[i]]$genes <- NULL
  df <- da_peaks_by_time_list_6v0[[i]][, c("peak", "AveExpr", "logFC", "P.Value", "adj.P.Val")]
  df$significance <- "nonsig"
  df$significance[df$adj.P.Val < 0.05 & df$logFC > 0] <- "pos"
  df$significance[df$adj.P.Val < 0.05 & df$logFC < 0] <- "neg"
  df$neg_log_p_val <- (0 - log10(df$P.Value))
  #df$labelling <- ""
  # for (j in 1:length(labels)){
  #   df$labelling[df$peak %in% labels[j]] <- names(labels[j])
  # }
  #print(unique(df$labelling))
  
  #  file_name <- paste0("./DApeaks_plots/volcano_plots/", name1, ".pdf")
  file_name <- paste0("./pseudobulk_da_volcano_plots_updated/", name6, ".pdf")
  pdf(file = file_name)
  p <- ggplot(data = df,
              aes(x = logFC, 
                  y = neg_log_p_val, 
                  #label = labelling,
                  color = significance)) + 
    geom_point() +
    scale_color_manual(values = colors) +
    theme_minimal() +
    ggtitle(str_replace_all(name6, "_", " ")) +
    # geom_text_repel(box.padding = unit(0.05, 'lines'),
    #                 min.segment.length = 0,
    #                 force = 0,
    #                 force_pull = 10,
    #                 max.overlaps = 100000) +
    xlim(-5, 5) 
  print(p)
  dev.off()
}
 # summary table ####
dapeak_counts <- matrix(nrow = 6, ncol = 4, 
                        dimnames = list(c("GABA 1v0", "GABA 6v0", 
                                          "nmglut 1v0", "nmglut 6v0",
                                          "npglut 1v0", "npglut 6v0"),
                                        c("total", "significant", "up", "down")))
for (i in 1:length(da_peaks_by_time_list_6v0)){
  print(i)
  dapeak_counts[2*i-1, 1] <- nrow(da_peaks_by_time_list_1v0[[i]])
  dapeak_counts[2*i, 1] <- nrow(da_peaks_by_time_list_6v0[[i]])
  dapeak_counts[2*i-1, 2] <- sum(da_peaks_by_time_list_1v0[[i]]$adj.P.Val < 0.05)
  dapeak_counts[2*i, 2] <- sum(da_peaks_by_time_list_6v0[[i]]$adj.P.Val < 0.05)
  dapeak_counts[2*i-1, 3] <- sum(da_peaks_by_time_list_1v0[[i]]$adj.P.Val < 0.05 & 
                                   da_peaks_by_time_list_1v0[[i]]$logFC > 0)
  dapeak_counts[2*i, 3] <- sum(da_peaks_by_time_list_6v0[[i]]$adj.P.Val < 0.05 & 
                                 da_peaks_by_time_list_6v0[[i]]$logFC > 0)
  dapeak_counts[2*i-1, 4] <- sum(da_peaks_by_time_list_1v0[[i]]$adj.P.Val < 0.05 & 
                                   da_peaks_by_time_list_1v0[[i]]$logFC < 0)
  dapeak_counts[2*i, 4] <- sum(da_peaks_by_time_list_6v0[[i]]$adj.P.Val < 0.05 & 
                                 da_peaks_by_time_list_6v0[[i]]$logFC < 0)
}
write.table(dapeak_counts, file = "dapeak_counts_summary_minpct0.01.csv", 
            quote = F, sep = ",", row.names = T, col.names = T)

# dotplot dfs ####
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NEFM- glut", "NEFM+ glut"),
                             avg_log2FC = rep_len(NA, 3),
                             pct.1 = rep_len(NA, 3),
                             pct.2 = rep_len(NA, 3))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NEFM- glut", "NEFM+ glut"),
                             avg_log2FC = rep_len(NA, 3),
                             pct.1 = rep_len(NA, 3),
                             pct.2 = rep_len(NA, 3))

#ATP2A2 ####
#chr12   110269448       110270452
for (i in 1:3){
  print(i)
  one <- da_peaks_by_time_list_1v0[[i]][str_detect(da_peaks_by_time_list_1v0[[i]]$peak, "chr12-11026944"), ]
  six <- da_peaks_by_time_list_6v0[[i]][str_detect(da_peaks_by_time_list_6v0[[i]]$peak, "chr12-11026944"), ]
  if (length(one$p_val) != 0){
    df_to_plot_1v0$avg_log2FC[i] <- one$avg_log2FC
    df_to_plot_1v0$pct.1[i] <- one$pct.1
    df_to_plot_1v0$pct.2[i] <- one$pct.2
  }
  if (length(six$p_val) != 0){
    df_to_plot_6v0$avg_log2FC[i] <- six$avg_log2FC
    df_to_plot_6v0$pct.1[i] <- six$pct.1
    df_to_plot_6v0$pct.2[i] <- six$pct.2
  }
}

#BDNF ####
#chr11   27770214        27771025
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NEFM- glut", "NEFM+ glut"),
                             avg_log2FC = c(NA, NA, 0.0368462596774217),
                             pct.1 = c(NA, NA, 0.035),
                             pct.2 = c(NA, NA, 0.011))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NEFM- glut", "NEFM+ glut"),
                             avg_log2FC = c(NA, NA, 0.0503842655550294),
                             pct.1 = c(NA, NA, 0.044),
                             pct.2 = c(NA, NA, 0.011))

#EGR1 ####
#chr5    138459043       138459607
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NEFM- glut", "NEFM+ glut"),
                             avg_log2FC = c(NA, NA, 0.01327452),
                             pct.1 = c(NA, NA, 0.034),
                             pct.2 = c(NA, NA, 0.025))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NEFM- glut", "NEFM+ glut"),
                             avg_log2FC = c(NA, NA, 0.004688141),
                             pct.1 = c(NA, NA, 0.029),
                             pct.2 = c(NA, NA, 0.025))

#FOS ####
#chr14   75236947        75238215
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NEFM- glut", "NEFM+ glut"),
                             avg_log2FC = c(0.05834332, 0.1665399, 0.1444198),
                             pct.1 = c(0.170, 0.219, 0.2),
                             pct.2 = c(0.134, 0.116, 0.113))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NEFM- glut", "NEFM+ glut"),
                             avg_log2FC = c(0.05812602, 0.09849978, 0.1024791),
                             pct.1 = c(0.171),
                             pct.2 = c(0.134, 0.116, 0.113))

#ETF1 ####
#chr5    138547883       138548351
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NEFM- glut", "NEFM+ glut"),
                             avg_log2FC = c(NA, 0.006100954, 0.008033972),
                             pct.1 = c(NA, 0.023, 0.031),
                             pct.2 = c(NA, 0.019, 0.025))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NEFM- glut", "NEFM+ glut"),
                             avg_log2FC = c(NA, 0.002119181, 0.003675232),
                             pct.1 = c(NA, 0.021, 0.022),
                             pct.2 = c(NA, 0.019, 0.025))

#MMP1  ####
#chr11   102797848       102798589
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NEFM- glut", "NEFM+ glut"),
                             avg_log2FC = c(NA, NA, 0.01898081),
                             pct.1 = c(NA, NA, 0.032),
                             pct.2 = c(NA, NA, 0.019))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NEFM- glut", "NEFM+ glut"),
                             avg_log2FC = c(NA, NA, 0.0410645),
                             pct.1 = c(NA, NA, 0.047),
                             pct.2 = c(NA, NA, 0.019))
for (i in 1:3){
  print(i)
  print(da_peaks_by_time_list_1v0[[i]][str_detect(da_peaks_by_time_list_1v0[[i]]$peak, "10279784"), ])
  print(da_peaks_by_time_list_6v0[[i]][str_detect(da_peaks_by_time_list_6v0[[i]]$peak, "10279784"), ])
}

#NPAS4 ####
#chr11   66416627        66417402 
for (i in 1:3){
  print(i)
  one <- da_peaks_by_time_list_1v0[[i]][str_detect(da_peaks_by_time_list_1v0[[i]]$peak, "chr11-66416626"), ]
  six <- da_peaks_by_time_list_6v0[[i]][str_detect(da_peaks_by_time_list_6v0[[i]]$peak, "chr11-66416626"), ]
  if (length(one$p_val) != 0){
    df_to_plot_1v0$avg_log2FC[i] <- one$avg_log2FC
    df_to_plot_1v0$pct.1[i] <- one$pct.1
    df_to_plot_1v0$pct.2[i] <- one$pct.2
  }
  if (length(six$p_val) != 0){
    df_to_plot_6v0$avg_log2FC[i] <- six$avg_log2FC
    df_to_plot_6v0$pct.1[i] <- six$pct.1
    df_to_plot_6v0$pct.2[i] <- six$pct.2
  }
}

#NR4A1 ####
#chr12   52036490        52037843
for (i in 1:3){
  print(i)
  one <- da_peaks_by_time_list_1v0[[i]][str_detect(da_peaks_by_time_list_1v0[[i]]$peak, "chr12-52036489"), ]
  six <- da_peaks_by_time_list_6v0[[i]][str_detect(da_peaks_by_time_list_6v0[[i]]$peak, "chr12-52036489"), ]
  if (length(one$p_val) != 0){
    df_to_plot_1v0$avg_log2FC[i] <- one$avg_log2FC
    df_to_plot_1v0$pct.1[i] <- one$pct.1
    df_to_plot_1v0$pct.2[i] <- one$pct.2
  }
  if (length(six$p_val) != 0){
    df_to_plot_6v0$avg_log2FC[i] <- six$avg_log2FC
    df_to_plot_6v0$pct.1[i] <- six$pct.1
    df_to_plot_6v0$pct.2[i] <- six$pct.2
  }
}

#RGS2 ####
#chr1    192788011       192788589 
for (i in 1:3){
  print(i)
  one <- da_peaks_by_time_list_1v0[[i]][str_detect(da_peaks_by_time_list_1v0[[i]]$peak, "chr1-192788010"), ]
  six <- da_peaks_by_time_list_6v0[[i]][str_detect(da_peaks_by_time_list_6v0[[i]]$peak, "chr1-192788010"), ]
  if (length(one$p_val) != 0){
    df_to_plot_1v0$avg_log2FC[i] <- one$avg_log2FC
    df_to_plot_1v0$pct.1[i] <- one$pct.1
    df_to_plot_1v0$pct.2[i] <- one$pct.2
  }
  if (length(six$p_val) != 0){
    df_to_plot_6v0$avg_log2FC[i] <- six$avg_log2FC
    df_to_plot_6v0$pct.1[i] <- six$pct.1
    df_to_plot_6v0$pct.2[i] <- six$pct.2
  }
}

#ZSWIM6 ####
#chr5    61394120        61394736
for (i in 1:3){
  print(i)
  one <- da_peaks_by_time_list_1v0[[i]][str_detect(da_peaks_by_time_list_1v0[[i]]$peak, "chr5-61394119"), ]
  six <- da_peaks_by_time_list_6v0[[i]][str_detect(da_peaks_by_time_list_6v0[[i]]$peak, "chr5-61394119"), ]
  if (length(one$p_val) != 0){
    df_to_plot_1v0$avg_log2FC[i] <- one$avg_log2FC
    df_to_plot_1v0$pct.1[i] <- one$pct.1
    df_to_plot_1v0$pct.2[i] <- one$pct.2
  }
  if (length(six$p_val) != 0){
    df_to_plot_6v0$avg_log2FC[i] <- six$avg_log2FC
    df_to_plot_6v0$pct.1[i] <- six$pct.1
    df_to_plot_6v0$pct.2[i] <- six$pct.2
  }
}

#PPP1R13B ####
#chr14   103803430       103804242
for (i in 1:3){
  print(i)
  one <- da_peaks_by_time_list_1v0[[i]][str_detect(da_peaks_by_time_list_1v0[[i]]$peak, "chr14-1038034"), ]
  six <- da_peaks_by_time_list_6v0[[i]][str_detect(da_peaks_by_time_list_6v0[[i]]$peak, "chr14-1038034"), ]
  if (length(one$p_val) != 0){
    df_to_plot_1v0$avg_log2FC[i] <- one$avg_log2FC
    df_to_plot_1v0$pct.1[i] <- one$pct.1
    df_to_plot_1v0$pct.2[i] <- one$pct.2
  }
  if (length(six$p_val) != 0){
    df_to_plot_6v0$avg_log2FC[i] <- six$avg_log2FC
    df_to_plot_6v0$pct.1[i] <- six$pct.1
    df_to_plot_6v0$pct.2[i] <- six$pct.2
  }
}

# make dotplot ####
df_to_plot <- rbind(df_to_plot_1v0,
                    df_to_plot_6v0)
df_to_plot$source <- c(rep_len("1v0hr", length.out = length(df_to_plot_1v0$cell_type)),
                       rep_len("6v0hr", length.out = length(df_to_plot_6v0$cell_type)))

df_to_plot$cell_type <- factor(df_to_plot$cell_type, 
                               levels = rev(c("GABA", "NEFM- glut", "NEFM+ glut")))


ggplot(df_to_plot,
       aes(x = as.factor(source),
           y = as.factor(cell_type),
           size = pct.1 * 100,
           fill = 2 ^ avg_log2FC * pct.1 / pct.2)) +
  geom_point(shape = 21) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("#ffdbd9", "red", "darkred")) +
  scale_size(range = c(5, 7)) +
  theme_grey(base_size = 10) +
  ggtitle("ATP2A2")


