# Chuxuan Li 12/09/2022
# Plot volcano plots and dotplots for limma results with PRS

# init ####
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(stringr)
library(readr)

# load data ####
pathlist <- list.files("./sep_by_celltype_time_DE_results/unfiltered_by_padj", full.names = T,
                       pattern = "all_DEGs")
pathlist
reslist <- vector('list', length(pathlist))
for (i in 1:length(reslist)){
  reslist[[i]] <- as.data.frame(read_csv(pathlist[i]))
}

# make plots ####
namelist <- str_extract(pathlist, "[A-Za-z]+_[0|1|6]hr")
namelist
setwd("/nvmefs/scARC_Duan_018/SZ_PRS_DE_analysis/sep_by_celltype_time_DE_results/volcano_plots/")
for (i in 1:length(reslist)){
  print(namelist[i])
  res <- reslist[[i]]
  res$significance <- "nonsignificant"
  res$significance[res$P.Value < 0.05 & res$logFC > 0] <- "beta > 0 (p-value < 0.05)"
  res$significance[res$P.Value < 0.05 & res$logFC < 0] <- "beta < 0 (p-value < 0.05)"
  unique(res$significance)
  res$significance <- factor(res$significance, levels = c("beta > 0 (p-value < 0.05)", 
                                                          "nonsignificant", 
                                                          "beta < 0 (p-value < 0.05)"))
  res$neg_log_pval <- (0 - log2(res$P.Value))
  res$labelling <- ""
  
  png(paste0(namelist[i], "_volcano_plot.png"))
  p <- ggplot(data = as.data.frame(res),
              aes(x = logFC, 
                  y = neg_log_pval, 
                  color = significance,
                  label = labelling)) + 
    geom_point(size = 0.2) +
    scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
    theme_minimal() +
    geom_text_repel(box.padding = unit(0.05, 'lines'),
                    min.segment.length = 0,
                    force = 0.2,
                    max.overlaps = Inf,
                    force_pull = 0.1,
                    show.legend = F) +
    xlim(c(-2.2, 2.2)) +
    ggtitle(str_replace_all(namelist[i], "_", " "))
  print(p)
  dev.off()
}

