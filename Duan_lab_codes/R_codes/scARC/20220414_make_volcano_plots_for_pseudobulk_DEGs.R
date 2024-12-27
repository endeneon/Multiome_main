# Chuxuan Li 04/14/2022
# Use 18-line only pseudobulk DEG lists to plot volcano plots

# init ####
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(stringr)
library(readr)

# load data ####
pathlist <- list.files("./pseudobulk_DE/filtered_by_basemean_only//", full.names = T,
                       pattern = "1v0|6v0")
pathlist
reslist <- vector('list', length(pathlist))
for (i in 1:length(reslist)){
  reslist[[i]] <- read_csv(pathlist[i])
}
# make plots ####
namelist <- str_extract(pathlist, "[A-Za-z]+_[1|6]v0")
namelist
setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/pseudobulk_DE/volcano_plots")

for (i in 1:length(reslist)){
  print(namelist[i])
  res <- reslist[[i]]
  #res$gene.symbol <- rownames(res)
  res$significance <- "nonsignificant"
  res$significance[res$q_value < 0.05 & res$log2FoldChange > 0] <- "up"
  res$significance[res$q_value < 0.05 & res$log2FoldChange < 0] <- "down"
  unique(res$significance)
  res$significance <- factor(res$significance, levels = c("up", "nonsignificant", "down"))
  res$neg_log_pval <- (0 - log2(res$pvalue))
  res$labelling <- ""
  for (j in c("FOS", "NPAS4", "IGF1", "NR4A1","FOSB", "EGR1", "BDNF")){
    res$labelling[res$gene %in% j] <- j
  }
  jpeg(paste0(namelist[i], "_volcano_plot.jpeg"))
  p <- ggplot(data = as.data.frame(res),
              aes(x = log2FoldChange, 
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
    xlim(c(-10, 10)) +
    ggtitle(str_replace_all(namelist[i], "_", " "))
  print(p)
  dev.off()
}
