# Chuxuan Li 12/12/2022
# Plot volcano plots and dotplots for DESeq2 results with 025 data

# init ####
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(stringr)
library(readr)

# load data ####
pathlist <- list.files("./Analysis_part2_GRCh38/DESeq2_DE_results/full_no_filtering", full.names = T,
                       pattern = "1v0|6v0")
pathlist
reslist <- vector('list', length(pathlist))
for (i in 1:length(reslist)){
  reslist[[i]] <- as.data.frame(read_csv(pathlist[i]))
}

# volcano plots ####
namelist <- str_extract(pathlist, "[A-Za-z]+_[1|6]v0")
namelist
setwd("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/DESeq2_DE_results/volcano_plots")
for (i in 1:length(reslist)){
  print(namelist[i])
  res <- reslist[[i]]
  res$significance <- "nonsignificant"
  res$significance[res$padj < 0.05 & res$log2FoldChange > 0] <- "up"
  res$significance[res$padj < 0.05 & res$log2FoldChange < 0] <- "down"
  unique(res$significance)
  res$significance <- factor(res$significance, levels = c("up", "nonsignificant", "down"))
  res$neg_log_pval <- (0 - log2(res$pvalue))
  res$labelling <- ""
  for (j in c("FOS", "NPAS4", "IGF1", "NR4A1","FOSB", "EGR1", "BDNF")){
    res$labelling[res$gene %in% j] <- j
  }
  png(paste0(namelist[i], "_volcano_plot.png"))
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

# bar graph ####
early_list <- c("FOS", "FOSB", "NPAS4", "NR4A1")
late_list <- c("BDNF", "VGF")

df_to_plot <- rbind(res_GABA_1v0[late_list, ],
                    res_GABA_6v0[late_list, ],
                    res_npglut_1v0[late_list, ],
                    res_npglut_6v0[late_list, ],
                    res_nmglut_1v0[late_list, ],
                    res_nmglut_6v0[late_list, ])
df_to_plot <- as.data.frame(df_to_plot)

df_to_plot$time <- c(rep_len("1v0", length(late_list)), rep_len("6v0", length(late_list)), 
                     rep_len("1v0", length(late_list)), rep_len("6v0", length(late_list)),
                     rep_len("1v0", length(late_list)), rep_len("6v0", length(late_list)))                    
df_to_plot$cell.type <- c(rep_len("GABA", length(late_list) * 2), 
                          rep_len("NEFM- glut", length(late_list) * 2),
                          rep_len("NEFM+ glut", length(late_list) * 2))
df_to_plot$gene.name <- c(late_list, late_list, 
                          late_list, late_list, 
                          late_list, late_list)

# for late, mark if they are excitatory or inhibitory
gene_labeller <- as_labeller(c('BDNF' = "BDNF (excitatory)", 
                               'IGF1' = "IGF1 (inhibitory)", 
                               'VGF' = "VGF (shared)"))
ggplot(df_to_plot, aes(x = cell.type, 
                       y = log2FoldChange,
                       color = cell.type,
                       fill = time,
                       group = time,
                       ymax = log2FoldChange - 1/2*lfcSE, 
                       ymin = log2FoldChange + 1/2*lfcSE
)) + 
  xlab("") +
  ylab("log2(FC)") +
  geom_col(position = position_dodge(0.6),
           width = 0.5,
           color = "black") +
  geom_errorbar(color = "black",
                position = position_dodge(0.6),
                width = 0.5) +
  facet_grid(labeller = gene_labeller,
    cols = vars(gene.name)) +
  scale_fill_manual(values = c("indianred", "steelblue")) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("differential expression of late response genes - DESeq2 results")
