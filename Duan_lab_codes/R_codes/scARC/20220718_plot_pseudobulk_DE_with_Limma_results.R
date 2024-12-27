# Chuxuan Li 07/18/2022
# Plot the results for Limma generated 18-line DEG results on pseudobulk count
#matrices (raw and adjusted by combat)

# init ####
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(colorRamps)
library(readr)
library(stringr)
library(dplyr)

# read data ####
# setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_pseudobulk/20line_limma/")

raw_res <- list.files(path = ".", 
                      pattern = "all_DEGs.csv", full.names = T, recursive = T, include.dirs = F)
# raw_res <- list.files(path = "./pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm_9in18_res/", 
#            pattern = "all_DEGs.csv", full.names = T, recursive = T, include.dirs = F)
# adj_res <- list.files(path = "./pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm_12in18_res/", 
#                       pattern = "all_DEGs.csv", full.names = T, recursive = T, include.dirs = F)
deg_lst_raw <- vector("list", length(raw_res))
# deg_lst_adj <- vector("list", length(adj_res))

for (i in 1:length(deg_lst_raw)) {
  deg_lst_raw[[i]] <- read_csv(raw_res[i], show_col_types = F)
  # deg_lst_adj[[i]] <- read_csv(adj_res[i], show_col_types = F)
  print(str_extract(raw_res[i], "[A-Za-z]+_[A-Za-z]+_[1|6]v[0|1|6]"))
}
length(deg_lst_raw[[6]]$genes)

# volcano plots ####
namelist <- str_extract(raw_res, "[A-Za-z]+_[A-Za-z]+_[1|6]v[0|1|6]")
namelist <- str_remove(namelist, "_res")
namelist

for (i in 1:length(deg_lst_raw)){
  res <- deg_lst_raw[[i]]
  #res$gene.symbol <- rownames(res)
  res$significance <- "nonsignificant"
  res$significance[res$adj.P.Val < 0.05 & res$logFC > 0] <- "up"
  res$significance[res$adj.P.Val < 0.05 & res$logFC < 0] <- "down"
  unique(res$significance)
  res$significance <- factor(res$significance, levels = c("up", "nonsignificant", "down"))
  res$neg_log_pval <- (0 - log2(res$P.Value))
  res$labelling <- ""
  for (j in c("FOS", "NPAS4", "IGF1", "NR4A1","FOSB", "EGR1", "BDNF")){
    res$labelling[res$gene %in% j] <- j
  }
  # pdf(paste0("./pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm_9in18_res/",
  #             namelist[i], "_volcano_plot.pdf"))
  pdf(paste0(namelist[i], "_volcano_plot.pdf"))
  
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
    xlim(c(-10, 10)) +
    ggtitle(str_replace_all(namelist[i], "_", " "))
  print(p)
  dev.off()
}

# bargraph ####
early_list <- c("FOS", "FOSB", "NPAS4", "NR4A1")
late_list <- c("BDNF", "IGF21", "VGF")
jubao_list <- c("BCL11B", "IMMP2L", "SGCD", "SNAP91")
df_to_plot <- rbind(deg_lst_raw[[1]][deg_lst_raw[[1]]$genes %in% jubao_list, ],
                    deg_lst_raw[[2]][deg_lst_raw[[2]]$genes %in% jubao_list, ],
                    deg_lst_raw[[4]][deg_lst_raw[[4]]$genes %in% jubao_list, ],
                    deg_lst_raw[[5]][deg_lst_raw[[5]]$genes %in% jubao_list, ],
                    deg_lst_raw[[7]][deg_lst_raw[[7]]$genes %in% jubao_list, ],
                    deg_lst_raw[[8]][deg_lst_raw[[8]]$genes %in% jubao_list, ])
df_to_plot <- as.data.frame(df_to_plot)
df_to_plot$SE <- df_to_plot$logFC/df_to_plot$t
df_to_plot$time <- rep(rep(c("1v0", "6v0"), each = sum(deg_lst_raw[[1]]$genes %in% jubao_list)), times = 3)                  
df_to_plot$cell.type <- c(rep_len("GABA", sum(deg_lst_raw[[1]]$genes %in% jubao_list) * 2), 
                          rep_len("NEFM- glut", sum(deg_lst_raw[[1]]$genes %in% jubao_list) * 2),
                          rep_len("NEFM+ glut", sum(deg_lst_raw[[1]]$genes %in% jubao_list) * 2))
df_to_plot <- df_to_plot %>%
  group_by(genes) %>%
  arrange(df_to_plot, .by_group = TRUE)

# for jubao, mark if they are excitatory or inhibitory
gene_labeller <- as_labeller(c('BDNF' = "BDNF (excitatory)", 
                               'IGF1' = "IGF1 (inhibitory)", 
                               'VGF' = "VGF (shared)"))
ggplot(df_to_plot, aes(x = cell.type, 
                       y = logFC,
                       color = cell.type,
                       fill = time,
                       #group = time,
                       ymax = logFC -  SE, 
                       ymin = logFC + SE
)) + 
  xlab("") +
  ylab("log2(FC)") +
  geom_col(position = position_dodge(0.6),
           width = 0.5,
           color = "black") +
  geom_errorbar(color = "black",
                position = position_dodge(0.6),
                width = 0.5) +
  facet_grid(#labeller = gene_labeller,
             cols = vars(genes)
  ) +
  scale_fill_manual(values = c("indianred", "steelblue")) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) #+
  ggtitle("differential expression of early response genes")

