library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(Signac)

load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")
load("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/multiomic_obj.RData")
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000_with_covars_from_30Sep2022_single_cell_DE.RData")

# RNAseq dimplot ####
obj$cell.type.forplot <- str_replace_all(str_replace_all(str_replace_all(obj$cell.type, "_", " "), " pos", "+"), " neg", "-")
DimPlot(obj, label.size = 4, group.by = "cell.type.forplot", 
        cols = rev(brewer.pal(6, "Set3")), label = T) +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10)) +
  NoLegend() +
  ggtitle("Gene expression")

# RNAseq cell type markers feature plot ####
plots <- FeaturePlot(obj,
            features = c("MAP2", "GAD1", "GAD2", "SLC17A6", "SLC17A7", "NEFM"),
            ncol = 3, cols = c("grey", "royalblue4"), combine = F)
plots <- lapply(X = plots, FUN = function(x) x+theme(text = element_text(size = 10), 
                                                     axis.text = element_text(size = 10),
                                                     axis.title = element_text(size = 10),
                                                     legend.text = element_text(size = 10)))
CombinePlots(plots = plots)

# ATACseq dimplot ####
multiomic_obj$cell.type.forplot <- 
  str_replace_all(str_replace_all(str_replace_all(multiomic_obj$cell.type, "_", " "), " pos", "+"), " neg", "-")
DimPlot(multiomic_obj, label.size = 4, group.by = "cell.type.forplot", 
        cols = rev(brewer.pal(6, "Set3")), label = T) +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10)) +
  NoLegend() +
  ggtitle("Chromatin accessibility")

# cellular composition ####
lines <- sort(unique(obj$cell.line.ident))
time_type_sum <- vector(mode = "list", length = 3L)
numlines <- length(lines)
types <- as.vector(unique(obj$cell.type.forplot))
numtypes <- length(types)
times <- unique(obj$time.ident)

for (k in 1:length(times)){
  print(times[k])
  df <- data.frame(cell.type = rep_len(NA, numlines*numtypes),
                   cell.line = rep_len(NA, numlines*numtypes),
                   counts = rep_len(0, numlines*numtypes))
  for (i in 1:length(lines)){
    subobj <- subset(obj, 
                     subset = cell.line.ident == lines[i] & time.ident == times[k])
    for (j in 1:numtypes){
      print(numtypes * i - numtypes + j)
      df$cell.type[numtypes * i - numtypes + j] <- types[j]
      df$cell.line[numtypes * i - numtypes + j] <- lines[i]
      df$counts[numtypes * i - numtypes + j] <- sum(subobj$cell.type.forplot == types[j])
    }
  }
  time_type_sum[[k]] <- df
}

for (i in 1:3){
  df_to_plot <- time_type_sum[[i]]
  pdf(paste0("../../ASHG 2022/cellular_compos_", times[i], ".pdf"))
  p <- ggplot(data = df_to_plot,
              aes(x = cell.line, 
                  y = counts, 
                  fill = cell.type)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = (brewer.pal(11, "Set3"))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) + 
    xlab("Cell Line") +
    ylab("Cellular composition (%)") +
    labs(fill = "Cell Types") +
    ggtitle(times[i])
  print(p)
  dev.off()
}

# volcano plots for DESeq2 DE ####
library(ggrepel)
library(stringr)
library(readr)

pathlist <- list.files("./pseudobulk_DE/res_use_combatseq_mat/filtered_by_basemean_only/", full.names = T,
                       pattern = "1v0|6v0")
pathlist
reslist <- vector('list', length(pathlist))
for (i in 1:length(reslist)){
  reslist[[i]] <- read_csv(pathlist[i])
}
# make plots ####
namelist <- str_extract(pathlist, "[A-Za-z]+_[1|6]v0")
namelist
setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/pseudobulk_DE/res_use_combatseq_mat/volcano_plots")

for (i in 1:length(reslist)){
  print(namelist[i])
  res <- reslist[[i]]
  #res$gene.symbol <- rownames(res)
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
  pdf(paste0(namelist[i], "_volcano_plot.pdf"))
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

# response gene DE result dotplot ####
early_genes <- c("FOS", "FOSB", "NPAS4", "NR4A1")
late_genes <- c("BDNF", "VGF", "IGF1")
for (i in 1:length(reslist)) {
  df_early <- reslist[[i]][reslist[[i]]$gene %in% early_genes, ]
  df_late <- reslist[[i]][reslist[[i]]$gene %in% late_genes, ]
  
  df_early$neg_log_pval <- (-1) * log10(df_early$pvalue)
  df_early$time <- str_extract(namelist[i], "[1|6]v0")
  df_early$cell.type <- str_extract(namelist[i], "[A-Za-z]+")
  
  df_late$neg_log_pval <- (-1) * log10(df_late$pvalue)
  df_late$time <- str_extract(namelist[i], "[1|6]v0")
  df_late$cell.type <- str_extract(namelist[i], "[A-Za-z]+")
  if (df_early$padj > 0.05) {
    df_early$log2FoldChange <- NA
    df_early$neg_log_pval <- NA
  }
  if (df_late$padj > 0.05) {
    df_late$log2FoldChange <- NA
    df_late$neg_log_pval <- NA
  }
  if (i == 1) {
    df_to_plot_early <- df_early
    df_to_plot_late <- df_late
  } else {
    df_to_plot_early <- rbind(df_to_plot_early, df_early)
    df_to_plot_late <- rbind(df_to_plot_late, df_late)
  }
}
df_to_plot_early$neg_log_pval[is.infinite(df_to_plot_early$neg_log_pval)] <- 500
df_to_plot_late$neg_log_pval[is.infinite(df_to_plot_late$neg_log_pval)] <- 500

ggplot(df_to_plot_early,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = neg_log_pval,
           fill = log2FoldChange)) +
  geom_point(shape = 21) +
  ylab("Cell Type") +
  xlab("") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  facet_grid(cols = vars(gene)) + 
  theme_bw() +
  labs(size = "-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 

ggplot(df_to_plot_late,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = neg_log_pval,
           fill = log2FoldChange)) +
  geom_point(shape = 21) +
  ylab("Cell Type") +
  xlab("") +
  scale_fill_gradientn(colours = c("white", "royalblue4")) +
  facet_grid(cols = vars(gene)) + 
  theme_bw() +
  labs(size = "-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
