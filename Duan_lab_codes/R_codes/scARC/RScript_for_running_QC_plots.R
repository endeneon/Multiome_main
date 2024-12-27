# Chuxuan Li 05/15/2023
# RScript for 20230515_QC_plots_for_018-030_data.R part 2 

# init ####
library(Seurat)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
library(readr)

load("/nvmefs/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/018-030_RNA_integrated_labeled_with_harmony.RData")

lines <- sort(unique(integrated_labeled$cell.line.ident))

time_type_sum <- vector(mode = "list", length = 3L)
numlines <- length(lines)
types <- as.vector(unique(integrated_labeled$fine.cell.type))
numtypes <- length(types)
times <- unique(integrated_labeled$time.ident)

for (k in 1:length(times)){
  print(times[k])
  df <- data.frame(cell.type = rep_len(NA, numlines*numtypes),
                   cell.line = rep_len(NA, numlines*numtypes),
                   counts = rep_len(0, numlines*numtypes))
  for (i in 1:length(lines)){
    subobj <- subset(integrated_labeled, 
                     subset = cell.line.ident == lines[i] & time.ident == times[k])
    for (j in 1:numtypes){
      print(numtypes * i - numtypes + j)
      df$cell.type[numtypes * i - numtypes + j] <- types[j]
      df$cell.line[numtypes * i - numtypes + j] <- lines[i]
      df$counts[numtypes * i - numtypes + j] <- sum(subobj$fine.cell.type == types[j])
    }
  }
  df$cell.type <- factor(df$cell.type, levels = c("GABA", "SEMA3E+ GABA", "SST+ GABA", "CUX2+ GABA",
                                                  "NEFM- glut", "NEFM+ glut", "glut?",
                                                  "immature neuron", "VIM+ cells", "unknown"))
  time_type_sum[[k]] <- df
}

for (i in 1:3){
  df_to_plot <- time_type_sum[[i]]
  png(filename = paste0("./QC_plots/full_cellular_compos_", times[i], ".png"),
      width = 750, height = 750)
  p <- ggplot(data = df_to_plot,
              aes(x = cell.line, 
                  y = counts, 
                  fill = cell.type)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("#B33E52", #GABA
                                 "#CC7A88", #SEMA3E GABA
                                 "#E6B8BF", #SST GABA
                                 "#f5bfc8", #CUX2 GABA
                                 "#E6D2B8", #nmglut,
                                 "#CCAA7A", #npglut
                                 "#ad6e15", #glut?
                                 "#0075DC", #immature neuron
                                 "#11800d", #VIM+ cells
                                 "#6200e3" #unknown
    )) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) + 
    xlab("Cell Line") +
    ylab("Cellular composition (%)") +
    labs(fill = "Cell Types") +
    ggtitle(times[i])
  print(p)
  dev.off()
}

# General cell types (no specific GABA)
time_type_sum <- vector(mode = "list", length = 3L)
numlines <- length(lines)
types <- as.vector(unique(integrated_labeled$cell.type))
numtypes <- length(types)
times <- unique(integrated_labeled$time.ident)
for (k in 1:length(times)){
  print(times[k])
  df <- data.frame(cell.type = rep_len(NA, numlines*numtypes),
                   cell.line = rep_len(NA, numlines*numtypes),
                   counts = rep_len(0, numlines*numtypes))
  for (i in 1:length(lines)){
    subobj <- subset(integrated_labeled, 
                     subset = cell.line.ident == lines[i] & time.ident == times[k])
    for (j in 1:numtypes){
      print(types[j])
      print(numtypes * i - numtypes + j)
      df$cell.type[numtypes * i - numtypes + j] <- types[j]
      df$cell.line[numtypes * i - numtypes + j] <- lines[i]
      df$counts[numtypes * i - numtypes + j] <- sum(subobj$cell.type == types[j])
    }
  }
  df$cell.type <- factor(df$cell.type, levels = c("GABA", "nmglut", "npglut", "glut?",
                                                  "unidentified"))
  time_type_sum[[k]] <- df
}

for (i in 1:3){
  df_to_plot <- time_type_sum[[i]]
  png(filename = paste0("./QC_plots/general_cellular_compos_", times[i], ".png"),
      width = 750, height = 750)
  p <- ggplot(data = df_to_plot,
              aes(x = cell.line, 
                  y = counts, 
                  fill = cell.type)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("#B33E52", #GABA
                                 "#E6D2B8", #nmglut
                                 "#CCAA7A", #npglut
                                 "#ad6e15", #glut?
                                 "#4293db" #unknown
    )) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) + 
    xlab("Cell Line") +
    ylab("Cellular composition (%)") +
    labs(fill = "Cell Types") +
    ggtitle(times[i])
  print(p)
  dev.off()
}

