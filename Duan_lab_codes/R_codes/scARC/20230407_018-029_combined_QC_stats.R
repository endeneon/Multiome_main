# Chuxuan Li 04/07/2023
# Summarize QC stats for 018-029 combined data

# init ####
library(Seurat)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
library(readr)

load("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/029_RNA_integrated_labeled.RData")
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
  df$cell.type <- factor(df$cell.type, levels = c("GABA", "SEMA3E+ GABA",
                                                  "NEFM+ glut", "NEFM- glut", "unknown glut",
                                                  "stem cell-like", "unknown"))
  time_type_sum[[k]] <- df
}

for (i in 1:3){
  df_to_plot <- time_type_sum[[i]]
  png(filename = paste0("./QC/QC_plots/full_cellular_compos_", times[i], ".png"),
      width = 750, height = 750)
  p <- ggplot(data = df_to_plot,
              aes(x = cell.line, 
                  y = counts, 
                  fill = cell.type)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("#B33E52", #GABA
                                 "#CC7A88", #SEMA3E GABA
                                 "#E6D2B8", #nmglut
                                 "#CCAA7A", #npglut
                                 "#f29116", #unknown glut
                                 "#4C005C", #stem cell-like
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
  df$cell.type <- factor(df$cell.type, levels = c("GABA", "nmglut", "npglut", 
                                                  "unidentified"))
  time_type_sum[[k]] <- df
}

for (i in 1:3){
  df_to_plot <- time_type_sum[[i]]
  png(filename = paste0("./QC/QC_plots/general_cellular_compos_", times[i], ".png"),
      width = 750, height = 750)
  p <- ggplot(data = df_to_plot,
              aes(x = cell.line, 
                  y = counts, 
                  fill = cell.type)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("#B33E52", #GABA
                                 "#E6D2B8", #nmglut
                                 "#CCAA7A", #npglut
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


# 5. response gene expression by line ####
integrated_labeled$linextime.ident <- ""
for (l in lines) {
  for (t in times) {
    integrated_labeled$linextime.ident[integrated_labeled$cell.line.ident == l &
                                         integrated_labeled$time.ident == t] <- 
      paste(l, t, sep = " ")
  }
}
unique(integrated_labeled$linextime.ident)
DotPlot(integrated_labeled, assay = "SCT", cols = c("white", "red3"),
        features = c("EGR1", "NR4A1", "NPAS4", "FOS"), group.by = "linextime.ident") +
  coord_flip() +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10))

DotPlot(integrated_labeled, assay = "SCT", cols = c("white", "red3"),
        features = c("IGF1", "VGF", "BDNF"), group.by = "linextime.ident") +
  coord_flip() +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10))

# 6. cellular composition with only GABA and glut ####
time_type_sum <- vector(mode = "list", length = 3L)
numlines <- length(lines)
integrated_labeled <- subset(integrated_labeled, cell.type != "unidentified")
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
  df$cell.type <- factor(df$cell.type, levels = c("GABA", "SEMA3E+ GABA",
                                                  "NEFM+ glut", "NEFM- glut"))
  time_type_sum[[k]] <- df
}

for (i in 1:3){
  df_to_plot <- time_type_sum[[i]]
  png(filename = paste0("./QC/QC_plots/cellular_compos_after_filtering_", times[i], ".png"),
      width = 750, height = 750)
  p <- ggplot(data = df_to_plot,
              aes(x = cell.line, 
                  y = counts, 
                  fill = cell.type)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("#B33E52", #GABA
                                 "#CC7A88", #SEMA3E GABA
                                 "#E6D2B8", #nmglut
                                 "#CCAA7A" #npglut
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
for (i in 1:3){
  write.table(time_type_sum[[i]], 
              file = paste0("./QC/cell_counts_filtered_cell_type_", times[i], ".csv"),
              quote = F, sep = ",", row.names = F,
              col.names = T)
}
