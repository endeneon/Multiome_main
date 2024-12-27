# Chuxuan Li 05/15/2023
# generate QC statistics (cell counts by line after QC, cellular composition) for
#018-030 integrated data

# init ####
library(Seurat)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
library(readr)

load("/nvmefs/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/018-030_RNA_integrated_labeled_with_harmony.RData")

# 1. number of cells and counts per line ####
lines <- sort(unique(integrated_labeled$cell.line.ident))
cell_counts_reads_df <- data.frame(line = lines,
                                   cell.counts = rep_len(0, length(lines)),
                                   n.reads = rep_len(0, length(lines)),
                                   reads.per.cell = rep_len(0, length(lines)))

for (i in 1:length(lines)){
  print(lines[i])
  tempobj <- subset(integrated_labeled, cell.line.ident == lines[i])
  cell_counts_reads_df$cell.counts[i] <- length(tempobj@assays$RNA@counts@Dimnames[[2]])
  cell_counts_reads_df$n.reads[i] <- sum(tempobj@assays$RNA@counts@x)
  cell_counts_reads_df$reads.per.cell[i] <- 
    sum(tempobj@assays$RNA@counts@x)/cell_counts_reads_df$cell.counts[i]
}

write.table(cell_counts_reads_df, 
            file = "./QC/cell_counts_UMI_counts_per_line_after_QC_in_018_030_combined_data.csv",
            quote = F, sep = ",", row.names = F,
            col.names = T)
cell_counts_reads_df$line <- factor(cell_counts_reads_df$line, 
                                    levels = rev(lines))
# plot the cell count first
ggplot(cell_counts_reads_df,
       aes(x = line,
           y = rep_len(2, length(lines)),
           fill = cell.counts)) +
  geom_col(width = 1, 
           position = position_dodge(width = 0.5)) +
  geom_text(aes(label = cell.counts), 
            nudge_y = -1,
            hjust = "middle") + 
  scale_fill_continuous(low = "#bee8c2", high = "#7bb881") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) + 
  ggtitle("# Nuclei") +
  coord_flip() +
  NoLegend()

# plot UMI per cell next
ggplot(cell_counts_reads_df,
       aes(x = line,
           y = reads.per.cell)) +
  ylab("# UMIs per nuclei") +
  geom_col(position = "dodge",
           color = "gray27",
           fill = "gray90") +
  #geom_text(aes(label = reads.per.cell)) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 11)) + 
  coord_flip()

# plot UMI total for each cell line
ggplot(cell_counts_reads_df,
       aes(x = line,
           y = n.reads)) +
  ylab("# total UMIs") +
  geom_col(position = "dodge",
           color = "gray27",
           fill = "gray90") +
  #geom_text(aes(label = reads.per.cell)) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 11),
        axis.text.x = element_text(size = 8)) + 
  coord_flip()

# 2. cellular composition ####
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

# 3. response gene expression by line ####
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
