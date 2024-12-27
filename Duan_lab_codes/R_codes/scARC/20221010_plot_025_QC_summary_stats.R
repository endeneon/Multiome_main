# Chuxuan Li 11/10/2022
# Summarize QC stats (cell counts before and after rat identification, demux, QC)

# init ####
{
  library(Seurat)
  library(RColorBrewer)
  library(viridis)
  library(ggplot2)
  library(stringr)
}

plan("multisession", workers = 4)
set.seed(42)
options(future.globals.maxSize = 229496729600)




load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/GRCh38_mapped_removed_rat_cells_assigned_demux_bc_list.RData")
h5list <- sort(list.files(path = "/nvmefs/scARC_Duan_025_GRCh38",
                          pattern = "filtered_feature_bc_matrix.h5",
                          recursive = T,
                          include.dirs = F, 
                          full.names = T))
library(readr)
summary <- 
  read_csv("Analysis_part2_GRCh38/cell_counts_based_solely_on_GRCh38_mapped_data.csv")

# 1. proportions of rat, unmapped, mapped cells ####
cell_counts_df <- data.frame(library = rep(str_extract(string = h5list,
                                                       pattern = "[0-9]+-[0-6]"),
                                           each = 3),
                             count = rep_len(0, 3 * length(objlist)),
                             type_of_count = rep(c("rat astrocytes", 
                                                   "human unmapped to demuxed barcodes",
                                                   "human cells"),
                                                 times = length(objlist)))
cell_counts_df$type_of_count <- factor(cell_counts_df$type_of_count, 
                                       levels = c("rat astrocytes", 
                                                  "human unmapped to demuxed barcodes",
                                                  "human cells"))
cell_counts_df$library <- factor(cell_counts_df$library,
                                 levels = (unique(cell_counts_df$library)))
for (i in 1:length(objlist)){
  cell_counts_df$count[i * 3] <- summary$human[summary$library == cell_counts_df$library[i * 3]]
  cell_counts_df$count[i * 3 - 1] <- 
    summary$nonhuman[summary$library == cell_counts_df$library[i * 3]]
  cell_counts_df$count[i * 3 - 2] <- summary$rat[summary$library == cell_counts_df$library[i * 3]]
}

ggplot(data = cell_counts_df,
       aes(x = library, 
           y = count, 
           fill = type_of_count)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = rev(brewer.pal(3, "Set3"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) + 
  xlab("Library") +
  ylab("Cell counts") +
  labs(fill = "Cell Types") +
  ggtitle("Number of rat astrocytes, unmapped cells,\nand human cells per library")

ggplot(data = cell_counts_df,
       aes(x = library, 
           y = count, 
           fill = type_of_count)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = rev(brewer.pal(3, "Set3"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) + 
  xlab("Library") +
  ylab("Percentage of ") +
  labs(fill = "Cell Types") +
  ggtitle("Percentage of rat astrocytes, unmapped cells,\nand human cells per library")

# 2. summarize number of cells and counts per library ####
plotCellCountsRectangles <- function(df) {
  p <- ggplot(df, aes(x = library, y = rep_len(2, length(libnames)), fill = cell.counts)) +
    geom_col(width = 1, position = position_dodge(width = 0.5)) +
    geom_text(aes(label = cell.counts), nudge_y = -1, hjust = "middle") + 
    scale_fill_continuous(low = "#a2c1f2", high = "#5c96f2") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          axis.title = element_blank(), axis.text.x = element_blank(),
          axis.ticks = element_blank()) + 
    ggtitle("# Nuclei") +
    coord_flip() +
    NoLegend()
  print(p)
}
plotUMIBars <- function(df) {
  p <- ggplot(df, aes(x = library, y = reads.per.cell)) +
    ylab("# UMIs per nuclei") +
    geom_col(position = "dodge", color = "gray27", fill = "gray90") +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 8)) + 
    coord_flip()
  print(p)
}

plotTotalUMIBars <- function(df) {
  p <- ggplot(df, aes(x = library, y = n.reads)) +
    ylab("# total UMIs") +
    geom_col(position = "dodge", color = "gray27", fill = "gray90") +
    theme_minimal() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 8)) + 
    coord_flip()
  print(p)
}

libnames <- str_extract(string = h5list,
                        pattern = "[0-9]+-[0-6]")
df_raw <- data.frame(library = libnames,
                     cell.counts = rep_len(0, length(libnames)),
                     n.reads = rep_len(0, length(libnames)),
                     reads.per.cell = rep_len(0, length(libnames)))
df_qc <- data.frame(library = libnames,
                     cell.counts = rep_len(0, length(libnames)),
                     n.reads = rep_len(0, length(libnames)),
                     reads.per.cell = rep_len(0, length(libnames)))
df_human <- data.frame(library = libnames,
                       cell.counts = rep_len(0, length(libnames)),
                       n.reads = rep_len(0, length(libnames)),
                       reads.per.cell = rep_len(0, length(libnames)))
load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/GRCh38_mapped_raw_list.RData")
load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/GRCh38_mapped_demux_matched_human_only_list.RData")
load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/GRCh38_mapped_after_QC_list.RData")
# raw
for (i in 1:length(libnames)){
  print(libnames[i])
  df_raw$cell.counts[i] <- length(objlist[[i]]@assays$RNA@counts@Dimnames[[2]])
  df_raw$n.reads[i] <- sum(objlist[[i]]@assays$RNA@counts@x)
  df_raw$reads.per.cell[i] <- 
    sum(objlist[[i]]@assays$RNA@counts@x)/df_raw$cell.counts[i]
}
write.table(df_raw, 
            file = "./Analysis_part2_GRCh38/QC/raw_cell_counts_UMI_counts_per_lib.csv",
            quote = F, sep = ",", row.names = F,
            col.names = T)
# human only
for (i in 1:length(libnames)){
  print(libnames[i])
  df_human$cell.counts[i] <- length(human_only_lst[[i]]@assays$RNA@counts@Dimnames[[2]])
  df_human$n.reads[i] <- sum(human_only_lst[[i]]@assays$RNA@counts@x)
  df_human$reads.per.cell[i] <- 
    sum(human_only_lst[[i]]@assays$RNA@counts@x)/df_human$cell.counts[i]
}
write.table(df_human, 
            file = "./Analysis_part2_GRCh38/QC/human_cell_counts_UMI_counts_per_lib.csv",
            quote = F, sep = ",", row.names = F,
            col.names = T)
# after QC
for (i in 1:length(libnames)){
  print(libnames[i])
  df_qc$cell.counts[i] <- length(QCed_lst[[i]]@assays$RNA@counts@Dimnames[[2]])
  df_qc$n.reads[i] <- sum(QCed_lst[[i]]@assays$RNA@counts@x)
  df_qc$reads.per.cell[i] <- 
    sum(QCed_lst[[i]]@assays$RNA@counts@x)/df_qc$cell.counts[i]
}
write.table(df_qc, 
            file = "./Analysis_part2_GRCh38/QC/QC_cell_counts_UMI_counts_per_lib.csv",
            quote = F, sep = ",", row.names = F,
            col.names = T)

df_raw$library <- factor(df_raw$library, levels = rev(libnames))
plotCellCountsRectangles(df_qc)
plotUMIBars(df_qc)
plotTotalUMIBars(df_qc)

# 3. number of cells and counts per line ####
load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/integrated_obj_nfeature_8000.RData")
lines <- sort(unique(integrated$cell.line.ident))
cell_counts_reads_df <- data.frame(line = lines,
                                   cell.counts = rep_len(0, length(lines)),
                                   n.reads = rep_len(0, length(lines)),
                                   reads.per.cell = rep_len(0, length(lines)))

for (i in 1:length(lines)){
  print(lines[i])
  tempobj <- subset(integrated, cell.line.ident == lines[i])
  cell_counts_reads_df$cell.counts[i] <- length(tempobj@assays$RNA@counts@Dimnames[[2]])
  cell_counts_reads_df$n.reads[i] <- sum(tempobj@assays$RNA@counts@x)
  cell_counts_reads_df$reads.per.cell[i] <- 
    sum(tempobj@assays$RNA@counts@x)/cell_counts_reads_df$cell.counts[i]
}

write.table(cell_counts_reads_df, 
            file = "./Analysis_part2_GRCh38/QC/cell_counts_UMI_counts_per_line.csv",
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

# 4. cellular composition ####
load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/integrated_labeled.RData")
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
  df$cell.type <- factor(df$cell.type, levels = c("GABA", "SEMA3E+ GABA", "SST+ GABA",
                                                  "NEFM+ glut", "NEFM- glut", 
                                                  "immature neuron", "unknown neuron", "unknown"))
  time_type_sum[[k]] <- df
}

for (i in 1:3){
  df_to_plot <- time_type_sum[[i]]
  png(filename = paste0("./Analysis_part2_GRCh38/QC/full_cellular_compos_", times[i], ".png"),
      width = 560, height = 560)
  p <- ggplot(data = df_to_plot,
              aes(x = cell.line, 
                  y = counts, 
                  fill = cell.type)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("#B33E52", #GABA
                                 "#CC7A88", #SEMA3E GABA
                                 "#E6B8BF", #SST GABA
                                 "#E6D2B8", #nmglut
                                 "#CCAA7A", #npglut
                                 "#0075DC", #immature neuron
                                 "#4C005C", #unknown neuron
                                 "#993F00" #unknown
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
  png(filename = paste0("./Analysis_part2_GRCh38/QC/general_cellular_compos_", times[i], ".png"),
      width = 560, height = 560)
  p <- ggplot(data = df_to_plot,
              aes(x = cell.line, 
                  y = counts, 
                  fill = cell.type)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("#B33E52", #GABA
                                 "#E6D2B8", #nmglut
                                 "#CCAA7A", #npglut
                                 "#0075DC", #immature neuron
                                 "#4C005C", #unknown neuron
                                 "#993F00" #unknown
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

lines <- unique(integrated_labeled$cell.line.ident)
times <- unique(integrated_labeled$time.ident)
# cell.line.ident == lines[i] & time.ident == times[k]

for (l in lines) {
  for (t in times) {
    integrated_labeled$linextime.ident[integrated_labeled$cell.line.ident == l &
                                         integrated_labeled$time.ident == t] <- 
      paste(l, t, sep = " ")
  }
}

sort(unique(integrated_labeled$linextime.ident))


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

DotPlot(integrated_labeled, 
        assay = "SCT", 
        cols = c("white", 
                 "red3"),
        features = c("BDNF", 
                     "VGF", 
                     "FOS", 
                     "NPAS4"), 
        group.by = "linextime.ident") +
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
  df$cell.type <- factor(df$cell.type, levels = c("GABA", "SEMA3E+ GABA", "SST+ GABA",
                                                  "NEFM+ glut", "NEFM- glut"))
  time_type_sum[[k]] <- df
}

for (i in 1:3){
  df_to_plot <- time_type_sum[[i]]
  png(filename = paste0("./Analysis_part2_GRCh38/QC/cellular_compos_after_filtering_", times[i], ".png"),
      width = 560, height = 560)
  p <- ggplot(data = df_to_plot,
              aes(x = cell.line, 
                  y = counts, 
                  fill = cell.type)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c("#B33E52", #GABA
                                 "#CC7A88", #SEMA3E GABA
                                 "#E6B8BF", #SST GABA
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
              file = paste0("./Analysis_part2_GRCh38/QC/cell_counts_filtered_cell_type_", times[i], ".csv"),
              quote = F, sep = ",", row.names = F,
              col.names = T)
}
