# Chuxuan Li 03/24/2022
# Get QC stats and make graphs for the group05-31-33-51-63 batch of 20-line data

# init ####
library(Seurat)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
setwd("/data/FASTQ/Duan_Project_024/hybrid_output/")
# note this h5 file contains both atac-seq and gex information
# Read files, separate the ATACseq data from the .h5 matrix
h5list <- list.files(path = ".",
                     pattern = "filtered_feature_bc_matrix.h5",
                     recursive = T,
                     include.dirs = F, 
                     full.names = T)

setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate")
# cell counts: total, rat, human unmapped ####
load("raw_data_list.RData")
load("removed_rat_list.RData")

# count number of cells unmatched for each lib, build data frame for 3 counts
# raw number of cells, rat astrocytes, unmatched to barcodes
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
for (i in 1:length(cleanobj_lst)){
  cell_counts_df$count[i * 3] <- sum(cleanobj_lst[[i]]$cell.line.ident != "unmatched")
  cell_counts_df$count[i * 3 - 2] <- 
    length(objlist[[i]]@assays$RNA@counts@Dimnames[[2]]) - length(cleanobj_lst[[i]]$rat.ident)
  cell_counts_df$count[i * 3 - 1] <- sum(cleanobj_lst[[i]]$cell.line.ident == "unmatched")
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


# raw cell counts and raw read by library ####
libnames <- str_extract(string = h5list,
                        pattern = "[0-9]+-[0-6]")
cell_counts_reads_df <- data.frame(library = libnames,
                                   cell.counts = rep_len(0, length(libnames)),
                                   n.reads = rep_len(0, length(libnames)),
                                   reads.per.cell = rep_len(0, length(libnames)))

for (i in 1:length(libnames)){
  print(libnames[i])
  cell_counts_reads_df$cell.counts[i] <- length(objlist[[i]]@assays$RNA@counts@Dimnames[[2]])
  cell_counts_reads_df$n.reads[i] <- sum(objlist[[i]]@assays$RNA@counts@x)
  cell_counts_reads_df$reads.per.cell[i] <- 
    sum(objlist[[i]]@assays$RNA@counts@x)/cell_counts_reads_df$cell.counts[i]
}

write.table(cell_counts_reads_df, 
            file = "./QC_dataframes/total_unfiltered_cell_counts_UMI_counts_per_lib.csv",
            quote = F, sep = ",", row.names = F,
            col.names = T)
cell_counts_reads_df$library <- factor(cell_counts_reads_df$library, 
                                       levels = rev(libnames))
# plot the cell count first
ggplot(cell_counts_reads_df,
       aes(x = library,
           y = rep_len(2, length(libnames)),
           fill = cell.counts)) +
  geom_col(width = 1, 
           position = position_dodge(width = 0.5)) +
  geom_text(aes(label = cell.counts), 
            nudge_y = -1,
            hjust = "middle") + 
  scale_fill_continuous(low = "#a2c1f2", high = "#5c96f2") +
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
            aes(x = library,
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
            aes(x = library,
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


# human only cell counts and raw read by library ####
libnames <- str_extract(string = h5list,
                        pattern = "[0-9]+-[0-6]")
cell_counts_reads_df <- data.frame(library = libnames,
                                   cell.counts = rep_len(0, length(libnames)),
                                   n.reads = rep_len(0, length(libnames)),
                                   reads.per.cell = rep_len(0, length(libnames)))

for (i in 1:length(libnames)){
  print(libnames[i])
  cell_counts_reads_df$cell.counts[i] <- length(cleanobj_lst[[i]]@assays$RNA@counts@Dimnames[[2]])
  cell_counts_reads_df$n.reads[i] <- sum(cleanobj_lst[[i]]@assays$RNA@counts@x)
  cell_counts_reads_df$reads.per.cell[i] <- 
    sum(cleanobj_lst[[i]]@assays$RNA@counts@x)/cell_counts_reads_df$cell.counts[i]
}

write.table(cell_counts_reads_df, 
            file = "./QC_dataframes/filtered_cell_counts_UMI_counts_per_lib.csv",
            quote = F, sep = ",", row.names = F,
            col.names = T)
cell_counts_reads_df$library <- factor(cell_counts_reads_df$library, 
                                       levels = rev(libnames))
# plot the cell count first
ggplot(cell_counts_reads_df,
       aes(x = library,
           y = rep_len(2, length(libnames)),
           fill = cell.counts)) +
  geom_col(width = 1, 
           position = position_dodge(width = 0.5)) +
  geom_text(aes(label = cell.counts), 
            nudge_y = -1,
            hjust = "middle") + 
  scale_fill_continuous(low = "#a2c1f2", high = "#5c96f2") +
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
            aes(x = library,
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
            aes(x = library,
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

# cell counts and raw read by library after QC filtering by pct.mt, nfeat, ncount ####
sum(integrated_labeled$percent.mt > 15 | integrated_labeled$nFeature_RNA > 8000 | integrated_labeled$nCount_RNA > 40000)
integrated_labeled <- subset(integrated_labeled, 
                             subset = nFeature_RNA > 400 & 
                               nFeature_RNA < 8000 & 
                               nCount_RNA > 500 &
                               nCount_RNA < 40000 &
                               percent.mt < 15)
libnames <- sort(unique(integrated_labeled$orig.ident))
libnames <- c(libnames[7:9], libnames)
libnames <- c(libnames[1:9], libnames[13:18])
cell_counts_reads_df <- data.frame(library = libnames,
                                   cell.counts = rep_len(0, length(libnames)),
                                   n.reads = rep_len(0, length(libnames)),
                                   reads.per.cell = rep_len(0, length(libnames)))

for (i in 1:length(libnames)){
  print(libnames[i])
  tempobj <- subset(integrated_labeled, orig.ident == libnames[i])
  cell_counts_reads_df$cell.counts[i] <- length(tempobj@assays$RNA@counts@Dimnames[[2]])
  cell_counts_reads_df$n.reads[i] <- sum(tempobj@assays$RNA@counts@x)
  cell_counts_reads_df$reads.per.cell[i] <- 
    sum(tempobj@assays$RNA@counts@x)/cell_counts_reads_df$cell.counts[i]
}

cell_counts_reads_df$library <- factor(cell_counts_reads_df$library, 
                                       levels = rev(libnames))
# plot the cell count first
ggplot(cell_counts_reads_df,
       aes(x = library,
           y = rep_len(2, length(libnames)),
           fill = cell.counts)) +
  geom_col(width = 1, 
           position = position_dodge(width = 0.5)) +
  geom_text(aes(label = cell.counts), 
            nudge_y = -1,
            hjust = "middle") + 
  scale_fill_continuous(low = "#a2c1f2", high = "#5c96f2") +
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
       aes(x = library,
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
       aes(x = library,
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


# cell counts and raw counts for each cell line ####
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
            file = "./QC_dataframes/cell_counts_UMI_counts_per_line.csv",
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


# cell type and time point composition ####
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
      print(numtypes * i - numtypes + j)
      df$cell.type[numtypes * i - numtypes + j] <- types[j]
      df$cell.line[numtypes * i - numtypes + j] <- lines[i]
      df$counts[numtypes * i - numtypes + j] <- sum(subobj$cell.type == types[j])
    }
  }
  time_type_sum[[k]] <- df
}

for (i in 1:3){
  df_to_plot <- time_type_sum[[i]]
  jpeg(filename = paste0("./integrated_5000vargene_plots/full_cellular_compos_", times[i], ".jpeg"))
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

for (i in 1:3){
  write.table(time_type_sum[[i]], 
              file = paste0("./cell_counts_per_4_cell_type_", times[i], ".csv"),
              quote = F, sep = ",", row.names = F,
              col.names = T)
}

# plot time point composition in each library ####
df <- data.frame(time.point = rep_len(NA, 5*3),
                 cell.line = rep_len(NA, 5*3),
                 counts = rep_len(0, 5*3))

for (i in 1:length(lines)){
  for (j in 1:length(times)){
    print(3 * i - 3 + j)
    sub <- subset(human_only, 
                  subset = cell.line.ident == lines[i])
    df$time.point[3 * i - 3 + j] <- times[j]
    df$cell.line[3 * i - 3 + j] <- lines[i]
    df$counts[3 * i - 3 + j] <- sum(sub$time.ident == times[j])
  }
}
write.table(df, file = "time_point_compos.csv",
            quote = F, sep = ",", row.names = F,
            col.names = T)

p <- ggplot(data = df,
            aes(x = cell.line, 
                y = counts, 
                fill = time.point)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = brewer.pal(3, "Set2")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) + 
  xlab("Cell Line") +
  ylab("Time point composition (%)") +
  labs(fill = "Time Points")
print(p)

# expression per cell line and time point for response genes ####
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")
late <- c("BDNF", "VGF", "IGF1")
early <- c("FOS", "NPAS4", "NR4A1", "EGR1")
types <- unique(integrated_labeled$cell.type)
lines <- unique(integrated_labeled$cell.line.ident)
times <- unique(integrated_labeled$time.ident)
integrated_labeled$timextype.ident <- "NA"
integrated_labeled$timextype.ident.forplot <- "NA"

for (i in 1:length(types)){
  print(types[i])
  for (j in 1:length(times)){
    print(times[j])
    integrated_labeled$timextype.ident[integrated_labeled$time.ident == times[j] &
                                         integrated_labeled$cell.type == types[i]] <- 
      paste(types[i], times[j], sep = "_")
    integrated_labeled$timextype.ident.forplot[integrated_labeled$time.ident == times[j] &
                                         integrated_labeled$cell.type == types[i]] <- 
      paste(str_replace_all(types[i], "_", " "), times[j], sep = " ")
  }
}
unique(integrated_labeled$timextype.ident)
unique(integrated_labeled$timextype.ident.forplot)
integrated_labeled$timextype.ident.forplot <-
  str_replace(integrated_labeled$timextype.ident.forplot, " pos", "+")
integrated_labeled$timextype.ident.forplot <-
  str_replace(integrated_labeled$timextype.ident.forplot, " neg", "-")

integrated_labeled$timexline.ident <- "NA"
integrated_labeled$timexline.ident.forplot <- "NA"
for (i in 1:length(lines)){
  print(lines[i])
  for (j in 1:length(times)){
    print(times[j])
    integrated_labeled$timexline.ident[integrated_labeled$time.ident == times[j] &
                                         integrated_labeled$cell.line.ident == lines[i]] <- 
      paste(lines[i], times[j], sep = "_")
    integrated_labeled$timexline.ident.forplot[integrated_labeled$time.ident == times[j] &
                                                 integrated_labeled$cell.line.ident == lines[i]] <- 
      paste(lines[i], times[j], sep = " ")
  }
}
unique(integrated_labeled$timexline.ident)
unique(integrated_labeled$timexline.ident.forplot)

DPnew(integrated_labeled, 
      features = late,
      cols = rev(brewer.pal(n = 4, name = "Reds")),
      group.by = "timexline.ident.forplot") + 
  RotatedAxis() +
  coord_flip() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11))
DPnew(integrated_labeled, 
      features = early,
      cols = rev(brewer.pal(n = 4, name = "Reds")),
      group.by = "timexline.ident.forplot") + 
  RotatedAxis() +
  coord_flip() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11))
save(integrated_labeled, file = "labeled_nfeat3000.RData")
