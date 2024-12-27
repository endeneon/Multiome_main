# Chuxuan Li 02/17/2022
# getting the cell type composition for each cell line on the 20line RNAseq data
#and raw counts and cell counts for each library, then plot them

# init ####
library(Seurat)
library(ggplot2)
library(stringr)
library(RColorBrewer)

obj <- subset(RNAseq_integrated_labeled, cell.line.ident %in% c("CD_54", "CD_27",
                                                                "CD_26", "CD_08",
                                                                "CD_25"))
# cell counts and raw counts for each cell line ####
lines <- unique(human_only$cell.line.ident)
cell_counts_reads_df <- data.frame(cell.line = rep_len(NA, length(lines)),
                                   cell.counts = rep_len(0, length(lines)),
                                   n.reads = rep_len(0, length(lines)),
                                   reads.per.cell = rep_len(0, length(lines)))

for (i in 1:length(lines)){
  print(lines[i])
  cell_counts_reads_df$cell.line <- lines[i]
  cell_counts_reads_df$cell.counts[i] <- sum(filtered_obj$cell.line.ident == lines[i])
  temp <- subset(filtered_obj, cell.line.ident == lines[i])
  cell_counts_reads_df$n.reads[i] <- sum(temp@assays$RNA@counts@x)
  cell_counts_reads_df$reads.per.cell[i] <- sum(temp@assays$RNA@counts@x)/cell_counts_reads_df$cell.counts[i]
}

write.table(cell_counts_reads_df, file = "renormalized_cell_counts_UMI_counts.csv",
            quote = F, sep = ",", row.names = F,
            col.names = T)

# plot the cell count first
p <- ggplot(cell_counts_reads_df,
            aes(x = lines,
                y = rep_len(2, length(lines)),
                fill = cell.counts)) +
  geom_col(width = 1, 
           position = position_dodge(width = 0.5)) +
  geom_text(aes(label = cell.counts), 
            nudge_y = -1,
            hjust = "middle") + 
  scale_fill_continuous(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) + 
  ggtitle("# Nuclei") +
  coord_flip() +
  NoLegend()

p

# plot UMI per cell next
p <- ggplot(cell_counts_reads_df,
            aes(x = lines,
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

p

# plot UMI total for each cell line
p <- ggplot(cell_counts_reads_df,
            aes(x = lines,
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

p

# cell counts and raw counts by library ####
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/normalize_by_6libs_integrated.RData")
libnames <- c("g_2_0", "g_2_1", "g_2_6", "g_8_0", "g_8_1", "g_8_6")
libs <- (unique(integrated$lib.ident))
cell_counts_reads_df <- data.frame(library = rep_len(NA, length(libs)),
                                   cell.counts = rep_len(0, length(lines)),
                                   n.reads = rep_len(0, length(lines)),
                                   reads.per.cell = rep_len(0, length(lines)))

for (i in 1:length(libs)){
  print(libnames[i])
  cell_counts_reads_df$library <- libnames[i]
  cell_counts_reads_df$cell.counts[i] <- sum(filtered_obj$lib.ident == libs[i])
  temp <- subset(integrated, lib.ident == libs[i])
  cell_counts_reads_df$n.reads[i] <- sum(temp@assays$RNA@counts@x)
  cell_counts_reads_df$reads.per.cell[i] <- sum(temp@assays$RNA@counts@x)/cell_counts_reads_df$cell.counts[i]
}

write.table(cell_counts_reads_df, file = "renormalized_cell_counts_UMI_counts.csv",
            quote = F, sep = ",", row.names = F,
            col.names = T)
cell_counts_reads_df$library <- factor(cell_counts_reads_df$library, 
                                       levels = rev(libnames))
# plot the cell count first
p <- ggplot(cell_counts_reads_df,
            aes(x = library,
                y = rep_len(2, length(libnames)),
                fill = cell.counts)) +
  geom_col(width = 1, 
           position = position_dodge(width = 0.5)) +
  geom_text(aes(label = cell.counts), 
            nudge_y = -1,
            hjust = "middle") + 
  scale_fill_continuous(low = "lightskyblue1", high = "lightskyblue3") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()) + 
  ggtitle("# Nuclei") +
  coord_flip() +
  NoLegend()

p

# plot UMI per cell next
p <- ggplot(cell_counts_reads_df,
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

p

# plot UMI total for each cell line
p <- ggplot(cell_counts_reads_df,
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

p

# cell type and time point composition ####
time_type_sum <- vector(mode = "list", length = 3L)
filtered_obj$cell.type <- factor(filtered_obj$cell.type, levels = c("NEFM_pos_glut",
                                                                    "NEFM_neg_glut",
                                                                    "GABA", "NPC"))
types <- as.vector(unique(filtered_obj$cell.type))
numtypes <- length(types)
times <- unique(filtered_obj$time.ident)

for (k in 1:length(times)){
  print(times[k])
  df <- data.frame(cell.type = rep_len(NA, 5*numtypes),
                   cell.line = rep_len(NA, 5*numtypes),
                   counts = rep_len(0, 5*numtypes))
  for (i in 1:length(lines)){
    subobj <- subset(filtered_obj, 
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
  jpeg(filename = paste0("full_cellular_compos_", times[i], ".jpeg"))
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
