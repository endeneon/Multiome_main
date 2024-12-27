# Chuxuan Li 03/29/2022
# Get QC stats from 18-line ATACseq data

# init ####
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)

library(ggplot2)
library(patchwork)

library(stringr)
library(future)

set.seed(2022)
plan("multisession", workers = 1)

# load data ####
pathlist <- list.files(path = "/data/FASTQ/Duan_Project_024/hybrid_output/",
                       pattern = "filtered_feature_bc_matrix.h5", full.names = T,
                       recursive = T)

load("./raw_obj_list.RData")
load("./no_rat_no_unmapped_no_rat_gene_clean_list.RData")

# 1. raw data: total cell counts and raw read by library ####
libnames <- str_extract(string = pathlist,
                        pattern = "[0-9]+-[0-6]")
cell_counts_reads_df <- data.frame(library = libnames,
                                   cell.counts = rep_len(0, length(libnames)),
                                   n.reads = rep_len(0, length(libnames)),
                                   reads.per.cell = rep_len(0, length(libnames)))

for (i in 1:length(libnames)){
  print(libnames[i])
  cell_counts_reads_df$cell.counts[i] <- length(objlist[[i]]@assays$ATAC@counts@Dimnames[[2]])
  cell_counts_reads_df$n.reads[i] <- sum(objlist[[i]]@assays$ATAC@counts@x)
  cell_counts_reads_df$reads.per.cell[i] <- 
    sum(objlist[[i]]@assays$ATAC@counts@x)/cell_counts_reads_df$cell.counts[i]
}

write.table(cell_counts_reads_df, 
            file = "./unfiltered_cell_counts_UMI_counts_per_lib.csv",
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


# 2. human only: cell counts and raw read by library ####
libnames <- str_extract(string = pathlist,
                        pattern = "[0-9]+-[0-6]")
cell_counts_reads_df <- data.frame(library = libnames,
                                   cell.counts = rep_len(0, length(libnames)),
                                   n.reads = rep_len(0, length(libnames)),
                                   reads.per.cell = rep_len(0, length(libnames)))

for (i in 1:length(libnames)){
  print(libnames[i])
  cell_counts_reads_df$cell.counts[i] <- length(cleanobj_lst[[i]]@assays$ATAC@counts@Dimnames[[2]])
  cell_counts_reads_df$n.reads[i] <- sum(cleanobj_lst[[i]]@assays$ATAC@counts@x)
  cell_counts_reads_df$reads.per.cell[i] <- 
    sum(cleanobj_lst[[i]]@assays$ATAC@counts@x)/cell_counts_reads_df$cell.counts[i]
}

write.table(cell_counts_reads_df, 
            file = "./filtered_cell_counts_UMI_counts_per_lib.csv",
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

# 3. by line: cell counts and raw counts for each cell line ####
lines <- c()
for (i in 1:length(cleanobj_lst)){
  lines <- c(lines, unique(cleanobj_lst[[i]]$cell.line.ident))
}
lines <- sort(unique(lines))
lines
cell_counts_reads_df <- data.frame(line = factor(lines, levels = rev(lines)),
                                   cell.counts = rep_len(0, length(lines)),
                                   n.reads = rep_len(0, length(lines)),
                                   reads.per.cell = rep_len(0, length(lines)))
for (i in 1:length(lines)){
  print(lines[i])
  for (j in cleanobj_lst){
    if (lines[i] %in% unique(j$cell.line.ident)){
      tempobj <- subset(j, cell.line.ident == lines[i])
      cell_counts_reads_df$cell.counts[i] <- cell_counts_reads_df$cell.counts[i] +
        length(tempobj@assays$ATAC@counts@Dimnames[[2]])
      cell_counts_reads_df$n.reads[i] <- cell_counts_reads_df$n.reads[i] +
        sum(tempobj@assays$ATAC@counts@x)
    }
  }
}
cell_counts_reads_df$reads.per.cell <- cell_counts_reads_df$n.reads /
  cell_counts_reads_df$cell.counts
write.table(cell_counts_reads_df, 
            file = "./cell_counts_UMI_counts_per_line.csv",
            quote = F, sep = ",", row.names = F,
            col.names = T)
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

# 4. TSS and nucleosome signal ####

# add annotation
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")
for (i in 1:length(cleanobj_lst)){
  Annotation(cleanobj_lst[[i]]) <- ens_use
}
# compute QC metrics
for (i in 1:length(cleanobj_lst)){
  cleanobj_lst[[i]] <- NucleosomeSignal(cleanobj_lst[[i]])
  cleanobj_lst[[i]] <- TSSEnrichment(object = cleanobj_lst[[i]], fast = FALSE)
  
  jpeg(filename = paste0(str_extract(pathlist[[i]], "[0-9]+-[0|1|6]"),
                         "TSS_plot.jpeg"),
       height = 800, width = 800)
  p <- TSSPlot(cleanobj_lst[[i]]) + 
    NoLegend() +
    theme(text = element_text(size = 10))
  print(p)
  dev.off()
  
  jpeg(filename = paste0(str_extract(pathlist[[i]], "[0-9]+-[0|1|6]"),
                         "nucleosome_signal_hist.jpeg"),
       height = 800, width = 800)
  p <- FragmentHistogram(cleanobj_lst[[i]]) +
    theme(text = element_text(size = 10))
  print(p)
  dev.off()
  
}
TSS_sum <- 0
num_cells <- 0
for (i in cleanobj_lst){
  TSS_sum <- TSS_sum + sum(i$TSS.enrichment)
  num_cells <- num_cells + length(colnames(i))
}
TSS_sum / num_cells

nuc_sig_sum <- 0
num_cells <- 0
for (i in cleanobj_lst){
  print(sum(is.infinite(i$nucleosome_signal)))
  print(sum(is.nan(i$nucleosome_signal)))
  ns <- i$nucleosome_signal[!is.infinite(i$nucleosome_signal)]
  ns <- ns[!is.nan(ns)]
  nuc_sig_sum <- nuc_sig_sum + sum(ns)
  num_cells <- num_cells + length(colnames(i))
}
nuc_sig_sum / num_cells

save(cleanobj_lst, file = "added_TSS_no_rat_no_unmapped_cleanobj_list.RData")


# 5. cellular composition ####
time_type_sum <- vector(mode = "list", length = 3L)
ATAC$cell.type <- factor(ATAC$cell.type, levels = c("NEFM_pos_glut",
                                                                    "NEFM_neg_glut",
                                                                    "GABA", "unknown"))
types <- as.vector(unique(ATAC$cell.type))
numtypes <- length(types)
times <- unique(ATAC$time.ident)
lines <- unique(ATAC$cell.line.ident)
for (k in 1:length(times)){
  print(times[k])
  df <- data.frame(cell.type = rep_len(NA, 18*numtypes),
                   cell.line = rep_len(NA, 18*numtypes),
                   counts = rep_len(0, 18*numtypes))
  for (i in 1:length(lines)){
    print(lines[i])
    subobj <- subset(ATAC, 
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
  jpeg(filename = paste0("cellular_compos_", times[i], ".jpeg"))
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
