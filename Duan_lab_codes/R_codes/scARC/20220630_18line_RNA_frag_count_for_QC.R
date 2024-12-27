# Chuxuan Li 06/30/2022
# using fragment file instead of filtered bc matrix and raw bc matrix h5 files
#to count number of reads in each library - this is for 18line RNAseq
# based on /nvmefs/scARC_10x_PBMC_10K/calc_read_counts.R

# init ####
library(Seurat)
library(Signac)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(stringr)
library(readr)

setwd(dir = "/data/FASTQ/Duan_Project_024/hybrid_output")

# read files ####
frag_files <-
  list.files(path = ".",
             pattern = "atac_fragments.tsv.gz$",
             full.names = T,
             recursive = T,
             include.dirs = F)
frag_files <-
  frag_files[str_detect(string = frag_files,
                        pattern = "group_")]


raw_bcmatrix_files <-
  list.files(path = ".",
             pattern = "raw_feature_bc_matrix.h5$",
             full.names = T,
             recursive = T,
             include.dirs = F)
raw_bcmatrix_files <-
  raw_bcmatrix_files[str_detect(string = raw_bcmatrix_files,
                                pattern = "group_")]

filtered_bcmatrix_files <-
  list.files(path = ".",
             pattern = "filtered_feature_bc_matrix.h5$",
             full.names = T,
             recursive = T,
             include.dirs = F)
filtered_bcmatrix_files <-
  filtered_bcmatrix_files[str_detect(string = filtered_bcmatrix_files,
                                     pattern = "group_")]
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/demux_mapped_list.RData")

# make a table ####
# using cell barcode index from the corresponding indices
sample_index <- str_extract(string = frag_files,
                            pattern = "[0-9]+-[0-6]")

df_raw_filtered_table <-
  data.frame(library = sample_index,
             raw.cell.counts = rep_len(0, length(sample_index)),
             raw.fragment.count = rep_len(0, length(sample_index)),
             raw.fragment.per.cell = rep_len(0, length(sample_index)),
             human.cell.counts = rep_len(0, length(sample_index)),
             human.fragment.count = rep_len(0, length(sample_index)),
             human.fragment.per.cell = rep_len(0, length(sample_index)))

for (i in 1:length(frag_files)) {
  print(df_raw_filtered_table$library[i])
  barcodes_h5_filtered <-
    Read10X_h5(filename = filtered_bcmatrix_files[i],
               use.names = T)$`Gene Expression`
  barcodes_h5_raw <-
    Read10X_h5(filename = raw_bcmatrix_files[i],
               use.names = T)$`Gene Expression`
  
  # raw
  df_raw_filtered_table$raw.cell.counts[i] <- ncol(barcodes_h5_raw)
  df_raw_filtered_table$raw.fragment.count[i] <-
    sum(CountFragments(fragments = frag_files[i],
                       cells = colnames(barcodes_h5_raw),
                       verbose = T)$frequency_count)
  df_raw_filtered_table$raw.fragment.per.cell[i] <- 
    df_raw_filtered_table$raw.fragment.count[i]/df_raw_filtered_table$raw.cell.counts[i]
  
  # filtered
  df_raw_filtered_table$filtered.cell.counts[i] <- ncol(barcodes_h5_filtered)
  df_raw_filtered_table$filtered.fragment.count[i] <-
    sum(CountFragments(fragments = frag_files[i],
                       cells = colnames(barcodes_h5_filtered),
                       verbose = T)$frequency_count)
  df_raw_filtered_table$filtered.fragment.per.cell[i] <- 
    df_raw_filtered_table$filtered.fragment.count[i]/df_raw_filtered_table$filtered.cell.counts[i]
}
# human only cells
for (i in 1:length(sample_index)){
  print(sample_index[i])
  print(unique(cleanobj_lst[[i]]$orig.ident))
  df_raw_filtered_table$human.cell.counts[i] <- ncol(cleanobj_lst[[i]])
  df_raw_filtered_table$human.fragment.count[i] <- sum(CountFragments(fragments = frag_files[i],
                                                                      cells = colnames(cleanobj_lst[[i]]),
                                                                      verbose = T)$frequency_count)
  df_raw_filtered_table$human.fragment.per.cell[i] <- 
    df_raw_filtered_table$human.fragment.count[i]/df_raw_filtered_table$human.cell.counts[i]
}

setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/")
write.table(df_raw_filtered_table, 
            file = "./Analysis_v1_normalize_by_lib_then_integrate/QC_dataframes/RNA_read_counts_from_fragment_file_raw_filtered.csv", 
            quote = F, sep = ",", row.names = F, col.names = T)

# plot the cell count ####
count_table <- read_csv("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/QC_dataframes/RNA_read_counts_from_fragment_file_raw_filtered.csv")

# raw
ggplot(count_table,
       aes(x = factor(library, levels = rev(c("5-0", "5-1", "5-6", "33-0", "33-1", "33-6",
                                              "35-0", "35-1", "35-6", "51-0", "51-1", "51-6",
                                              "63-0", "63-1", "63-6"))),
           y = rep_len(2, length(sample_index)),
           fill = raw.cell.counts)) +
  geom_col(width = 1, 
           position = position_dodge(width = 0.5)) +
  geom_text(aes(label = raw.cell.counts), 
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
# filtered
ggplot(count_table,
       aes(x = factor(library, levels = rev(c("5-0", "5-1", "5-6", "33-0", "33-1", "33-6",
                                              "35-0", "35-1", "35-6", "51-0", "51-1", "51-6",
                                              "63-0", "63-1", "63-6"))),
           y = rep_len(2, length(sample_index)),
           fill = filtered.cell.counts)) +
  geom_col(width = 1, 
           position = position_dodge(width = 0.5)) +
  geom_text(aes(label = filtered.cell.counts), 
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
# human
ggplot(count_table,
       aes(x = factor(library, levels = rev(c("5-0", "5-1", "5-6", "33-0", "33-1", "33-6",
                                              "35-0", "35-1", "35-6", "51-0", "51-1", "51-6",
                                              "63-0", "63-1", "63-6"))),
           y = rep_len(2, length(sample_index)),
           fill = human.cell.counts)) +
  geom_col(width = 1, 
           position = position_dodge(width = 0.5)) +
  geom_text(aes(label = human.cell.counts), 
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

# plot UMI per cell ####
# raw
ggplot(count_table,
       aes(x = library,
           y = raw.fragment.per.cell)) +
  ylab("# reads per nuclei") +
  geom_col(position = "dodge",
           color = "gray27",
           fill = "gray90") +
  #geom_text(aes(label = reads.per.cell)) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 11)) + 
  coord_flip()
# filtered
ggplot(count_table,
       aes(x = library,
           y = filtered.fragment.per.cell)) +
  ylab("# reads per nuclei") +
  geom_col(position = "dodge",
           color = "gray27",
           fill = "gray90") +
  #geom_text(aes(label = reads.per.cell)) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 11)) + 
  coord_flip()
# human only
ggplot(count_table,
       aes(x = library,
           y = human.fragment.per.cell)) +
  ylab("# reads per nuclei") +
  geom_col(position = "dodge",
           color = "gray27",
           fill = "gray90") +
  #geom_text(aes(label = reads.per.cell)) + 
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 11)) + 
  coord_flip()

# plot UMI total ####
ggplot(count_table,
       aes(x = library,
           y = raw.fragment.count)) +
  ylab("# total reads") +
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
# filtered
ggplot(count_table,
       aes(x = library,
           y = filtered.fragment.count)) +
  ylab("# total reads") +
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
# human only
ggplot(count_table,
       aes(x = library,
           y = human.fragment.count)) +
  ylab("# total reads") +
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
