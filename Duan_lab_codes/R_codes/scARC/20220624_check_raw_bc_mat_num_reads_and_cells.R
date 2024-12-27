# Chuxuan Li 06/24/2022
# check how many cells and reads are in each library in raw_bc_matrix

# init ####
library(Seurat)
library(glmGamPoi)
library(sctransform)

library(ggplot2)
library(RColorBrewer)
library(viridis)
library(graphics)
library(ggrepel)

library(dplyr)
library(readr)
library(stringr)
library(readxl)

library(future)
# set threads and parallelization
plan("multisession", workers = 2)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

# load data ####
setwd("/data/FASTQ/Duan_Project_024/hybrid_output/")
# note this h5 file contains both atac-seq and gex information
# Read files, separate the ATACseq data from the .h5 matrix
h5list <- list.files(path = ".",
                     pattern = "raw_feature_bc_matrix.h5",
                     recursive = T,
                     include.dirs = F, 
                     full.names = T)
h5list <- h5list[1:15]
objlist <- vector(mode = "list", length = length(h5list))

for (i in 1:length(objlist)){
  
  h5file <- Read10X_h5(filename = h5list[i])
  print(str_extract(string = h5list[i],
                    pattern = "[0-9]+-[0-6]"))
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(string = h5list[i],
                                                  pattern = "[0-9]+-[0-6]"))
  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list[i],
                                                  pattern = "-[0-6]"), "-"), "hr")
  print(unique(obj$time.ident))
  # check number of genes and cells
  print(paste0(i, " number of reads: ", sum(colSums(obj)),
               ", number of cells: ", ncol(obj)))
  objlist[[i]] <- obj
}

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



# ATAC ####
# load data ####
setwd("/data/FASTQ/Duan_Project_024/hybrid_output/")
# note this h5 file contains both atac-seq and gex information
# Read files, separate the ATACseq data from the .h5 matrix
h5list <- list.files(path = ".",
                     pattern = "raw_feature_bc_matrix.h5",
                     recursive = T,
                     include.dirs = F, 
                     full.names = T)
h5list <- h5list[1:15]
objlist <- vector(mode = "list", length = length(h5list))

for (i in 1:length(objlist)){
  
  h5file <- Read10X_h5(filename = h5list[i])
  print(str_extract(string = h5list[i],
                    pattern = "[0-9]+-[0-6]"))
  obj <- CreateSeuratObject(counts = h5file$Peaks,
                            project = str_extract(string = h5list[i],
                                                  pattern = "[0-9]+-[0-6]"))
  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list[i],
                                                  pattern = "-[0-6]"), "-"), "hr")
  print(unique(obj$time.ident))
  # check number of genes and cells
  print(paste0(i, " number of reads: ", sum(colSums(obj)),
               ", number of cells: ", ncol(obj)))
  objlist[[i]] <- obj
}
