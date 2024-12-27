# Chuxuan Li 03/04/2022
# Apply the new RNAseq cell types on ATACseq data from ATAcseq analysis v2 (called
#new peaks using Signac embedded MACS2) and downstream analysis

# init ####
library(Seurat)
library(stringr)
library(future)

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(304)
# set threads and parallelization
plan("multisession", workers = 1)

# load data
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/ATACseq_5line_objonly_labeled.RData")
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/human_only_obj_with_cell_type_labeling.RData")

# project RNAseq cell type ident onto the ATACseq data ####
obj_complete$RNA.cell.type.ident <- "Unmatched"
# find all barcodes corresponding to each cell type
obj_complete$trimmed.barcodes <- str_sub(obj_complete@assays$peaks@data@Dimnames[[2]], end = -3L)
for (i in unique(human_only$cell.type)){
  barcodes <- str_sub(human_only@assays$RNA@counts@Dimnames[[2]][human_only$cell.type %in% i], 
                      end = -5L)
  
  obj_complete$RNA.cell.type.ident[obj_complete$trimmed.barcodes %in% barcodes] <- i
}

obj_complete$RNA.cluster.ident <- "Unmatched"
for (i in unique(human_only$seurat_clusters)){
  barcodes <- str_sub(human_only@assays$RNA@counts@Dimnames[[2]][human_only$seurat_clusters %in% i], 
                      end = -5L)
  
  obj_complete$RNA.cluster.ident[obj_complete$trimmed.barcodes %in% barcodes] <- i
}

unique(obj_complete$RNA.cluster.ident)

DimPlot(obj_complete,
        group.by = "seurat_clusters",
        label = F) +
  ggtitle("ATACseq clustering results") +
  theme(text = element_text(size = 10))
DimPlot(obj_complete,
        group.by = "RNA.cluster.ident",
        label = T,
        repel = T) +
  ggtitle("ATACseq data projected by RNAseq cluster numbers") +
  theme(text = element_text(size = 10))
DimPlot(obj_complete, group.by = "RNA.cell.type.ident", 
        label = T, 
        repel = T) +
  ggtitle("ATACseq clusters labeled by RNAseq cell types") +
  theme(text = element_text(size = 10))

# project ATACseq clusters onto RNAseq clusters
human_only$ATAC.cluster.ident <- "Unmatched"
human_only$trimmed.barcodes <- str_sub(human_only@assays$RNA@counts@Dimnames[[2]], 
        end = -5L)
for (i in unique(obj_complete$seurat_clusters)){
  barcodes <- obj_complete$trimmed.barcodes[obj_complete$seurat_clusters%in% i]
  human_only$ATAC.cluster.ident[human_only$trimmed.barcodes %in% barcodes] <- i
}
unique(human_only$ATAC.cluster.ident)
DimPlot(human_only,
        group.by = "ATAC.cluster.ident",
        label = T,
        repel = T) +
  ggtitle("RNAseq data projected by ATACseq cluster numbers") +
  theme(text = element_text(size = 10))
# check cell counts for each cell type
obj_complete$RNA.cell.type.ident[obj_complete$RNA.cell.type.ident == "Unmatched"] <- "unknown"
types <- unique(obj_complete$RNA.cell.type.ident)
cell_counts_reads_df <- data.frame(cell.type = rep_len(NA, length(types)),
                                   cell.counts = rep_len(0, length(types)),
                                   n.reads = rep_len(0, length(types)),
                                   reads.per.cell = rep_len(0, length(types)))

for (i in 1:length(types)){
  print(types[i])
  cell_counts_reads_df$cell.line <- types[i]
  cell_counts_reads_df$cell.counts[i] <- sum(obj_complete$RNA.cell.type.ident == types[i])
  temp <- subset(obj_complete, RNA.cell.type.ident == types[i])
  cell_counts_reads_df$n.reads[i] <- sum(temp@assays$peaks@counts@x)
  cell_counts_reads_df$reads.per.cell[i] <- 
    sum(temp@assays$peaks@counts@x)/cell_counts_reads_df$cell.counts[i]
}

# plot the cell count first
ggplot(cell_counts_reads_df,
            aes(x = types,
                y = rep_len(2, length(types)),
                fill = cell.counts)) +
  geom_col(width = 1, 
           position = position_dodge(width = 0.5)) +
  geom_text(aes(label = cell.counts), 
            nudge_y = -1,
            hjust = "middle") + 
  scale_fill_continuous(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0, size = 11, vjust = -1),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  NoLegend()

# plot UMI per cell next
ggplot(cell_counts_reads_df,
            aes(x = types,
                y = reads.per.cell)) +
  ylab("# reads per cell") +
  geom_col(position = "dodge",
           color = "gray27",
           fill = "gray90") +
  #geom_text(aes(label = reads.per.cell)) + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 11)) 


# plot UMI total for each cell line
ggplot(cell_counts_reads_df,
            aes(x = types,
                y = n.reads)) +
  ylab("# total reads") +
  geom_col(position = "dodge",
           color = "gray27",
           fill = "gray90") +
  #geom_text(aes(label = reads.per.cell)) + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 11),
        axis.text.y = element_text(size = 8)) 

# remove other cell types
obj_part <- subset(obj_complete, RNA.cell.type.ident %in% c("GABA", 
                                                            "NEFM_pos_glut",
                                                            "NEFM_neg_glut",
                                                            "NPC"))
DimPlot(obj_part, group.by = "RNA.cell.type.ident", 
        label = T, 
        repel = T) +
  ggtitle("ATACseq clusters labeled by RNAseq cell types\nwith four main cell types") +
  theme(text = element_text(size = 10)) +
  NoLegend()


# 
