# Chuxuan Li 03/21/2022
# get link peak to gene results

# init ####
library(Seurat)
library(Signac)
library(edgeR)

library(readr)
library(stringr)
library(future)

library(ggplot2)
library(RColorBrewer)
library(pheatmap)

plan("multisession", workers = 1)

# read data ####
setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/link_peak_to_genes_results")
pathlist <- list.files("./", pattern = "*.csv")
celltype_time_sep_list <- vector(mode = "list", length = length(pathlist))
for (i in 1:length(pathlist)){
  celltype_time_sep_list[[i]] <- as.data.frame(read_csv(pathlist[i]),row.names = F)
}
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/new_peak_set_after_motif.RData")
obj_complete$celltype.time.ident[is.na(obj_complete$celltype.time.ident)] <- "NA"
obj_complete <- subset(obj_complete, celltype.time.ident != "NA")
idents_list <- unique(obj_complete$celltype.time.ident)
idents_list <- sort(idents_list)
idents_list <- idents_list[4:15]

# how many genes have links ####
linked_genes_list <- vector(mode = "list", length = length(celltype_time_sep_list))
for (i in 1:length(linked_genes_list)){
  linked_genes_list[[i]] <- unique(celltype_time_sep_list[[i]]$gene) 
}
names(linked_genes_list) <- idents_list

linked_gene_counts <- rep_len(0, length(linked_genes_list))
for (i in 1:length(linked_genes_list)){
  linked_gene_counts[i] <- length(linked_genes_list[[i]]) 
}
# use edgeR to calculate cpm for each cell type & time point's genes
# make pseudobulk dataframe
# load data
# Read files, separate the RNAseq data from the .h5 matrix
setwd("~/NVME/scARC_Duan_018/GRCh38_mapped_only")

h5list <- list.files(path = ".", 
                     pattern = "filtered_feature_bc_matrix.h5",
                     recursive = T,
                     include.dirs = T)
objlist <- vector(mode = "list", length = length(h5list))

for (i in 1:6){
  h5file <- Read10X_h5(filename = h5list[i])
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(string = h5list[i], 
                                                  pattern = "libraries_[0-9]_[0-6]"))
  print(paste0(i, " number of genes: ", nrow(obj), 
               ", number of cells: ", ncol(obj)))
  obj$lib.ident <- i
  obj$time.ident <- paste0(str_extract(string = h5list[i], 
                                       pattern = "[0-6]"), "hr")
  objlist[[i]] <- 
    PercentageFeatureSet(obj,
                         pattern = c("^MT-"),
                         col.name = "percent.mt")
  
}

# make pseudobulk dataframe
setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_pseudobulk")
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/27Jan2022_5line_added_cellline_idents.RData")
unique(RNAseq_integrated_labeled$cell.line.ident)
filtered_obj <- subset(RNAseq_integrated_labeled, cell.line.ident != "unmatched")

unique(filtered_obj$broad.cell.type)
GABA <- subset(filtered_obj, subset = cell.type == "GABA")
nmglut <- subset(filtered_obj, subset = cell.type == "NEFM_neg_glut")
npglut <- subset(filtered_obj, subset = cell.type == "NEFM_pos_glut")
NPC <- subset(filtered_obj, cell.type == "NPC")

for (i in 1:length(objlist)){
  print(i)
  obj <- objlist[[i]]
  obj$cell.type.ident <- "unmatched"
  obj$cell.type.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in%
                        str_sub(GABA@assays$RNA@counts@Dimnames[[2]], end = -5L)] <- "GABA"
  obj$cell.type.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in%
                        str_sub(nmglut@assays$RNA@counts@Dimnames[[2]], end = -5L)] <- "NmCp_glut"
  obj$cell.type.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in%
                        str_sub(npglut@assays$RNA@counts@Dimnames[[2]], end = -5L)] <- "NpCm_glut"
  obj$cell.type.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in%
                        str_sub(NPC@assays$RNA@counts@Dimnames[[2]], end = -5L)] <- "NPC"
  print(unique(obj$cell.type.ident))
  objlist[[i]] <- subset(obj, cell.type.ident != "unmatched")
}


genelist <- rownames(objlist[[1]]@assays$RNA@data)
filler <- rep(0, length(genelist))
df_pseudobulk <- data.frame(GABA = filler, nmglut = filler,
                            NPC = filler, npglut = filler, 
                            row.names = genelist)
types <- unique(str_remove_all(idents_list, "_[0|1|6]hr"))
types
for (i in 1:length(objlist)){
  obj <- objlist[[i]]
  j = 0
  for (type in types){
    j = j + 1
    subobj <- subset(obj, cell.type.ident == type)
    print(j)
    df_pseudobulk[j] <- df_pseudobulk[j] + as.array(rowSums(subobj@assays$RNA@counts))
  }
}
df_cpm <- cpm(df_pseudobulk)
hist(df_cpm[, 1], breaks = 100000, xlim = c(0, 50))
sum(df_cpm[, 1] > 1)

setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/link_peak_to_genes_results")
linked_gene_df_list <- vector(mode = "list", length = length(types))
for (i in 1:length(types)){
  df <- data.frame(time = c("0hr", "0hr", "1hr", "1hr", "6hr", "6hr"),
                   count_type = c("Without link", "With link", 
                                  "Without link", "With link",
                                  "Without link", "With link"),
                   count = c((sum((df_cpm[, i] > 1)) - linked_gene_counts[3 * i - 2]),
                             linked_gene_counts[3 * i - 2],
                             (sum((df_cpm[, i] > 1)) - linked_gene_counts[3 * i - 1]),
                             linked_gene_counts[3 * i - 1],
                             (sum((df_cpm[, i] > 1)) - linked_gene_counts[3 * i]),
                             linked_gene_counts[3 * i])
  )
  df$count_type <- factor(df$count_type, levels = c("Without link", "With link"))
  linked_gene_df_list[[i]] <- df
  # pdf(paste0("linked_gene_count_bargraph_", types[i], ".pdf"))
  # p <- ggplot(data = df,
  #             aes(x = time, 
  #                 y = count,
  #                 fill = count_type)) +
  #   geom_col(position = "stack") +
  #   scale_fill_manual(values = (brewer.pal(3, "Set2"))) +
  #   theme_bw() + 
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1),
  #         plot.title = element_text(hjust = 0.5)) + 
  #   xlab("Time") +
  #   ylab("Number of genes with/without links to peaks") +
  #   ggtitle(paste0("Genes linked to peaks in ", types[i]))
  # print(p)
  # dev.off()
}

# how many links per linked gene ####
link_per_gene_counts <- rep_len(0, length(linked_genes_list))
for (i in 1:length(linked_genes_list)){
  link_per_gene_counts[i] <- length(celltype_time_sep_list[[i]]$peak) / length(linked_genes_list[[i]])
}

df <- data.frame(time = c(c("0hr", "1hr", "6hr"), c("0hr", "1hr", "6hr"),
                          c("0hr", "1hr", "6hr"), c("0hr", "1hr", "6hr")),
                 cell_type = str_remove_all(idents_list, "_[0|1|6]hr"),
                 avg = link_per_gene_counts)
cell_type_labs <- list(
  "GABA" = "GABA",
  'NmCp_glut'= "NEFM- glut",
  'NpCm_glut'= "NEFM+ glut",
  'NPC'= "NPC"
)

plot.labeller <- function(variable,value){
    return(cell_type_labs[value])
}

ggplot(data = df,
       aes(x = time, 
           y = avg,
           fill = time)) +
  geom_col() +
  facet_grid(cols = vars(cell_type),
             labeller = plot.labeller) + 
  scale_fill_manual(values = (brewer.pal(3, "Set2"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
         plot.title = element_text(hjust = 0.5)) + 
  xlab("Time") +
  ylab("Average number of links per gene") +
  ggtitle(paste0("Average number of links per gene for each cell type"))

# number of links to a gene - distribution ####
num_links_for_each_gene <- vector(mode = "list", length = length(celltype_time_sep_list))
for (i in 1:length(linked_genes_list)){
  num_links_for_each_gene[[i]] <- table(celltype_time_sep_list[[i]]$gene)
}
hist(num_links_for_each_gene[[12]], xlim = c(0.5, 30))


# check number of up-regulated 1v0/6v0 genes in genes with links ####



# gene activity matrix and gene expression matrix heatmap ####
common_genes <- intersect(celltype_time_sep_list[[1]]$gene,
                          intersect(celltype_time_sep_list[[2]]$gene, 
                                    celltype_time_sep_list[[3]]$gene))
zero_ATAC <- subset(obj_complete, celltype.time.ident == "GABA_0hr")
one_ATAC <- subset(obj_complete, celltype.time.ident == "GABA_1hr")
six_ATAC <- subset(obj_complete, celltype.time.ident == "GABA_6hr")

# comb 1: RNA scale.data x peaks
DefaultAssay(filtered_obj) <- "RNA"
filtered_obj <- ScaleData(filtered_obj)
common_genes <- intersect(common_genes, filtered_obj@assays$RNA@data@Dimnames[[1]])
unique(obj_complete$celltype.time.ident)

zero_gact <- GeneActivity(zero_ATAC, assay = "peaks", features = common_genes, 
             extend.downstream = 0, extend.upstream = 1000)
one_gact <- GeneActivity(one_ATAC, assay = "peaks", features = common_genes, 
                          extend.downstream = 0, extend.upstream = 1000)
six_gact <- GeneActivity(six_ATAC, assay = "peaks", features = common_genes, 
                          extend.downstream = 0, extend.upstream = 1000)

# gene expression matrix
zero_gexp <- subset(filtered_obj, time.ident == "0hr")@assays$RNA@scale.data
zero_gexp <- zero_gexp[filtered_obj@assays$RNA@counts@Dimnames[[1]] %in% zero_gact@Dimnames[[1]], ]
one_gexp <- subset(filtered_obj, time.ident == "1hr")@assays$RNA@scale.data
one_gexp <- one_gexp[filtered_obj@assays$RNA@counts@Dimnames[[1]] %in% one_gact@Dimnames[[1]], ]
six_gexp <- subset(filtered_obj, time.ident == "6hr")@assays$RNA@scale.data
six_gexp <- six_gexp[filtered_obj@assays$RNA@counts@Dimnames[[1]] %in% six_gact@Dimnames[[1]], ]


# comb 2: RNA scale.data x ATAC
DefaultAssay(filtered_obj) <- "RNA"
filtered_obj <- ScaleData(filtered_obj)
common_genes <- intersect(common_genes, filtered_obj@assays$RNA@data@Dimnames[[1]])
unique(obj_complete$celltype.time.ident)
zero_gact <- GeneActivity(zero_ATAC, assay = "ATAC", features = common_genes, 
                          extend.downstream = 0, extend.upstream = 1000)
one_gact <- GeneActivity(one_ATAC, assay = "ATAC", features = common_genes, 
                         extend.downstream = 0, extend.upstream = 1000)
six_gact <- GeneActivity(six_ATAC, assay = "ATAC", features = common_genes, 
                         extend.downstream = 0, extend.upstream = 1000)

# gene expression matrix
zero_gexp <- subset(filtered_obj, time.ident == "0hr")@assays$RNA@scale.data
zero_gexp <- zero_gexp[filtered_obj@assays$RNA@counts@Dimnames[[1]] %in% zero_gact@Dimnames[[1]], ]
one_gexp <- subset(filtered_obj, time.ident == "1hr")@assays$RNA@scale.data
one_gexp <- one_gexp[filtered_obj@assays$RNA@counts@Dimnames[[1]] %in% one_gact@Dimnames[[1]], ]
six_gexp <- subset(filtered_obj, time.ident == "6hr")@assays$RNA@scale.data
six_gexp <- six_gexp[filtered_obj@assays$RNA@counts@Dimnames[[1]] %in% six_gact@Dimnames[[1]], ]


# comb 3: SCT x peaks
DefaultAssay(filtered_obj) <- "SCT"
common_genes <- intersect(common_genes, filtered_obj@assays$SCT@data@Dimnames[[1]])
unique(obj_complete$celltype.time.ident)
zero_gact <- GeneActivity(zero_ATAC, assay = "peaks", features = common_genes, 
                          extend.downstream = 0, extend.upstream = 1000)
one_gact <- GeneActivity(one_ATAC, assay = "peaks", features = common_genes, 
                         extend.downstream = 0, extend.upstream = 1000)
six_gact <- GeneActivity(six_ATAC, assay = "peaks", features = common_genes, 
                         extend.downstream = 0, extend.upstream = 1000)

# gene expression matrix
zero_gexp <- subset(filtered_obj, time.ident == "0hr")@assays$SCT
zero_gexp <- zero_gexp[filtered_obj@assays$RNA@counts@Dimnames[[1]] %in% zero_gact@Dimnames[[1]], ]
one_gexp <- subset(filtered_obj, time.ident == "1hr")@assays$SCT
one_gexp <- one_gexp[filtered_obj@assays$RNA@counts@Dimnames[[1]] %in% one_gact@Dimnames[[1]], ]
six_gexp <- subset(filtered_obj, time.ident == "6hr")@assays$SCT
six_gexp <- six_gexp[filtered_obj@assays$RNA@counts@Dimnames[[1]] %in% six_gact@Dimnames[[1]], ]

# make matrices to plot
gene_activity_matrix <- data.frame(rbind(rowSums(zero_gact),
                                         rowSums(one_gact),
                                         rowSums(six_gact)),
                                   stringsAsFactors = F)
rownames(gene_activity_matrix) <- c("0hr", "1hr", "6hr")
gene_activity_matrix <- as.matrix(t(gene_activity_matrix))
colnames(gene_activity_matrix)


gene_expression_matrix <- data.frame(rbind(rowSums(zero_gexp),
                                           rowSums(one_gexp),
                                           rowSums(six_gexp)),
                                     stringsAsFactors = F)
rownames(gene_expression_matrix) <- c("0hr", "1hr", "6hr")
gene_expression_matrix <- as.matrix(t(gene_expression_matrix))
colnames(gene_expression_matrix)

# plot
annot_df <- rownames(gene_expression_matrix)
annot_df[!(annot_df %in% c("FOS", "NPAS4", "EGR1", "FOSB", "NR4A1", "IGF1", "BDNF"))] <- ""
test_gact <- gene_activity_matrix[p_gexp$tree_row$order, ]
p_gact <- pheatmap(gene_activity_matrix, scale = "row",
                   cluster_rows = F, cluster_cols = F, 
                   treeheight_row = 0, labels_row = annot_df,
                   border_color = NA, fontsize = 8, main = "Gene Activity")


test_gexp <- gene_expression_matrix[rowSums(gene_expression_matrix) > 0, ]
sum(is.na(rowSums(gene_expression_matrix)))
sum(is.infinite((rowSums(gene_expression_matrix))))
test_gexp <- gene_expression_matrix[apply(gene_expression_matrix, 1, var) != 0, ]
test_gexp <- test_gexp[p_gact$tree_row$order, ]

p_gexp <- pheatmap(test_gexp, scale = "row",
              cluster_rows = T, cluster_cols = F, 
              treeheight_row = 0, labels_row = annot_df,
              border_color = NA, fontsize = 8, main = "Gene Expression")

