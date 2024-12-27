# Chuxuan Li 03/03/2022
# create gene activity matrix using 5000 variable genes from RNAseq analysis v3,
#ATACseq object from ATAcseq analysis v2 (new_peak_set_after_motif.RData),
#then plot gene activity and gene expression for 5000 genes by time

# init ####
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(future)
set.seed(1234)

plan("multisession", workers = 1)

load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/mapped_to_demuxed_barcodes_pure_human_obj.RData")
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/ATACseq_5line_objonly_labeled.RData")

unique(pure_human$lib.ident)

# preprocess data ####
# get Seurat object with only peaks assay
DefaultAssay(obj_complete) <- "ATAC"
zero_ATAC <- subset(obj_complete, time.ident == "0hr")
one_ATAC <- subset(obj_complete, time.ident == "1hr")
six_ATAC <- subset(obj_complete, time.ident == "6hr")
zero_ATAC[["RNA"]] <- NULL
zero_ATAC[["SCT"]] <- NULL
zero_ATAC[["peaks"]] <- NULL
one_ATAC[["RNA"]] <- NULL
one_ATAC[["SCT"]] <- NULL
one_ATAC[["peaks"]] <- NULL
six_ATAC[["RNA"]] <- NULL
six_ATAC[["SCT"]] <- NULL
six_ATAC[["peaks"]] <- NULL
zero_gact <- GeneActivity(zero_ATAC, 
                          features = pure_human@assays$integrated@data@Dimnames[[1]],
                          extend.upstream = 1000,
                          extend.downstream = 0)
one_gact <- GeneActivity(one_ATAC, 
                          features = pure_human@assays$integrated@data@Dimnames[[1]],
                          extend.upstream = 1000,
                          extend.downstream = 0)
six_gact <- GeneActivity(six_ATAC, 
                          features = pure_human@assays$integrated@data@Dimnames[[1]],
                          extend.upstream = 1000,
                          extend.downstream = 0)

# gene expression matrix
zero_gexp <- subset(pure_human, time.ident == "0hr")@assays$integrated@data
zero_gexp <- zero_gexp[zero_gexp@Dimnames[[1]] %in% zero_gact@Dimnames[[1]], ]
one_gexp <- subset(pure_human, time.ident == "1hr")@assays$integrated@data
one_gexp <- one_gexp[one_gexp@Dimnames[[1]] %in% one_gact@Dimnames[[1]], ]
six_gexp <- subset(pure_human, time.ident == "6hr")@assays$integrated@data
six_gexp <- six_gexp[six_gexp@Dimnames[[1]] %in% six_gact@Dimnames[[1]], ]

DefaultAssay(pure_human) <- "RNA"
ScaleData(pure_human)

zero_gexp <- subset(pure_human, time.ident == "0hr")@assays$RNA@scale.data
zero_gexp <- zero_gexp[zero_gexp@Dimnames[[1]] %in% zero_gact@Dimnames[[1]], ]
one_gexp <- subset(pure_human, time.ident == "1hr")@assays$RNA@scale.data
one_gexp <- one_gexp[one_gexp@Dimnames[[1]] %in% one_gact@Dimnames[[1]], ]
six_gexp <- subset(pure_human, time.ident == "6hr")@assays$RNA@scale.data
six_gexp <- six_gexp[six_gexp@Dimnames[[1]] %in% six_gact@Dimnames[[1]], ]

save(zero_gexp, one_gexp, six_gexp, zero_gact, one_gact, six_gact, 
     file = "gact_gexp_matrices.RData")

# plot heatmap ####
library(gplots)
library(colorspace)
library(pheatmap)

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

# plot gene activity
test <- gene_activity_matrix[rowSums(gene_activity_matrix) > 0, ]
sum(is.na(rowSums(gene_activity_matrix)))
sum(is.infinite((rowSums(gene_activity_matrix))))
test <- gene_activity_matrix[apply(gene_activity_matrix, 1, var) != 0, ]

annot_df <- rownames(gene_activity_matrix)
annot_df[!(annot_df %in% c("FOS", "NPAS4", "EGR1", "FOSB", "NR4A1", "IGF1", "BDNF"))] <- ""
p <- pheatmap(test, scale = "row", 
         cluster_rows = T, cluster_cols = F, 
         clustering_distance_rows = "correlation", 
         treeheight_row = 0, labels_row = annot_df,
         border_color = NA, fontsize = 8, main = "Gene Activity")

add.flag(p,
         kept.labels = annot_df,
         repel.degree = 0.1)

# plot gene expression
test <- gene_expression_matrix[rowSums(gene_expression_matrix) > 0, ]
sum(is.na(rowSums(gene_expression_matrix)))
sum(is.infinite((rowSums(gene_expression_matrix))))
test <- gene_expression_matrix[apply(gene_expression_matrix, 1, var) != 0, ]

annot_df <- rownames(gene_expression_matrix)
annot_df[!(annot_df %in% c("FOS", "NPAS4", "EGR1", "FOSB", "NR4A1", "IGF1", "BDNF"))] <- ""
p <- pheatmap(test, scale = "row",
              cluster_rows = T, cluster_cols = F, 
              clustering_distance_rows = "correlation", 
              treeheight_row = 0, labels_row = annot_df,
              border_color = NA, fontsize = 8, main = "Gene Expression")

add.flag(p,
         kept.labels = annot_df,
         repel.degree = 0.05)
hist(test)
gene_activity_matrix[rownames(gene_activity_matrix) == "NPAS4", ]
