# 27 Jan 2022 Siwei
# run pseudo-bulk RNA-seq to check some key gene expression

# init
library(edgeR)


load("/nvmefs/scARC_Duan_018/R_pseudo_bulk_RNA_seq/data_frames_for_pseudobulk_new.RData")

##### GABA ####

GABA_DGE <- DGEList(counts = as.matrix(GABA_pseudobulk), 
                    samples = colnames(GABA_pseudobulk),
                    genes = rownames(GABA_pseudobulk),
                    remove.zeros = T)

# filter out very low expression genes
# genes with cpm < 1 in at least 2/3 of the samples
GABA_DGE <- GABA_DGE[rowSums(cpm(GABA_DGE) > 1) >= 10, , 
                     keep.lib.sizes = F]

GABA_targets <- data.frame(cell_line = c(rep_len("CD_08", length.out = 3),
                                         rep_len("CD_25", length.out = 3),
                                         rep_len("CD_26", length.out = 3),
                                         rep_len("CD_27", length.out = 3),
                                         rep_len("CD_54", length.out = 3)),
                           time_point = c(rep_len(c("0hr", "1hr", "6hr"),
                                                  length.out = 5)))
rownames(GABA_targets) <- colnames(GABA_pseudobulk)

# normalise data series
GABA_DGE <- calcNormFactors(GABA_DGE)
GABA_DGE <- estimateDisp(GABA_DGE)

hist(cpm(GABA_DGE), breaks = 1000, xlim = c(0, 500))

# make design matrix, add blocking factor (cell_line)
GABA_design_matrix <- model.matrix(~ 0 + 
                                     GABA_targets$time_point +
                                     GABA_targets$cell_line)

# Run GLMQLFit with design matrix
GABA_fit <- glmQLFit(GABA_DGE, 
                      design = GABA_design_matrix)
GABA_1vs0 <- glmQLFTest(GABA_fit,
                        contrast = c(-1, 1, 0, 0, 0, 0, 0))

GABA_results_1vs0 <- GABA_1vs0$table

GABA_6vs0 <- glmQLFTest(GABA_fit,
                        contrast = c(-1, 0, 1, 0, 0, 0, 0))

GABA_results_6vs0 <- GABA_6vs0$table

GABA_results_1vs0[rownames(GABA_results_1vs0) %in% "BDNF", ]
GABA_results_6vs0[rownames(GABA_results_6vs0) %in% "BDNF", ]
GABA_results_1vs0[rownames(GABA_results_1vs0) %in% "GAPDH", ]
GABA_results_6vs0[rownames(GABA_results_6vs0) %in% "GAPDH", ]

library(ggplot2)

GABA_results_6vs0$FDR <- p.adjust(GABA_results_6vs0$PValue,
                                  method = "fdr")
GABA_results_6vs0$sign <- "black"
GABA_results_6vs0$sign[GABA_results_6vs0$FDR < 0.05] <- "red"

sum(GABA_results_6vs0$FDR < 0.05)

sum(GABA_results_6vs0$logFC > 0)
sum(GABA_results_6vs0$logFC < 0)


ggplot(GABA_results_6vs0,
       aes(x = logFC,
           y = (0 - log(PValue)))) +
  geom_point(size = 0.2) +
  theme_classic()

ggplot(Glut_results_6vs0,
       aes(x = logFC,
           y = (0 - log(PValue)))) +
  geom_point(color = Glut_results_6vs0$sign) +
  theme_classic()


##### Glut ####

Glut_DGE <- DGEList(counts = as.matrix(glut_pseudobulk), 
                    samples = colnames(glut_pseudobulk),
                    genes = rownames(glut_pseudobulk),
                    remove.zeros = T)

Glut_DGE <- Glut_DGE[rowSums(cpm(Glut_DGE) > 1) >= 10, , 
                     keep.lib.sizes = F]

Glut_targets <- data.frame(cell_line = c(rep_len("CD_08", length.out = 3),
                                         rep_len("CD_25", length.out = 3),
                                         rep_len("CD_26", length.out = 3),
                                         rep_len("CD_27", length.out = 3),
                                         rep_len("CD_54", length.out = 3)),
                           time_point = c(rep_len(c("0hr", "1hr", "6hr"),
                                                  length.out = 5)))
rownames(Glut_targets) <- colnames(glut_pseudobulk)

Glut_DGE <- calcNormFactors(Glut_DGE)
Glut_DGE <- estimateDisp(Glut_DGE)

hist(cpm(Glut_DGE), breaks = 1000, xlim = c(0, 500))

Glut_design_matrix <- model.matrix(~ 0 +
                                     Glut_targets$time_point +
                                     Glut_targets$cell_line)
Glut_design_matrix <- model.matrix(~ 0 +
                                     Glut_targets$time_point )

Glut_fit <- glmQLFit(Glut_DGE, 
                     design = Glut_design_matrix)
Glut_1vs0 <- glmQLFTest(Glut_fit,
                        contrast = c(-1, 1, 0, 0, 0, 0, 0))

Glut_results_1vs0 <- Glut_1vs0$table

Glut_6vs0 <- glmQLFTest(Glut_fit,
                        contrast = c(-1, 0, 1, 0, 0, 0, 0))

Glut_results_6vs0 <- Glut_6vs0$table

Glut_results_1vs0[rownames(Glut_results_1vs0) %in% "VGF", ]
Glut_results_6vs0[rownames(Glut_results_6vs0) %in% "VGF", ]

Glut_results_1vs0[rownames(Glut_results_1vs0) %in% "BDNF", ]
Glut_results_6vs0[rownames(Glut_results_6vs0) %in% "BDNF", ]

Glut_results_1vs0[rownames(Glut_results_1vs0) %in% "GAPDH", ]
Glut_results_6vs0[rownames(Glut_results_6vs0) %in% "GAPDH", ]


Glut_results_1vs0$FDR <- p.adjust(Glut_results_1vs0$PValue, 
                                  method = "fdr")
Glut_results_6vs0$FDR <- p.adjust(Glut_results_6vs0$PValue, 
                                  method = "fdr")
GABA_results_1vs0$FDR <- p.adjust(GABA_results_1vs0$PValue, 
                                  method = "fdr")
GABA_results_6vs0$FDR <- p.adjust(GABA_results_6vs0$PValue, 
                                  method = "fdr")

write.table(GABA_results_1vs0, 
            file = "GABA_results_pseudobulk_1vs0.tsv", 
            row.names = T, col.names = T, 
            sep = "\t", quote = F)
write.table(GABA_results_6vs0, 
            file = "GABA_results_pseudobulk_6vs0.tsv", 
            row.names = T, col.names = T, 
            sep = "\t", quote = F)
write.table(Glut_results_1vs0, 
            file = "Glut_results_pseudobulk_1vs0.tsv", 
            row.names = T, col.names = T, 
            sep = "\t", quote = F)
write.table(Glut_results_6vs0, 
            file = "Glut_results_pseudobulk_6vs0.tsv", 
            row.names = T, col.names = T, 
            sep = "\t", quote = F)


maPlot(0-log10(Glut_results_1vs0$PValue),
       y = Glut_results_1vs0$logFC,
       na.rm = T)

###### use DeSeq2

library(DESeq2)
library(BiocParallel)

##
load("/nvmefs/scARC_Duan_018/R_pseudo_bulk_RNA_seq/data_frames_for_pseudobulk_new.RData")

Glut_raw <- as.matrix(glut_pseudobulk)

Glut_targets <- data.frame(cell_line = c(rep_len("CD_08", length.out = 3),
                                         rep_len("CD_25", length.out = 3),
                                         rep_len("CD_26", length.out = 3),
                                         rep_len("CD_27", length.out = 3),
                                         rep_len("CD_54", length.out = 3)),
                           time_point = c(rep_len(c("0hr", "1hr", "6hr"),
                                                  length.out = 5)))
rownames(Glut_targets) <- colnames(glut_pseudobulk)

Glut_dds <- DESeqDataSetFromMatrix(countData = Glut_raw,
                                   colData = Glut_targets,
                                   design = ~ cell_line + time_point)
# filter all 0 rows
Glut_dds <- Glut_dds[rowSums(counts(Glut_dds)) > 0, ]

Glut_dds$time_point <- relevel(Glut_dds$time_point, ref = "0hr")

Glut_dds <- DESeq(Glut_dds,
                  test = "Wald")

Glut_results_6vs0 <- results(Glut_dds,
                             contrast = c("time_point", "6hr", "0hr"),
                             pAdjustMethod = "fdr")

Glut_results_6vs0 <- as.data.frame(Glut_results_6vs0)

Glut_results_1vs0 <- results(Glut_dds,
                             contrast = c("time_point", "1hr", "0hr"),
                             pAdjustMethod = "fdr")

Glut_results_1vs0 <- as.data.frame(Glut_results_1vs0)

Glut_results_1vs0[rownames(Glut_results_1vs0) %in% "BDNF", ]
Glut_results_6vs0[rownames(Glut_results_6vs0) %in% "BDNF", ]

Glut_results_1vs0[rownames(Glut_results_1vs0) %in% "VGF", ]
Glut_results_6vs0[rownames(Glut_results_6vs0) %in% "VGF", ]

Glut_results_1vs0[rownames(Glut_results_1vs0) %in% "FOS", ]
Glut_results_6vs0[rownames(Glut_results_6vs0) %in% "FOS", ]

Glut_results_1vs0[rownames(Glut_results_1vs0) %in% "GAPDH", ]
Glut_results_6vs0[rownames(Glut_results_6vs0) %in% "GAPDH", ]

library(ggplot2)

ggplot(Glut_results_1vs0,
       aes(x = log2FoldChange,
           y = (0 - log10(pvalue)))) +
  geom_point(size = 0.2) +
  theme_classic()

ggplot(Glut_results_6vs0,
       aes(x = log2FoldChange,
           y = (0 - log10(pvalue)))) +
  geom_point(size = 0.2) +
  xlim(-5, 5) +
  ylim(0, 5) +
  theme_classic()

sum(Glut_results_1vs0$log2FoldChange < 0)

sum(Glut_results_1vs0$log2FoldChange > 0)

sum(Glut_results_6vs0$log2FoldChange < 0)

sum(Glut_results_6vs0$log2FoldChange > 0)
