# Chuxuan Li 09/07/2022
# Use 18-line combat-corrected pseudobulk count matrices and DESeq2, adding
#covariates to do DE analysis

# init ####
library(DESeq2)
library(Seurat)
library(sva)

library(stringr)
library(future)

library(ggplot2)
library(RColorBrewer)
library(corrplot)
library(ggrepel)
library(reshape2)

load("../../../covariates_pooled_together_for_all_linextime_5+18+20.RData")
load("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")
RNA_18line <- integrated_labeled
rm(integrated_labeled)
unique(RNA_18line$cell.type)

GABA <- subset(RNA_18line, cell.type %in% c("GABA", "SST_pos_GABA", "SEMA3E_pos_GABA"))
npglut <- subset(RNA_18line, cell.type == "NEFM_pos_glut")
nmglut <- subset(RNA_18line, cell.type == "NEFM_neg_glut")

# make raw count matrix ####
makeMat4CombatseqLine <- function(typeobj, linetimes) {
  colnames <- rep_len(NA, (length(linetimes)))
  for (i in 1:length(linetimes)) {
    print(linetimes[i])
    colnames[i] <- linetimes[i]
    obj <- subset(typeobj, timexline.ident == linetimes[i])
    count_matrix <- obj 
    if (i == 1) {
      mat <- as.array(rowSums(obj@assays$RNA@counts))
    } else {
      to_bind <- as.array(rowSums(obj@assays$RNA@counts))
      mat <- cbind(mat, to_bind)
    }
  }
  colnames(mat) <- colnames
  return(mat)
}

groups <- sort(unique(RNA_18line$orig.ident))
lts <- sort(unique(RNA_18line$timexline.ident))
GABA_mat <- makeMat4CombatseqLine(GABA, lts)
nmglut_mat <- makeMat4CombatseqLine(nmglut, lts)
npglut_mat <- makeMat4CombatseqLine(npglut, lts)

# extract batch information directly from the covar table
celllines <- str_remove(colnames(GABA_mat), "_[0|1|6]hr$")

for (i in 1:length(celllines)) {
  print(celllines[i])
  if (i == 1) {
    batch = unique(covar_table$group[covar_table$cell_line == celllines[i]])
  } else {
    batch = c(batch, unique(covar_table$group[covar_table$cell_line == celllines[i]]))
  }
}
batch <- as.factor(batch)

# combat-seq
GABA_mat_adj <- ComBat_seq(GABA_mat, batch = batch, group = NULL)
nmglut_mat_adj <- ComBat_seq(nmglut_mat, batch = batch,group = NULL)
npglut_mat_adj <- ComBat_seq(npglut_mat, batch = batch,group = NULL)

# make colData df####
for (i in 1:length(celllines)) {
  print(celllines[i])
  if (i == 1) {
    age = unique(covar_table$age[covar_table$cell_line == celllines[i]])
    sex = unique(covar_table$sex[covar_table$cell_line == celllines[i]])
    disease = unique(covar_table$disease[covar_table$cell_line == celllines[i]])
  } else {
    age = c(age, unique(covar_table$age[covar_table$cell_line == celllines[i]]))
    sex = c(sex, unique(covar_table$sex[covar_table$cell_line == celllines[i]]))
    disease = c(disease, unique(covar_table$disease[covar_table$cell_line == celllines[i]]))
  }
}
for (i in 1:length(colnames(GABA_mat))) {
  l <- str_remove(colnames(GABA_mat)[i], "_[0|1|6]hr$")
  t <- str_extract(colnames(GABA_mat)[i], "[0|1|6]hr$")
  print(paste(l, t))
  if (i == 1) {
    GABA_fraction = covar_table$GABA_fraction[covar_table$cell_line == l & covar_table$time == t]
    nmglut_fraction = covar_table$nmglut_fraction[covar_table$cell_line == l & covar_table$time == t]
    npglut_fraction = covar_table$npglut_fraction[covar_table$cell_line == l & covar_table$time == t]
  } else {
    GABA_fraction = c(GABA_fraction, covar_table$GABA_fraction[covar_table$cell_line == l & covar_table$time == t])
    nmglut_fraction = c(nmglut_fraction, covar_table$nmglut_fraction[covar_table$cell_line == l & covar_table$time == t])
    npglut_fraction = c(npglut_fraction, covar_table$npglut_fraction[covar_table$cell_line == l & covar_table$time == t])
  }
}

coldata_GABA <- data.frame(time = str_extract(string = colnames(GABA_mat_adj),
                                                   pattern = "[0|1|6]"),
                           batch = batch,
                           age = age,
                           aff = disease,
                           sex = sex,
                           GABA_fraction = GABA_fraction)
rownames(coldata_GABA) <- colnames(GABA_mat_adj)

# check
all(rownames(coldata_GABA) %in% colnames(GABA_mat_adj)) # TRUE
all(rownames(coldata_GABA) == colnames(GABA_mat_adj)) #TRUE

# make deseq2 obj
dds_GABA <- DESeqDataSetFromMatrix(countData = GABA_mat_adj,
                                   colData = coldata_GABA,
                                   design = ~ time + batch + age + aff + sex + GABA_fraction)
dds_GABA
p <- hist(rowSums(counts(dds_GABA)), xlim = c(0, 200), breaks = 1000000)
p$counts[1]/sum(p$counts)
p$breaks[2]
keep <- rowSums(counts(dds_GABA)) > 200
dds_GABA <- dds_GABA[keep,]

# de analysis
dds_GABA <- DESeq(dds_GABA)
res_GABA_0v1 <- results(dds_GABA, contrast = c("time", "1", "0"), independentFiltering = F)
res_GABA_0v6 <- results(dds_GABA, contrast = c("time", "6", "0"), independentFiltering = F)
res_GABA_1v6 <- results(dds_GABA, contrast = c("time", "6", "1"), independentFiltering = F)


# make colData df
coldata_nmglut <- data.frame(time = str_extract(string = colnames(GABA_mat_adj),
                                                pattern = "[0|1|6]"),
                             batch = batch,
                             age = age,
                             aff = disease,
                             sex = sex,
                             nmglut_fraction = nmglut_fraction)
rownames(coldata_nmglut) <- colnames(nmglut_mat_adj)

# check
all(rownames(coldata_nmglut) %in% colnames(nmglut_mat_adj)) # TRUE
all(rownames(coldata_nmglut) == colnames(nmglut_mat_adj)) #TRUE

# make deseq2 obj
dds_nmglut <- DESeqDataSetFromMatrix(countData = nmglut_mat_adj,
                                     colData = coldata_nmglut,
                                     design = ~ time + batch + age + aff + sex + nmglut_fraction)
dds_nmglut
keep <- rowSums(counts(dds_nmglut)) > 200
dds_nmglut <- dds_nmglut[keep,]
# de analysis
dds_nmglut <- DESeq(dds_nmglut)
res_nmglut_0v1 <- results(dds_nmglut, contrast = c("time", "1", "0"), independentFiltering = F)
res_nmglut_0v6 <- results(dds_nmglut, contrast = c("time", "6", "0"), independentFiltering = F)
res_nmglut_1v6 <- results(dds_nmglut, contrast = c("time", "6", "1"), independentFiltering = F)

# make colData df
coldata_npglut <- data.frame(time = str_extract(string = colnames(GABA_mat_adj),
                                                pattern = "[0|1|6]"),
                             batch = batch,
                             age = age,
                             aff = disease,
                             sex = sex,
                             npglut_fraction = npglut_fraction)
rownames(coldata_npglut) <- colnames(npglut_mat_adj)

# check
all(rownames(coldata_npglut) %in% colnames(npglut_mat_adj)) # TRUE
all(rownames(coldata_npglut) == colnames(npglut_mat_adj)) #TRUE

# make deseq2 obj
dds_npglut <- DESeqDataSetFromMatrix(countData = npglut_mat_adj,
                                     colData = coldata_npglut,
                                     design = ~ time + batch + age + aff + sex + npglut_fraction)
dds_npglut
p <- hist(rowSums(counts(dds_npglut)), xlim = c(0, 1000), breaks = 200000)
keep <- rowSums(counts(dds_npglut)) > 200
dds_npglut <- dds_npglut[keep,]

# de analysis
dds_npglut <- DESeq(dds_npglut)
res_npglut_0v1 <- results(dds_npglut, contrast =c ("time", "1", "0"), independentFiltering = F)
res_npglut_0v6 <- results(dds_npglut, contrast =c ("time", "6", "0"), independentFiltering = F)
res_npglut_1v6 <- results(dds_npglut, contrast =c ("time", "6", "1"), independentFiltering = F)

# look at specific genes
genelist <- c("BDNF", "IGF1", "VGF", "FOS", "FOSB", "NPAS4", "NR4A1")
res_GABA_0v1[genelist, ]
res_GABA_0v6[genelist, ]
res_GABA_1v6[genelist, ]
res_npglut_0v1[genelist, ]
res_npglut_0v6[genelist, ]
res_npglut_1v6[genelist, ]
res_nmglut_0v1[genelist, ]
res_nmglut_0v6[genelist, ]
res_nmglut_1v6[genelist, ]

sum(res_nmglut_0v1$pvalue < 0.05, na.rm = T)
hist(res_nmglut_0v1$padj, breaks = 1000)

# bar graph ####
early_list <- c("FOS", "FOSB", "NPAS4", "NR4A1")
late_list <- c("BDNF", "IGF1", "VGF")

df_to_plot <- rbind(res_GABA_0v1[late_list, ],
                    res_GABA_0v6[late_list, ],
                    res_npglut_0v1[late_list, ],
                    res_npglut_0v6[late_list, ],
                    res_nmglut_0v1[late_list, ],
                    res_nmglut_0v6[late_list, ])
df_to_plot <- as.data.frame(df_to_plot)

df_to_plot$time <- c(rep_len("0v1", length(late_list)), rep_len("0v6", length(late_list)), 
                     rep_len("0v1", length(late_list)), rep_len("0v6", length(late_list)),
                     rep_len("0v1", length(late_list)), rep_len("0v6", length(late_list)))                    
df_to_plot$cell.type <- c(rep_len("GABA", length(late_list) * 2), 
                          rep_len("NEFM- glut", length(late_list) * 2),
                          rep_len("NEFM+ glut", length(late_list) * 2))
df_to_plot$gene.name <- c(late_list, late_list, late_list, late_list, late_list, late_list)

# for late, mark if they are excitatory or inhibitory
gene_labeller <- as_labeller(c('BDNF' = "BDNF (excitatory)", 
                               'IGF1' = "IGF1 (inhibitory)", 
                               'VGF' = "VGF (shared)"))
ggplot(df_to_plot, aes(x = cell.type, 
                       y = log2FoldChange,
                       color = cell.type,
                       fill = time,
                       group = time,
                       ymax = log2FoldChange - 1/2*lfcSE, 
                       ymin = log2FoldChange + 1/2*lfcSE
)) + 
  xlab("") +
  ylab("log2(FC)") +
  geom_col(position = position_dodge(0.6),
           width = 0.5,
           color = "black") +
  geom_errorbar(color = "black",
                position = position_dodge(0.6),
                width = 0.5) +
  facet_grid(labeller = gene_labeller,
    cols = vars(gene.name)) +
  scale_fill_manual(values = c("indianred", "steelblue")) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("differential expression of late response genes")


# summary table ####
cnames <- c("1v0 up", "1v0 down", "6v0 up", "6v0 down")
types <- c("GABA", "nmglut", "npglut")
res_0v1_list <- list(res_GABA_0v1, res_nmglut_0v1, res_npglut_0v1)
res_0v6_list <- list(res_GABA_0v6, res_nmglut_0v6, res_npglut_0v6)
reslists <- list(res_0v1_list, res_0v6_list)
deg_counts <- array(dim = c(6, 5), 
                    dimnames = list(rep(types, times = 2),
                                    c("up", "down", "nonsig", "total", "NA"))) 
for (j in 1:length(reslists)) {
  lst <- reslists[[j]]
  for (i in 1:length(lst)){                                                               
    res <- lst[[i]]
    deg_counts[j * 3 - 3 + i, 1] <- sum(res$log2FoldChange > 0 & res$padj < 0.05, na.rm = T)
    deg_counts[j * 3 - 3 + i, 2] <- sum(res$log2FoldChange < 0 & res$padj < 0.05, na.rm = T)     
    deg_counts[j * 3 - 3 + i, 3] <- sum(res$padj > 0.05, na.rm = T)
    deg_counts[j * 3 - 3 + i, 4] <- nrow(res)
    deg_counts[j * 3 - 3 + i, 5] <- sum(is.na(res$padj))
  }
}
write.table(deg_counts, file = "./pseudobulk_DE/res_use_combatseq_mat/redo_DESeq2_summary_table.csv", quote = F, sep = ",", 
            row.names = T, col.names = T)
