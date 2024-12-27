# Chuxuan Li 07/08/2022
# Use the combat-corrected matrix for pseudobulk DE analysis for 18line RNAseq

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

load("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")
RNA_18line <- integrated_labeled
rm(integrated_labeled)
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/renormalized_4celltype_obj.RData")
RNA_5line <- filtered_obj
rm(filtered_obj)

# make raw count matrix ####

makeMat4CombatseqGroup <- function(typeobj, groups) {
  rownames <- rep_len(NA, (length(groups)))
  for (i in 1:length(groups)) {
    print(groups[i])
    rownames[i] <- groups[i]
    obj <- subset(typeobj, orig.ident == groups[i])
    count_matrix <- obj
    if (i == 1) {
      mat <- as.array(rowSums(obj@assays$RNA@counts))
    } else {
      to_bind <- as.array(rowSums(obj@assays$RNA@counts))
      mat <- rbind(mat, to_bind)
    }
  }
  rownames(mat) <- rownames
  return(mat)
}
makeMat4CombatseqLine <- function(typeobj, linetimes) {
  rownames <- rep_len(NA, (length(linetimes)))
  for (i in 1:length(linetimes)) {
    print(linetimes[i])
    rownames[i] <- linetimes[i]
    obj <- subset(typeobj, timexline.ident == linetimes[i])
    count_matrix <- obj 
    if (i == 1) {
      mat <- as.array(rowSums(obj@assays$RNA@counts))
    } else {
      to_bind <- as.array(rowSums(obj@assays$RNA@counts))
      mat <- rbind(mat, to_bind)
    }
  }
  rownames(mat) <- rownames
  return(mat)
}

groups <- unique(RNA_18line$orig.ident)
lts <- unique(RNA_18line$timexline.ident)
GABA_mat <- t(makeMat4CombatseqLine(GABA, lts))
nmglut_mat <- t(makeMat4CombatseqLine(nmglut, lts))
npglut_mat <- t(makeMat4CombatseqLine(npglut, lts))
# batch <- str_split(string = groups, pattern = "-", n = 2, simplify = T)[,1]
# group <- str_split(string = groups, pattern = "-", n = 2, simplify = T)[,2]
batch <- str_split(string = lts, pattern = "_", n = 3, simplify = T)[,2]
group <- str_split(string = lts, pattern = "_", n = 3, simplify = T)[,3]

# combat-seq
GABA_mat_adj <- ComBat_seq(GABA_mat, batch = batch, group = group)
nmglut_mat_adj <- ComBat_seq(nmglut_mat, batch = batch, group = group)
npglut_mat_adj <- ComBat_seq(npglut_mat, batch = batch, group = group)

# DESeq2 ####

# make colData df
coldata_GABA <- data.frame(condition = str_extract(string = colnames(GABA_mat_adj),
                                                   pattern = "[0|1|6]"),
                           # type = str_extract(string = colnames(GABA_mat_adj),
                           #                    pattern = "^[0-9]+-")
                           type = str_extract(string = colnames(GABA_mat_adj),
                                              pattern = "^CD_[0-9][0-9]"))
rownames(coldata_GABA) <- colnames(GABA_mat_adj)

# check
all(rownames(coldata_GABA) %in% colnames(GABA_mat_adj)) # TRUE
all(rownames(coldata_GABA) == colnames(GABA_mat_adj)) #TRUE

# make deseq2 obj
dds_GABA <- DESeqDataSetFromMatrix(countData = GABA_mat_adj,
                                   colData = coldata_GABA,
                                   design = ~ condition)
dds_GABA
p <- hist(rowSums(counts(dds_GABA)), xlim = c(0, 1000), breaks = 100000)
p$counts[1]/sum(p$counts)
p$breaks[2]
keep <- rowSums(counts(dds_GABA)) >= 100
dds_GABA <- dds_GABA[keep,]
# de analysis
dds_GABA <- DESeq(dds_GABA)
res_GABA_0v1 <- results(dds_GABA, contrast =c ("condition", "1", "0"))
res_GABA_0v6 <- results(dds_GABA, contrast =c ("condition", "6", "0"))
res_GABA_1v6 <- results(dds_GABA, contrast =c ("condition", "6", "1"))


# make colData df
coldata_nmglut <- data.frame(condition = str_extract(string = colnames(nmglut_mat_adj),
                                                   pattern = "[0|1|6]"),
                             # type = str_extract(string = colnames(nmglut_mat_adj),
                             #                    pattern = "^[0-9]+-")
                             type = str_extract(string = colnames(nmglut_mat_adj),
                                                pattern = "^CD_[0-9][0-9]"))
rownames(coldata_nmglut) <- colnames(nmglut_mat_adj)

# check
all(rownames(coldata_nmglut) %in% colnames(nmglut_mat_adj)) # TRUE
all(rownames(coldata_nmglut) == colnames(nmglut_mat_adj)) #TRUE

# make deseq2 obj
dds_nmglut <- DESeqDataSetFromMatrix(countData = nmglut_mat_adj,
                                   colData = coldata_nmglut,
                                   design = ~ condition)
dds_nmglut
p <- hist(rowSums(counts(dds_nmglut)), xlim = c(0, 1000), breaks = 100000)
(p$counts[1]+p$counts[2])/sum(p$counts)
p$breaks[2]
keep <- rowSums(counts(dds_nmglut)) >= 50
dds_nmglut <- dds_nmglut[keep,]
# de analysis
dds_nmglut <- DESeq(dds_nmglut)
res_nmglut_0v1 <- results(dds_nmglut, contrast =c ("condition", "1", "0"))
res_nmglut_0v6 <- results(dds_nmglut, contrast =c ("condition", "6", "0"))
res_nmglut_1v6 <- results(dds_nmglut, contrast =c ("condition", "6", "1"))

# make colData df
coldata_npglut <- data.frame(condition = str_extract(string = colnames(npglut_mat_adj),
                                                     pattern = "[0|1|6]"),
                             # type = str_extract(string = colnames(npglut_mat_adj),
                             #                    pattern = "^[0-9]+-")
                             type = str_extract(string = colnames(npglut_mat_adj),
                                                pattern = "^CD_[0-9][0-9]"))
rownames(coldata_npglut) <- colnames(npglut_mat_adj)

# check
all(rownames(coldata_npglut) %in% colnames(npglut_mat_adj)) # TRUE
all(rownames(coldata_npglut) == colnames(npglut_mat_adj)) #TRUE

# make deseq2 obj
dds_npglut <- DESeqDataSetFromMatrix(countData = npglut_mat_adj,
                                     colData = coldata_npglut,
                                     design = ~ condition)
dds_npglut
p <- hist(rowSums(counts(dds_npglut)), xlim = c(0, 1000), breaks = 200000)
(p$counts[1] + p$counts[2] + p$counts[3] + p$counts[4])/sum(p$counts)
p$breaks[5]
keep <- rowSums(counts(dds_npglut)) >= 200
dds_npglut <- dds_npglut[keep,]
# de analysis
dds_npglut <- DESeq(dds_npglut)
res_npglut_0v1 <- results(dds_npglut, contrast =c ("condition", "1", "0"))
res_npglut_0v6 <- results(dds_npglut, contrast =c ("condition", "6", "0"))
res_npglut_1v6 <- results(dds_npglut, contrast =c ("condition", "6", "1"))

save(dds_GABA, dds_nmglut, dds_npglut, file = "pseudobulk_DE_08jul2022_cellline_as_covar_dds.RData")

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


# save results ####
# generate file names
types <- c("GABA", "npglut", "nmglut")
times <- c("1v0", "6v0", "1v6")
filenames <- c()
for (p in types){
  for (m in times){
    filenames <- c(filenames, paste0("./pseudobulk_DE/res_use_combatseq_mat/no_filter_full_df/", 
                                     paste(p, m, "DE_analysis_res_full_dataframe.csv", sep = "_")))
  }
}
filenames

# full results
j = 0
for (i in list(res_GABA_0v1, res_GABA_0v6, res_GABA_1v6,
               res_npglut_0v1, res_npglut_0v6, res_npglut_1v6, 
               res_nmglut_0v1, res_nmglut_0v6, res_nmglut_1v6)) {
  j = j + 1
  print(filenames[j])
  i$gene <- rownames(i)
  write.table(as.data.frame(i), file = filenames[j], quote = F,
              sep = ",", row.names = F, col.names = T)
}

# filter by basemean
res_0v1_list <- list(res_GABA_0v1,
                     res_npglut_0v1,
                     res_nmglut_0v1)
res_0v6_list <- list(res_GABA_0v6,
                     res_npglut_0v6,
                     res_nmglut_0v6)
res_1v6_list <- list(res_GABA_1v6,
                     res_npglut_1v6,
                     res_nmglut_1v6)
hist(res_GABA_0v1$baseMean, breaks = 100000, xlim = c(0, 200))
filter_res <- function(res){
  res <- res[res$baseMean > 10, ]
  res$q_value <- p.adjust(p = res$pvalue, method = "BH")
  res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
  res$pvalue[is.na(res$pvalue)] <- Inf
  res$padj[is.na(res$padj)] <- Inf
  res$q_value[is.na(res$q_value)] <- Inf
  res$gene <- rownames(res)
  res <- as.data.frame(res)
  return(res)
}

reslists <- list(res_0v1_list, res_0v6_list, res_1v6_list)
resorder <- c("1v0", "6v0", "6v1")

for (j in 1:length(reslists)) {
  lst <- reslists[[j]]
  for (i in 1:length(lst)){                                                               
    res <- lst[[i]]
    res <- filter_res(res)
    filename <- paste0("./pseudobulk_DE/res_use_combatseq_mat/filtered_by_basemean_only/", 
                       types[i], "_", resorder[j],
                       "_basemean_filtered_only.csv")
    print(filename)
    write.table(res, file = filename, quote = F, sep = ",", col.names = T, 
                row.names = F)
    filename <- paste0("./pseudobulk_DE/res_use_combatseq_mat/upreg_qval_filtered/", 
                       types[i], "_", resorder[j],
                       "_upregulated_significant.csv")
    print(filename)
    write.table(res[res$q_value < 0.05 & res$log2FoldChange > 0, ],
                file = filename, quote = F, sep = ",", col.names = T, 
                row.names = F)
    filename <- paste0("./pseudobulk_DE/res_use_combatseq_mat/downreg_qval_filtered/", 
                       types[i], "_", resorder[j],
                       "_downregulated_significant.csv")
    print(filename)
    write.table(res[res$q_value < 0.05 & res$log2FoldChange < 0, ],
                file = filename, quote = F, sep = ",", col.names = T, 
                row.names = F)
    }
}
  
# output gene lists only
for (j in 1:length(reslists)) {
  lst <- reslists[[j]]
  for (i in 1:length(lst)){                                                               
    res <- lst[[i]]
    res <- filter_res(res)
    filename <- paste0("./pseudobulk_DE/gene_only/", types[i], "_", resorder[j],
                       "_upregulated_gene_only.txt")
    print(filename)
    write.table(rownames(res[res$q_value < 0.05 & res$log2FoldChange > 0, ]),
                file = filename, quote = F, sep = "\t", col.names = F, row.names = F)
    filename <- paste0("./pseudobulk_DE/gene_only/", types[i], "_", resorder[j], 
                       "_downregulated_gene_only.txt")
    print(filename)
    write.table(rownames(res[res$q_value < 0.05 & res$log2FoldChange < 0, ]),
                file = filename, quote = F, sep = "\t", col.names = F, row.names = F)
  }
}

# bar graph ####
early_list <- c("FOS", "FOSB", "NPAS4", "NR4A1")
late_list <- c("BDNF", "IGF1", "VGF")

df_to_plot <- rbind(res_GABA_0v1[early_list, ],
                    res_GABA_0v6[early_list, ],
                    res_npglut_0v1[early_list, ],
                    res_npglut_0v6[early_list, ],
                    res_nmglut_0v1[early_list, ],
                    res_nmglut_0v6[early_list, ])
df_to_plot <- as.data.frame(df_to_plot)

df_to_plot$time <- c(rep_len("0v1", length(early_list)), rep_len("0v6", length(early_list)), 
                     rep_len("0v1", length(early_list)), rep_len("0v6", length(early_list)),
                     rep_len("0v1", length(early_list)), rep_len("0v6", length(early_list)))                    
df_to_plot$cell.type <- c(rep_len("GABA", length(early_list) * 2), 
                          rep_len("NEFM- glut", length(early_list) * 2),
                          rep_len("NEFM+ glut", length(early_list) * 2))
df_to_plot$gene.name <- c(early_list, early_list, early_list, early_list, early_list, early_list)

# for early, mark if they are excitatory or inhibitory
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
  facet_grid(#labeller = gene_labeller,
             cols = vars(gene.name)) +
  scale_fill_manual(values = c("indianred", "steelblue")) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("differential expression of early response genes")


# summary table ####
cnames <- c("1v0 up", "1v0 down", "6v0 up", "6v0 down")
deg_counts <- array(dim = c(3, 4), dimnames = list(types, cnames)) # 3 cell types are rows, 1/6 up/down are cols
reslists <- list(res_0v1_list, res_0v6_list)
for (j in 1:length(reslists)) {
  lst <- reslists[[j]]
  for (i in 1:length(lst)){                                                               
    res <- lst[[i]]
    res <- filter_res(res)
    deg_counts[i, 2*j-1] <- sum(res$padj < 0.05 & res$log2FoldChange > 0)
    deg_counts[i, 2*j] <- sum(res$padj < 0.05 & res$log2FoldChange < 0)
    #print(sum(res$log2FoldChange < 0))
  }
}
