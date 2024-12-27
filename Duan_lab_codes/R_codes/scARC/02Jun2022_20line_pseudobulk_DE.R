# Chuxuan Li 03/30/2022
# 5-line 20-line combined differential expression analysis using pseudobulk

# init ####
library(Seurat)
library(stringr)
library(readr)
library(ggplot2)
library(RColorBrewer)

# make df ####
load("./demux_20line_integrated_labeled_obj.RData")
line_time <- unique(integrated_labeled$orig.ident)
unique(integrated_labeled$cell.type)
GABA <- subset(integrated_labeled, cell.type %in% c("GABA", "SST_pos_GABA", "SEMA3E_pos_GABA"))
npglut <- subset(integrated_labeled, cell.type == "NEFM_pos_glut")
nmglut <- subset(integrated_labeled, cell.type == "NEFM_neg_glut")

genelist <- rownames(GABA@assays$RNA@counts)
filler <- rep(0, length(genelist))
GABA_pseudobulk <- data.frame(filler, filler, filler, filler, filler, filler, 
                              filler, filler, filler, filler, filler, filler, 
                              filler, filler, filler, row.names = genelist)
colnames(GABA_pseudobulk) <- line_time

genelist <- rownames(npglut@assays$RNA@counts)
filler <- rep(0, length(genelist))
npglut_pseudobulk <- data.frame(filler, filler, filler, filler, filler, filler, 
                                filler, filler, filler, filler, filler, filler, 
                                filler, filler, filler, row.names = genelist)
colnames(npglut_pseudobulk) <- line_time

genelist <- rownames(nmglut@assays$RNA@counts)
filler <- rep(0, length(genelist))
nmglut_pseudobulk <- data.frame(filler, filler, filler, filler, filler, filler, 
                                filler, filler, filler, filler, filler, filler, 
                                filler, filler, filler, row.names = genelist)
colnames(nmglut_pseudobulk) <- line_time

for (i in 1:length(line_time)){
  lt <- line_time[i]
  subobj <- subset(GABA, subset = orig.ident == lt)
  GABA_pseudobulk[, i] <- as.array(rowSums(subobj@assays$RNA@counts))
  subobj <- subset(npglut, subset = orig.ident == lt)
  npglut_pseudobulk[, i] <- as.array(rowSums(subobj@assays$RNA@counts))
  subobj <- subset(nmglut, subset = orig.ident == lt)
  nmglut_pseudobulk[, i] <- as.array(rowSums(subobj@assays$RNA@counts))
}


# DESeq2 ####
library("DESeq2")

# make colData df
coldata_GABA <- data.frame(condition = str_extract(string = colnames(GABA_pseudobulk),
                                                   pattern = "[0|1|6]"),
                           type = str_extract(string = colnames(GABA_pseudobulk),
                                              pattern = "^[0-9]+_"))
rownames(coldata_GABA) <- colnames(GABA_pseudobulk)

# check
all(rownames(coldata_GABA) %in% colnames(GABA_pseudobulk)) # TRUE
all(rownames(coldata_GABA) == colnames(GABA_pseudobulk)) #TRUE

# make deseq2 obj
dds_GABA <- DESeqDataSetFromMatrix(countData = GABA_pseudobulk,
                                   colData = coldata_GABA,
                                   design = ~ type + condition)
dds_GABA
# de analysis
dds_GABA <- DESeq(dds_GABA)
res_GABA_0v1 <- results(dds_GABA, contrast =c ("condition", "1", "0"))
res_GABA_0v6 <- results(dds_GABA, contrast =c ("condition", "6", "0"))
res_GABA_1v6 <- results(dds_GABA, contrast =c ("condition", "6", "1"))


# make colData df
coldata_npglut <- data.frame(condition = str_extract(string = colnames(npglut_pseudobulk),
                                                     pattern = "[0|1|6]"),
                             type = str_extract(string = colnames(npglut_pseudobulk),
                                                pattern = "^[0-9]+_"))
rownames(coldata_npglut) <- colnames(npglut_pseudobulk)

# check
all(rownames(coldata_npglut) %in% colnames(npglut_pseudobulk)) # TRUE
all(rownames(coldata_npglut) == colnames(npglut_pseudobulk)) #TRUE

# make deseq2 obj
dds_npglut <- DESeqDataSetFromMatrix(countData = npglut_pseudobulk,
                                     colData = coldata_npglut,
                                     design = ~ type + condition)
dds_npglut
# de analysis
dds_npglut <- DESeq(dds_npglut)
res_npglut_0v1 <- results(dds_npglut, contrast =c ("condition", "1", "0"))
res_npglut_0v6 <- results(dds_npglut, contrast =c ("condition", "6", "0"))
res_npglut_1v6 <- results(dds_npglut, contrast =c ("condition", "6", "1"))

# make colData df
coldata_nmglut <- data.frame(condition = str_extract(string = colnames(nmglut_pseudobulk),
                                                     pattern = "[0|1|6]"),
                             type = str_extract(string = colnames(nmglut_pseudobulk),
                                                pattern = "^[0-9]+_"))
rownames(coldata_nmglut) <- colnames(nmglut_pseudobulk)

# check
all(rownames(coldata_nmglut) %in% colnames(nmglut_pseudobulk)) # TRUE
all(rownames(coldata_nmglut) == colnames(nmglut_pseudobulk)) #TRUE

# make deseq2 obj
dds_nmglut <- DESeqDataSetFromMatrix(countData = nmglut_pseudobulk,
                                     colData = coldata_nmglut,
                                     design = ~ type + condition)
dds_nmglut
# de analysis
dds_nmglut <- DESeq(dds_nmglut)
res_nmglut_0v1 <- results(dds_nmglut, contrast =c ("condition", "1", "0"))
res_nmglut_0v6 <- results(dds_nmglut, contrast =c ("condition", "6", "0"))
res_nmglut_1v6 <- results(dds_nmglut, contrast =c ("condition", "6", "1"))

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


# save pseudobulk output ####
# # gct file for GSEA
# counts_GABA <- counts(dds_GABA, normalized = T)
# counts_GABA <- as.data.frame(counts_GABA)
# counts_GABA$Description <- 1
# counts_GABA$NAME <- rownames(counts_GABA)
# counts_GABA <- counts_GABA[, c(23, 22, 1:21)]
# write.table(counts_GABA, "./pseudobulk_DE/gct/GABA.txt", 
#             sep = "\t", quote = F, row.names = F, na = "na")
# length(counts_GABA$NAME)
# length(unique(counts_GABA$NAME))
# 
# counts_npglut <- counts(dds_npglut, normalized = T)
# counts_npglut <- as.data.frame(counts_npglut)
# counts_npglut$Description <- rownames(counts_npglut)
# counts_npglut$NAME <- rownames(counts_npglut)
# counts_npglut <- counts_npglut[, c(23, 22, 1:21)]
# write.table(counts_npglut, "./pseudobulk_DE/gct/npglut.gct", sep = "\t", quote = F, row.names = F)
# 
# counts_nmglut <- counts(dds_nmglut, normalized = T)
# counts_nmglut <- as.data.frame(counts_nmglut)
# counts_nmglut$Description <- rownames(counts_nmglut)
# counts_nmglut$NAME <- rownames(counts_nmglut)
# counts_nmglut <- counts_nmglut[, c(23, 22, 1:21)]
# write.table(counts_nmglut, "./pseudobulk_DE/gct/nmglut.gct", sep = "\t", quote = F, row.names = F)
# 

# generate file names
types <- c("GABA", "npglut", "nmglut")
times <- c("1v0", "6v0", "1v6")
filenames <- c()
for (p in types){
  for (m in times){
    filenames <- c(filenames, paste("./pseudobulk_DEresults/", 
                                    p, "_", m, 
                                    "DE_analysis_res_full_dataframe.csv", sep = ""))
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

filter_res <- function(res){
  res <- res[res$baseMean > 10, ]
  res$gene <- rownames(res)
  res$q_value <- p.adjust(p = res$pvalue, method = "BH")
  res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
  res$pvalue[is.na(res$pvalue)] <- Inf
  res$q_value[is.na(res$q_value)] <- Inf
  res <- as.data.frame(res)
  return(res)
}

# save filtered by basemean and by significance up/down files
for (i in 1:length(res_0v1_list)){
  res <- res_0v1_list[[i]]
  res <- filter_res(res)
  print(res$q_value[1])
  
  filename <- paste0("./pseudobulk_DEresults/basemean_filtered_full_res/", types[i], 
                     "_1v0_full_DEG_list_after_basemean_filtering_only.csv")
  print(filename)
  write.table(res, file = filename, quote = F, sep = ",", col.names = T,
              row.names = F)
  filename <- paste0("./pseudobulk_DEresults/significant_res/", types[i], 
                     "_1v0_upregulated_significant_DEG.csv")
  print(filename)
  write.table(res[res$q_value < 0.05 & res$log2FoldChange > 0, ],
              file = filename, quote = F, sep = ",", col.names = T,
              row.names = F)
  filename <- paste0("./pseudobulk_DEresults/significant_res/", types[i], 
                     "_1v0_downregulated_significant_DEG.csv")
  print(filename)
  write.table(res[res$q_value < 0.05 & res$log2FoldChange < 0, ],
              file = filename, quote = F, sep = ",", col.names = T,
              row.names = F)
  
  # 0v6
  res <- res_0v6_list[[i]]
  res <- filter_res(res)
  
  filename <- paste0("./pseudobulk_DEresults/basemean_filtered_full_res/", types[i], 
                     "_6v0_full_DEG_list_after_basemean_filtering_only.csv")
  print(filename)
  write.table(res, file = filename, quote = F, sep = ",", col.names = T,
              row.names = F)
  filename <- paste0("./pseudobulk_DEresults/significant_res/", types[i], 
                     "_6v0_upregulated_significant_DEG.csv")
  print(filename)
  write.table(res[res$q_value < 0.05 & res$log2FoldChange > 0, ],
              file = filename, quote = F, sep = ",", col.names = T,
              row.names = F)
  filename <- paste0("./pseudobulk_DEresults/significant_res/", types[i], 
                     "_6v0_downregulated_significant_DEG.csv")
  print(filename)
  write.table(res[res$q_value < 0.05 & res$log2FoldChange < 0, ],
              file = filename, quote = F, sep = ",", col.names = T,
              row.names = F)
  
  # 1v6
  res <- res_1v6_list[[i]]
  res <- filter_res(res)
  
  filename <- paste0("./pseudobulk_DEresults/basemean_filtered_full_res/", types[i], 
                     "_1v6_full_DEG_list_after_basemean_filtering_only.csv")
  print(filename)
  write.table(res, file = filename, quote = F, sep = ",", col.names = T,
              row.names = F)
  filename <- paste0("./pseudobulk_DEresults/significant_res/", types[i], 
                     "_1v6_upregulated_significant_DEG.csv")
  print(filename)
  write.table(res[res$q_value < 0.05 & res$log2FoldChange > 0, ],
              file = filename, quote = F, sep = ",", col.names = T,
              row.names = F)
  filename <- paste0("./pseudobulk_DEresults/significant_res/", types[i], 
                     "_1v6_downregulated_significant_DEG.csv")
  print(filename)
  write.table(res[res$q_value < 0.05 & res$log2FoldChange < 0, ],
              file = filename, quote = F, sep = ",", col.names = T,
              row.names = F)
}

# output gene lists only
for (i in 1:length(res_0v1_list)){
  res <- res_0v1_list[[i]]
  res <- filter_res(res)
  filename <- paste0(types[i], "_0v1_upregulated_gene_only.txt")
  print(filename)
  write.table(rownames(res[res$q_value < 0.05 & res$log2FoldChange > 0, ]),
              file = filename, quote = F, sep = "\t", col.names = F, row.names = F)
  filename <- paste0(types[i], "_0v1_downregulated_gene_only.txt")
  print(filename)
  write.table(rownames(res[res$q_value < 0.05 & res$log2FoldChange < 0, ]),
              file = filename, quote = F, sep = "\t", col.names = F, row.names = F)
  
  # 0v6
  res <- res_0v6_list[[i]]
  res <- filter_res(res)
  filename <- paste0(types[i], "_0v6_upregulated_gene_only.txt")
  print(filename)
  write.table(rownames(res[res$q_value < 0.05 & res$log2FoldChange > 0, ]),
              file = filename, quote = F, sep = "\t", col.names = F, row.names = F)
  filename <- paste0(types[i], "_0v6_downregulated_gene_only.txt")
  print(filename)
  write.table(rownames(res[res$q_value < 0.05 & res$log2FoldChange < 0, ]),
              file = filename, quote = F, sep = "\t", col.names = F, row.names = F)  
  # 1v6
  res <- res_1v6_list[[i]]
  res <- filter_res(res)
  filename <- paste0(types[i], "_1v6_upregulated_gene_only.txt")
  print(filename)
  write.table(rownames(res[res$q_value < 0.05 & res$log2FoldChange > 0, ]),
              file = filename, quote = F, sep = "\t", col.names = F, row.names = F)
  filename <- paste0(types[i], "_1v6_downregulated_gene_only.txt")
  print(filename)
  write.table(rownames(res[res$q_value < 0.05 & res$log2FoldChange < 0, ]),
              file = filename, quote = F, sep = "\t", col.names = F, row.names = F)
}

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
  facet_grid(#labeller = gene_labeller,
    cols = vars(gene.name)
  ) +
  scale_fill_manual(values = c("indianred", "steelblue")) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("differential expression of late response genes")


# dotplot ####
gene_list <- c("FOS", "FOSB", "NPAS4", "NR4A1")
gene_list <- c("BDNF", "IGF1", "VGF")
df_to_plot <- rbind(res_GABA_0v1[gene_list, ],
                    res_GABA_0v6[gene_list, ],
                    res_npglut_0v1[gene_list, ],
                    res_npglut_0v6[gene_list, ],
                    res_nmglut_0v1[gene_list, ],
                    res_nmglut_0v6[gene_list, ])
df_to_plot <- as.data.frame(df_to_plot)

df_to_plot$time <- factor(c(rep_len("0v1", length(gene_list)), rep_len("0v6", length(gene_list)), 
                            rep_len("0v1", length(gene_list)), rep_len("0v6", length(gene_list)),
                            rep_len("0v1", length(gene_list)), rep_len("0v6", length(gene_list))),
                          levels = c("0v1","0v6"))
df_to_plot$cell.type <- c(rep_len("GABA", length(gene_list) * 2), 
                          rep_len("NEFM- glut", length(gene_list) * 2),
                          rep_len("NEFM+ glut", length(gene_list) * 2))
df_to_plot$gene.name <- c(gene_list, gene_list, gene_list, gene_list, gene_list, gene_list)

ggplot(df_to_plot, mapping = aes(x = time, y = gene.name, 
                                 size = -log(pvalue))) +
  geom_point(mapping = aes(color = log2FoldChange)) +
  facet_grid(cols = vars(cell.type)) + 
  scale_color_gradientn(colors = brewer.pal(5, "Reds")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 8)) +
  ylab("Gene")

