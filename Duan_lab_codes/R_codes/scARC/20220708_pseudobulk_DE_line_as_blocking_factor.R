# Chuxuan Li 04/04/2022
# pseudobulk DE analysis on 18-line RNAseq data

# init ####
library(Seurat)
library(DESeq2)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(future)

# load data ####
load("./labeled_nfeat3000.RData")
line_time <- unique(RNA_18line$timexline.ident)
unique(RNA_18line$cell.type)
GABA <- subset(RNA_18line, cell.type %in% c("GABA", "SST_pos_GABA", "SEMA3E_pos_GABA"))
npglut <- subset(RNA_18line, cell.type == "NEFM_pos_glut")
nmglut <- subset(RNA_18line, cell.type == "NEFM_neg_glut")

genelist <- rownames(GABA@assays$RNA@counts)
GABA_pseudobulk <- as.data.frame(array(dim = c(length(genelist), length(line_time)), 
                                       dimnames = list(genelist, line_time)))

genelist <- rownames(npglut@assays$RNA@counts)
npglut_pseudobulk <- as.data.frame(array(dim = c(length(genelist), length(line_time)), 
                                         dimnames = list(genelist, line_time)))

genelist <- rownames(nmglut@assays$RNA@counts)
nmglut_pseudobulk <- as.data.frame(array(dim = c(length(genelist), length(line_time)), 
                                         dimnames = list(genelist, line_time)))

for (i in 1:length(line_time)){
  lt <- line_time[i]
  subobj <- subset(GABA, subset = timexline.ident == lt)
  GABA_pseudobulk[, i] <- as.array(rowSums(subobj@assays$RNA@counts))
  subobj <- subset(npglut, subset = timexline.ident == lt)
  npglut_pseudobulk[, i] <- as.array(rowSums(subobj@assays$RNA@counts))
  subobj <- subset(nmglut, subset = timexline.ident == lt)
  nmglut_pseudobulk[, i] <- as.array(rowSums(subobj@assays$RNA@counts))
}


# DESeq2 ####
# make colData df
coldata_GABA <- data.frame(condition = str_extract(string = colnames(GABA_pseudobulk),
                                                   pattern = "[0|1|6]"),
                           type = str_extract(string = colnames(GABA_pseudobulk),
                                              pattern = "^CD_[0-9][0-9]"))
rownames(coldata_GABA) <- colnames(GABA_pseudobulk)

# check
all(rownames(coldata_GABA) %in% colnames(GABA_pseudobulk)) # TRUE
all(rownames(coldata_GABA) == colnames(GABA_pseudobulk)) #TRUE

# make deseq2 obj
dds_GABA <- DESeqDataSetFromMatrix(countData = GABA_pseudobulk,
                                   colData = coldata_GABA,
                                   design = ~ type + condition)
dds_GABA
p <- hist(rowSums(counts(dds_GABA)), breaks = 100000, xlim = c(0, 10000))
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
coldata_npglut <- data.frame(condition = str_extract(string = colnames(npglut_pseudobulk),
                                                     pattern = "[0|1|6]"),
                             type = str_extract(string = colnames(npglut_pseudobulk),
                                                pattern = "^CD_[0-9][0-9]"))
rownames(coldata_npglut) <- colnames(npglut_pseudobulk)

# check
all(rownames(coldata_npglut) %in% colnames(npglut_pseudobulk)) # TRUE
all(rownames(coldata_npglut) == colnames(npglut_pseudobulk)) #TRUE

# make deseq2 obj
dds_npglut <- DESeqDataSetFromMatrix(countData = npglut_pseudobulk,
                                     colData = coldata_npglut,
                                     design = ~ type + condition)
dds_npglut
keep <- rowSums(counts(dds_npglut)) >= 100
dds_npglut <- dds_npglut[keep,]
# de analysis
dds_npglut <- DESeq(dds_npglut)
res_npglut_0v1 <- results(dds_npglut, contrast =c ("condition", "1", "0"))
res_npglut_0v6 <- results(dds_npglut, contrast =c ("condition", "6", "0"))
res_npglut_1v6 <- results(dds_npglut, contrast =c ("condition", "6", "1"))

# make colData df
coldata_nmglut <- data.frame(condition = str_extract(string = colnames(nmglut_pseudobulk),
                                                     pattern = "[0|1|6]"),
                             type = str_extract(string = colnames(nmglut_pseudobulk),
                                                pattern = "^CD_[0-9][0-9]"))
rownames(coldata_nmglut) <- colnames(nmglut_pseudobulk)

# check
all(rownames(coldata_nmglut) %in% colnames(nmglut_pseudobulk)) # TRUE
all(rownames(coldata_nmglut) == colnames(nmglut_pseudobulk)) #TRUE

# make deseq2 obj
dds_nmglut <- DESeqDataSetFromMatrix(countData = nmglut_pseudobulk,
                                     colData = coldata_nmglut,
                                     design = ~ type + condition)
dds_nmglut
p <- hist(rowSums(counts(dds_nmglut)), breaks = 100000, xlim = c(0, 1000))
keep <- rowSums(counts(dds_nmglut)) >= 100
dds_nmglut <- dds_nmglut[keep,]
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
  res$q_value[is.na(res$q_value)] <- Inf
  res$padj[is.na(res$padj)] <- Inf
  res$gene <- rownames(res)
  res <- as.data.frame(res)
  return(res)
}

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
    print(sum(res$log2FoldChange < 0))
  }
}
write.table(deg_counts, file = "./pseudobulk_DE/deg_summary_deseq_line_as_cov_min_count_200.csv", quote = F, sep = ",", row.names = T, col.names = T)
