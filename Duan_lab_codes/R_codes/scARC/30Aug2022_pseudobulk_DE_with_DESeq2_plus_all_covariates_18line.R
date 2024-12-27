# Chuxuan Li 08/30/2022
# use DESeq2 + more covariates (same as the ones used in Limma DE analysis),
#redo 18line DE analysis

# init ####
library(DESeq2)
library(readr)
library(stringr)
library(ggplot2)
# load combat-adjusted data frames for each cell type
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/pseudobulk_DE/res_use_combatseq_mat/by_line_pseudobulk_mat_adjusted.RData")

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
