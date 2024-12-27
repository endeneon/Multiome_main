# Chuxuan Li 01/28/2022
# Pseudobulk analysis on 20 line RNAseq data, using clustering from 10x

# init ####
library(Seurat)
library(stringr)
library(readr)
library(ggplot2)
library(RColorBrewer)

# load data ####
# read in the combined object obtained by the -aggr method from 10x
aggr_raw <- Read10X_h5(filename = "/home/cli/Data/FASTQ/Duan_Project_022_Reseq/Duan_022_reseq_aggregated/outs/filtered_feature_bc_matrix.h5")

aggr_gex <- CreateSeuratObject(counts = aggr_raw$`Gene Expression`)
rm(aggr_raw)

aggr_gex$time.group.ident <- str_extract(string = aggr_gex@assays$RNA@data@Dimnames[[2]], 
                                         pattern = "[0-9]+$")
unique(aggr_gex$time.group.ident)

# assign idents ####
aggr_gex$time.ident <- NA
aggr_gex$time.ident[aggr_gex$time.group.ident %in% as.character(seq(1, 13, 3))] <- "0hr"
aggr_gex$time.ident[aggr_gex$time.group.ident %in% as.character(seq(2, 14, 3))] <- "1hr"  
aggr_gex$time.ident[aggr_gex$time.group.ident %in% as.character(seq(3, 15, 3))] <- "6hr"
unique(aggr_gex$time.ident)


cluster_idents <- 
  read_csv("/home/cli/Data/FASTQ/Duan_Project_022_Reseq/Duan_022_reseq_aggregated/outs/analysis/clustering/gex/graphclust/clusters.csv",
  col_names = T)
         
aggr_gex$tenx.cluster <- NA
# test with one
aggr_gex$tenx.cluster[aggr_gex@assays$RNA@counts@Dimnames[[2]] == 
                        cluster_idents$Barcode[1]] <- cluster_idents$Cluster[1]
unique(aggr_gex$tenx.cluster)
# loop
for (i in 1:length(cluster_idents$Barcode)){
  aggr_gex$tenx.cluster[aggr_gex@assays$RNA@counts@Dimnames[[2]] == 
                          cluster_idents$Barcode[i]] <- cluster_idents$Cluster[i]
  
}
unique(aggr_gex$tenx.cluster)

aggr_gex$cell.type <- "others"
aggr_gex$cell.type[aggr_gex$tenx.cluster %in% c(1, 2, 3, 4, 17, 9, 19)] <- "GABA"
aggr_gex$cell.type[aggr_gex$tenx.cluster %in% c(7, 5, 14, 18)] <- "NEFM_pos_glut"
aggr_gex$cell.type[aggr_gex$tenx.cluster %in% c(7, 5, 14, 18)] <- "NEFM_neg_glut"
unique(aggr_gex$cell.type)

GABA <- subset(aggr_gex, subset = cell.type == "GABA")
glut <- subset(aggr_gex, subset = cell.type == "NEFM_neg_glut")

line_time <- as.character(sort(as.numeric(unique(aggr_gex$time.group.ident))))


genelist <- rownames(GABA@assays$RNA@data)
filler <- rep(0, length(genelist))
GABA_pseudobulk <- data.frame(lib22_0hr = filler, lib22_1hr = filler, lib22_6hr = filler,
                              lib36_0hr = filler, lib36_1hr = filler, lib36_6hr = filler,
                              lib39_0hr = filler, lib39_1hr = filler, lib39_6hr = filler,
                              lib44_0hr = filler, lib44_1hr = filler, lib44_6hr = filler,
                              lib49_0hr = filler, lib49_1hr = filler, lib49_6hr = filler, row.names = genelist)
genelist <- rownames(glut@assays$RNA@data)
filler <- rep(0, length(genelist))
glut_pseudobulk <- data.frame(lib22_0hr = filler, lib22_1hr = filler, lib22_6hr = filler,
                              lib36_0hr = filler, lib36_1hr = filler, lib36_6hr = filler,
                              lib39_0hr = filler, lib39_1hr = filler, lib39_6hr = filler,
                              lib44_0hr = filler, lib44_1hr = filler, lib44_6hr = filler,
                              lib49_0hr = filler, lib49_1hr = filler, lib49_6hr = filler, row.names = genelist)


for (i in 1:length(line_time)){
  lt <- line_time[i]
  subobj <- subset(GABA, subset = time.group.ident == lt)
  GABA_pseudobulk[, i] <- as.array(rowSums(subobj@assays$RNA@data))
  subobj <- subset(glut, subset = time.group.ident == lt)
  glut_pseudobulk[, i] <- as.array(rowSums(subobj@assays$RNA@data))
}


# DESeq2 ####
library("DESeq2")

# make colData df
coldata_GABA <- data.frame(condition = str_extract(string = colnames(GABA_pseudobulk),
                                                   pattern = "[0-6]hr"),
                           type = str_extract(string = colnames(GABA_pseudobulk),
                                              pattern = "lib[0-9][0-9]"))
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
res_GABA_0v1 <- results(dds_GABA, contrast =c ("condition", "1hr", "0hr"))
res_GABA_0v6 <- results(dds_GABA, contrast =c ("condition", "6hr", "0hr"))
res_GABA_1v6 <- results(dds_GABA, contrast =c ("condition", "6hr", "1hr"))


# make colData df
coldata_glut <- data.frame(condition = str_extract(string = colnames(glut_pseudobulk),
                                                   pattern = "[0-6]hr"),
                           type = str_extract(string = colnames(glut_pseudobulk),
                                              pattern = "lib[0-9][0-9]"))
rownames(coldata_glut) <- colnames(glut_pseudobulk)

# check
all(rownames(coldata_glut) %in% colnames(glut_pseudobulk)) # TRUE
all(rownames(coldata_glut) == colnames(glut_pseudobulk)) #TRUE

# make deseq2 obj
dds_glut <- DESeqDataSetFromMatrix(countData = glut_pseudobulk,
                                   colData = coldata_glut,
                                   design = ~ type + condition)
dds_glut
# de analysis
dds_glut <- DESeq(dds_glut)
res_glut_0v1 <- results(dds_glut, contrast =c ("condition", "1hr", "0hr"))
res_glut_0v6 <- results(dds_glut, contrast =c ("condition", "6hr", "0hr"))
res_glut_1v6 <- results(dds_glut, contrast =c ("condition", "6hr", "1hr"))


# look at specific genes
genelist <- c("GAPDH", "BDNF", "FOS", "VGF", "NRN1")
res_glut_0v1[genelist, ]
res_glut_0v6[genelist, ]
res_glut_1v6[genelist, ]

res_GABA_0v1[genelist, ]
res_GABA_0v6[genelist, ]
res_GABA_1v6[genelist, ]


# edgeR ####

library(edgeR)

Glut_DGEList <- DGEList(counts = as.matrix(glut_pseudobulk),
                        genes = rownames(glut_pseudobulk), 
                        samples = colnames(glut_pseudobulk),
                        remove.zeros = T)

Glut_DGEList <- calcNormFactors(Glut_DGEList)
Glut_DGEList <- estimateDisp(Glut_DGEList)

hist(cpm(Glut_DGEList), breaks = 10000, xlim = c(0, 10))
sum(cpm(Glut_DGEList) < 1)

# filter DGEList
Glut_DGEList <- Glut_DGEList[rowSums(cpm(Glut_DGEList) > 1) >= 5, ,
                             keep.lib.sizes = F]
# Glut_DGEList$genes

Glut_design_matrix <- coldata_glut
Glut_design_matrix$condition <- as.factor(Glut_design_matrix$condition)
Glut_design_matrix$type <- as.factor(Glut_design_matrix$type)

Glut_design <- model.matrix(~ 0 + 
                              Glut_design_matrix$condition +
                              Glut_design_matrix$type)

Glut_fit <- glmQLFit(Glut_DGEList,
                     design = Glut_design)
Glut_1vs0 <- glmQLFTest(Glut_fit, 
                        contrast = c(-1, 1, 0, 0, 0, 0, 0))
Glut_6vs0 <- glmQLFTest(Glut_fit, 
                        contrast = c(-1, 0, 1, 0, 0, 0, 0))
Glut_6vs1 <- glmQLFTest(Glut_fit, 
                        contrast = c(0, -1, 1, 0, 0, 0, 0))

Glut_1vs0 <- Glut_1vs0$table
Glut_6vs0 <- Glut_6vs0$table
Glut_6vs1 <- Glut_6vs1$table


Glut_1vs0[rownames(Glut_1vs0) %in% genelist, ]
Glut_6vs0[rownames(Glut_6vs0) %in% genelist, ]
Glut_6vs1[rownames(Glut_6vs1) %in% genelist, ]



GABA_DGEList <- DGEList(counts = as.matrix(GABA_pseudobulk),
                        genes = rownames(GABA_pseudobulk), 
                        samples = colnames(GABA_pseudobulk),
                        remove.zeros = T)

GABA_DGEList <- calcNormFactors(GABA_DGEList)
GABA_DGEList <- estimateDisp(GABA_DGEList)

hist(cpm(GABA_DGEList), breaks = 10000, xlim = c(0, 10))
sum(cpm(GABA_DGEList) < 1)

# filter DGEList
GABA_DGEList <- GABA_DGEList[rowSums(cpm(GABA_DGEList) > 1) >= 5, ,
                             keep.lib.sizes = F]
# GABA_DGEList$genes

GABA_design_matrix <- coldata_GABA
GABA_design_matrix$condition <- as.factor(GABA_design_matrix$condition)
GABA_design_matrix$type <- as.factor(GABA_design_matrix$type)

GABA_design <- model.matrix(~ 0 + 
                              GABA_design_matrix$condition +
                              GABA_design_matrix$type)

GABA_fit <- glmQLFit(GABA_DGEList,
                     design = GABA_design)
GABA_1vs0 <- glmQLFTest(GABA_fit, 
                        contrast = c(-1, 1, 0, 0, 0, 0, 0))
GABA_6vs0 <- glmQLFTest(GABA_fit, 
                        contrast = c(-1, 0, 1, 0, 0, 0, 0))
GABA_6vs1 <- glmQLFTest(GABA_fit, 
                        contrast = c(0, -1, 1, 0, 0, 0, 0))

GABA_1vs0 <- GABA_1vs0$table
GABA_6vs0 <- GABA_6vs0$table
GABA_6vs1 <- GABA_6vs1$table


GABA_1vs0[rownames(GABA_1vs0) %in% genelist, ]
GABA_6vs0[rownames(GABA_6vs0) %in% genelist, ]
GABA_6vs1[rownames(GABA_6vs1) %in% genelist, ]


# bar graph ####
genelist <- genelist[genelist != c("GAPDH")]
df_to_plot <- rbind(res_glut_0v1[genelist, ],
                    res_glut_0v6[genelist, ],
                    res_GABA_0v1[genelist, ],
                    res_GABA_0v6[genelist, ])
df_to_plot$time <- c(rep_len("0v1", 3), rep_len("0v6", 3), rep_len("0v1", 3), rep_len("0v6", 3))                    
df_to_plot$cell.type <- c(rep_len("glut", 6), rep_len("GABA", 6))
df_to_plot$gene.name <- c(genelist, genelist)
df_to_plot <- as.data.frame(df_to_plot)

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
  facet_grid(cols = vars(gene.name)) +
  scale_fill_manual(values = c("indianred", "steelblue")) +
  theme_bw() +
  theme(legend.title = element_blank())


df_to_plot$log2FoldChange - 1/2*df_to_plot$lfcSE
