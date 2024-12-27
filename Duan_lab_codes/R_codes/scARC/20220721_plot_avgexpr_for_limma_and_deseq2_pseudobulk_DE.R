# Chuxuan Li 07/21/2022
# plot the average expression level in the intermediate steps in Limma and Deseq2
# Limma uses combat adjusted counts and filter by cpm>1 in 9 samples in the 0hr
#samples

# init ####
library(Seurat)
library(DESeq2)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(readr)
load("18line_pseudobulk_limma_post_normfactor_voom_dfs.RData")
load("./pseudobulk_DE/res_use_combatseq_mat/by_line_pseudobulk_mat_adjusted.RData")
genes <- c("SGCD", "IMMP2L", "SNAP91", "BCL11B")

# cpm for limma ####
cexprs <- cpm(npglut_mat_adj, log = T)
cexprs <- cexprs[rownames(cexprs) %in% genes, ]

# Limma ####
lexprs <- v_npglut$E
lexprs <- lexprs[rownames(lexprs) %in% genes, ]

# deseq2 after adding covariates ####
dexprs <- dds_npglut@assays@data@listData$mu
rownames(dexprs) <- rownames(dds_npglut@assays@data@listData$counts)
colnames(dexprs) <- colnames(dds_npglut@assays@data@listData$counts)
dexprs <- dexprs[rownames(dexprs) %in% genes, ]
dexprs <- dexprs[, order(colnames(dexprs))]

# deseq2 using only cell line as covariate ####
load("pseudobulk_DE_08jul2022_cellline_as_covar_dds.RData")
dexprs0 <- dds_npglut@assays@data@listData$mu
rownames(dexprs0) <- rownames(dds_npglut@assays@data@listData$counts)
colnames(dexprs0) <- colnames(dds_npglut@assays@data@listData$counts)
dexprs0 <- dexprs0[rownames(dexprs0) %in% genes, ]
dexprs0 <- dexprs0[, order(colnames(dexprs0))]

# plot three together ####
ldf <- matrix(nrow = 54*4, ncol = 4, dimnames = list(rep_len("", 54*4), 
                                                     c("Time", "Sample", "Gene", 
                                                           "Expression")))
cdf <- matrix(nrow = 54*4, ncol = 4, dimnames = list(rep_len("", 54*4),
                                                     c("Time", "Sample", "Gene", 
                                                       "Expression")))
ddf <- matrix(nrow = 54*4, ncol = 4, dimnames = list(rep_len("", 54*4),
                                                     c("Time", "Sample", "Gene", 
                                                       "Expression")))
ddf0 <- matrix(nrow = 54*4, ncol = 4, dimnames = list(rep_len("", 54*4),
                                                      c("Time", "Sample", "Gene", 
                                                        "Expression")))
cnames <- (colnames(dexprs))
times <- str_extract(cnames, "[0|1|6]hr")
lines <- str_extract(cnames, "CD_[0-9]+")

for (i in 1:length(cnames)) {
  start = i * 4 - 3
  end = i * 4
  print(paste(start, end, sep = "-"))
  ldf[start:end, 1] <- rep_len(times[i], 4)
  ldf[start:end, 2] <- rep_len(lines[i], 4)
  ldf[start:end, 3] <- as.array(names(lexprs[,i]))
  ldf[start:end, 4] <- as.numeric(lexprs[,i])
  cdf[start:end, 1] <- rep_len(times[i], 4)
  cdf[start:end, 2] <- rep_len(lines[i], 4)
  cdf[start:end, 3] <- as.array(names(cexprs[,i]))
  cdf[start:end, 4] <- as.numeric(cexprs[,i])
  ddf[start:end, 1] <- rep_len(times[i], 4)
  ddf[start:end, 2] <- rep_len(lines[i], 4)
  ddf[start:end, 3] <- as.array(names(dexprs[,i]))
  ddf[start:end, 4] <- as.numeric(dexprs[,i])
  ddf0[start:end, 1] <- rep_len(times[i], 4)
  ddf0[start:end, 2] <- rep_len(lines[i], 4)
  ddf0[start:end, 3] <- as.array(names(dexprs0[,i]))
  ddf0[start:end, 4] <- as.numeric(dexprs0[,i])
  
}

cdf <- as.data.frame(cdf)
cdf$Expression <- as.numeric(cdf$Expression)
ddf0 <- as.data.frame(ddf0)
ddf0$Expression <- as.numeric(ddf0$Expression)
ldf$Expression <- exp(ldf$Expression)
ggplot(data = ldf, aes(x = Time, y = Expression)) + 
  geom_violin(aes(fill = Time), alpha = 0.5) +
  facet_grid(rows = vars(Gene), scales = "free_y") +
  geom_dotplot(stackdir = "center", binaxis = "y", #binwidth = 0.03,
               dotsize = 0.8) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  ggtitle("Limma corrected expression - exponential") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  NoLegend() 

ggplot(data = cdf, aes(x = Time, y = Expression)) + 
  geom_violin(aes(fill = Time), alpha = 0.5) +
  facet_grid(rows = vars(Gene), scales = "free_y") +
  geom_dotplot(stackdir = "center", binaxis = "y", #binwidth = 0.03,
               dotsize = 0.8) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  ggtitle("log2 transformed cpm") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  NoLegend() 

ggplot(data = ddf, aes(x = Time, y = Expression)) + 
  geom_violin(aes(fill = Time), alpha = 0.5) +
  facet_grid(rows = vars(Gene), scales = "free_y") +
  geom_dotplot(stackdir = "center", binaxis = "y", #binwidth = 0.03,
               dotsize = 0.8) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  ggtitle("DESeq2 corrected basemean \nusing multiple covariates") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  NoLegend() 

ggplot(data = ddf0, aes(x = Time, y = Expression)) + 
  geom_violin(aes(fill = Time), alpha = 0.5) +
  facet_grid(rows = vars(Gene), scales = "free_y") +
  geom_dotplot(stackdir = "center", binaxis = "y", #binwidth = 0.03,
               dotsize = 0.8) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  ggtitle("DESeq2 corrected basemean \nusing cell line as covariate") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  NoLegend() 
