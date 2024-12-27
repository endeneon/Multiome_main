# Chuxuan Li 01/05/2023
# plot SZ gene fold change and p-value from MAST-direct results with 5pt filter

# init ####
library(ggplot2)
library(ggrepel)
library(GGally)
library(GSEABase)
library(limma)
library(reshape2)
library(data.table)
library(knitr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(NMF)
library(rsvd)
library(RColorBrewer)
library(MAST)
library(Seurat)
library(readr)

respath <- list.files("./MAST_direct_DE/filter_by_5pct_in_either_time", pattern = "^result", full.names = T)
reslist <- vector("list", length(respath))
for (i in 1:length(respath)) {
  load(respath[i])
  reslist[[i]] <- fcHurdle
}
names(reslist) <- str_extract(respath, "[1|6]v0_[A-Za-z]+")
SZ_genes <- c("NAB2", "GRIN2A", "IMMP2L", "SNAP91", 
              "PTPRK", "SP4", "RORB", "TCF4", "DIP2A")

# plot ####
for (i in 1:length(reslist)) {
  if (i == 1) {
    df_to_plot <- reslist[[i]][reslist[[i]]$primerid %in% SZ_genes, ]
    df_to_plot$time <- str_extract(names(reslist)[i], "[1|6]v0")
    df_to_plot$cell.type <- str_remove(str_extract(names(reslist)[i], "_[A-Za-z]+"), "_")
    df_to_plot$neglogp <- (-1) * log10(df_to_plot$`Pr(>Chisq)`)
  } else {
    df_to_append <- reslist[[i]][reslist[[i]]$primerid %in% SZ_genes, ]
    df_to_append$time <- str_extract(names(reslist)[i], "[1|6]v0")
    df_to_append$cell.type <- str_remove(str_extract(names(reslist)[i], "_[A-Za-z]+"), "_")
    df_to_append$neglogp <- (-1) * log10(df_to_append$`Pr(>Chisq)`)
    df_to_plot <- rbind(df_to_plot, df_to_append)
  }
}
df_to_plot$neglogp[df_to_plot$fdr > 0.05] <- NA 
df_to_plot$coef[df_to_plot$fdr > 0.05] <- NA 
df_to_plot$neglogp[is.infinite(df_to_plot$neglogp)] <- 500

ggplot(df_to_plot,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = neglogp,
           fill = coef)) +
  geom_point(shape = 21) +
  ylab("Cell Type") +
  xlab("") +
  scale_fill_gradientn(colours = c("royalblue3", "white", "darkred")) +
  facet_grid(cols = vars(primerid)) + 
  theme_bw() +
  labs(size = "-log10(p-value)", fill = "coefficient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 

# compared to limma+combat+9 out of 18 sample cpm>1 ####
pathlist <- list.files("./pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm_9in18_res/unfiltered_by_padj", 
                       full.names = T,
                       pattern = "1v0|6v0")
pathlist
reslist <- vector('list', length(pathlist))
for (i in 1:length(reslist)){
  reslist[[i]] <- read_csv(pathlist[i])
}
namelist <- str_extract(pathlist, "[A-Za-z]+_res_[1|6]v0")
namelist
names(reslist) <- namelist
for (i in 1:length(reslist)) {
  if (i == 1) {
    df_to_plot <- reslist[[i]][reslist[[i]]$genes %in% SZ_genes, ]
    df_to_plot$time <- str_extract(names(reslist)[i], "[1|6]v0")
    df_to_plot$cell.type <- str_remove(str_extract(names(reslist)[i], "^[A-Za-z]+_res"), "_res")
    df_to_plot$neglogp <- (-1) * log10(df_to_plot$P.Value)
  } else {
    df_to_append <- reslist[[i]][reslist[[i]]$genes %in% SZ_genes, ]
    df_to_append$time <- str_extract(names(reslist)[i], "[1|6]v0")
    df_to_append$cell.type <- str_remove(str_extract(names(reslist)[i], "^[A-Za-z]+_res"), "_res")
    df_to_append$neglogp <- (-1) * log10(df_to_append$P.Value)
    df_to_plot <- rbind(df_to_plot, df_to_append)
  }
}
df_to_plot$neglogp[df_to_plot$adj.P.Val > 0.05] <- NA 
df_to_plot$logFC[df_to_plot$adj.P.Val > 0.05] <- NA 
df_to_plot$neglogp[is.infinite(df_to_plot$neglogp)] <- 500

ggplot(df_to_plot,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = neglogp,
           fill = logFC)) +
  geom_point(shape = 21) +
  ylab("Cell Type") +
  xlab("") +
  scale_fill_gradientn(colours = c("royalblue3", "white", "darkred")) +
  facet_grid(cols = vars(genes)) + 
  theme_bw() +
  labs(size = "-log10(p-value)", fill = "logFC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 

# plot early and late genes too ####
early_late_genes <- c("FOS", "FOSB", "NPAS4", "NR4A1", "BDNF", "VGF", "IGF1")
for (i in 1:length(reslist)) {
  if (i == 1) {
    df_to_plot <- reslist[[i]][reslist[[i]]$primerid %in% early_late_genes, ]
    df_to_plot$time <- str_extract(names(reslist)[i], "[1|6]v0")
    df_to_plot$cell.type <- str_remove(str_extract(names(reslist)[i], "_[A-Za-z]+"), "_")
    df_to_plot$neglogp <- (-1) * log10(df_to_plot$`Pr(>Chisq)`)
  } else {
    df_to_append <- reslist[[i]][reslist[[i]]$primerid %in% early_late_genes, ]
    df_to_append$time <- str_extract(names(reslist)[i], "[1|6]v0")
    df_to_append$cell.type <- str_remove(str_extract(names(reslist)[i], "_[A-Za-z]+"), "_")
    df_to_append$neglogp <- (-1) * log10(df_to_append$`Pr(>Chisq)`)
    df_to_plot <- rbind(df_to_plot, df_to_append)
  }
}
df_to_plot$neglogp[df_to_plot$fdr > 0.05] <- NA 
df_to_plot$coef[df_to_plot$fdr > 0.05] <- NA 
df_to_plot$neglogp[is.infinite(df_to_plot$neglogp)] <- 500

ggplot(df_to_plot,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = neglogp,
           fill = coef)) +
  geom_point(shape = 21) +
  ylab("Cell Type") +
  xlab("") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  facet_grid(cols = vars(factor(primerid, 
                                levels = c("FOS", "FOSB", "NPAS4", "NR4A1", "BDNF", "VGF", "IGF1")))) + 
  theme_bw() +
  labs(size = "-log10(p-value)", fill = "coefficient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 

# plot early and late again for limma ####
for (i in 1:length(reslist)) {
  if (i == 1) {
    df_to_plot <- reslist[[i]][reslist[[i]]$genes %in% early_late_genes, ]
    df_to_plot$time <- str_extract(names(reslist)[i], "[1|6]v0")
    df_to_plot$cell.type <- str_remove(str_extract(names(reslist)[i], "^[A-Za-z]+_res"), "_res")
    df_to_plot$neglogp <- (-1) * log10(df_to_plot$P.Value)
  } else {
    df_to_append <- reslist[[i]][reslist[[i]]$genes %in% early_late_genes, ]
    df_to_append$time <- str_extract(names(reslist)[i], "[1|6]v0")
    df_to_append$cell.type <- str_remove(str_extract(names(reslist)[i], "^[A-Za-z]+_res"), "_res")
    df_to_append$neglogp <- (-1) * log10(df_to_append$P.Value)
    df_to_plot <- rbind(df_to_plot, df_to_append)
  }
}
df_to_plot$neglogp[df_to_plot$adj.P.Val > 0.05] <- NA 
df_to_plot$logFC[df_to_plot$adj.P.Val > 0.05] <- NA 
df_to_plot$neglogp[is.infinite(df_to_plot$neglogp)] <- 500

ggplot(df_to_plot,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = neglogp,
           fill = logFC)) +
  geom_point(shape = 21) +
  ylab("Cell Type") +
  xlab("") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  facet_grid(cols = vars(factor(genes, 
                           levels = c("FOS", "FOSB", "NPAS4", "NR4A1", "BDNF", "VGF", "IGF1")))) + 
  theme_bw() +
  labs(size = "-log10(p-value)", fill = "logFC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 

