# Chuxuan Li 08/30/2022
# calculate the correlation between the logFC obtained by DESeq2 and Limma, and
#p-value obtained by these two methods.

# init ####
library(readr)
library(readxl)
library(stringr)
library(ggplot2)

DESeq2_paths <- sort(list.files(path = "./pseudobulk_DE/res_use_combatseq_mat/filtered_by_basemean_only", 
           pattern = "1v0|6v0", full.names = T))
limma_paths <- sort(list.files(path = "./pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm_9in18_res/unfiltered_by_padj",
                          pattern = "1v0|6v0", full.names = T))
DESeq2_dfs <- vector("list", length(DESeq2_paths))
limma_dfs <- vector("list", length(limma_paths))
for (i in 1:length(DESeq2_dfs)) {
  DESeq2_dfs[[i]] <- read_csv(DESeq2_paths[i])
  limma_dfs[[i]] <- read_csv(limma_paths[i])
  names(DESeq2_dfs)[i] <- paste("DESeq2", 
                                str_extract(DESeq2_paths[i], "[A-Za-z]+_[1-6]v0"), 
                                sep = "-")
  names(limma_dfs)[i] <- paste("limma", 
                                str_remove(str_extract(limma_paths[i], "[A-Za-z]+_res_[1-6]v0"), "res_"), 
                                sep = "-")
}

# logFC ####
szgenes <- read_excel("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                      sheet = "prioritized SZ single genes")
szgenes <- sort(szgenes$gene.symbol)
logFC_dfs <- vector("list", length(DESeq2_dfs))
for (i in 1:length(DESeq2_dfs)) {
  #genelist <- sort(intersect(DESeq2_dfs[[i]]$gene, limma_dfs[[i]]$genes))
  des <- DESeq2_dfs[[i]]$log2FoldChange[DESeq2_dfs[[i]]$gene %in% szgenes]
  lim <- limma_dfs[[i]]$logFC[limma_dfs[[i]]$genes %in% szgenes]
  df <- as.data.frame(cbind(des, lim))
  colnames(df) <- c(names(DESeq2_dfs)[i], names(limma_dfs)[i])
  logFC_dfs[[i]] <- df
}
names(logFC_dfs) <- str_remove(names(DESeq2_dfs), "DESeq2-")
for (i in 1:length(logFC_dfs)) {
  jpeg(filename = paste0(names(logFC_dfs)[i], "_logFC_scatter.jpeg"))
  p <- ggplot(logFC_dfs[[i]], aes(x = `DESeq2-npglut_6v0`, y = `limma-npglut_6v0`)) + 
    geom_point()
}
ggplot(logFC_dfs[[6]], aes(x = `DESeq2-npglut_6v0`, y = `limma-npglut_6v0`)) + 
  geom_point()

# pvalue ####
pval_dfs <- vector("list", length(DESeq2_dfs))
for (i in 1:length(DESeq2_dfs)) {
  #genelist <- sort(intersect(DESeq2_dfs[[i]]$gene, limma_dfs[[i]]$genes))
  des <- DESeq2_dfs[[i]]$pvalue[DESeq2_dfs[[i]]$gene %in% szgenes]
  lim <- limma_dfs[[i]]$P.Value[limma_dfs[[i]]$genes %in% szgenes]
  df <- as.data.frame(cbind(des, lim))
  colnames(df) <- c(names(DESeq2_dfs)[i], names(limma_dfs)[i])
  pval_dfs[[i]] <- df
}
names(pval_dfs) <- str_remove(names(DESeq2_dfs), "DESeq2-")
for (i in 1:length(logFC_dfs)) {
  jpeg(filename = paste0(names(logFC_dfs)[i], "_logFC_scatter.jpeg"))
  p <- ggplot(logFC_dfs[[i]], aes(x = `DESeq2-npglut_6v0`, y = `limma-npglut_6v0`)) + 
    geom_point()
}
ggplot(pval_dfs[[6]], aes(x = `DESeq2-npglut_6v0`, y = `limma-npglut_6v0`)) + 
  geom_point()
