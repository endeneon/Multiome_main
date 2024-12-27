# Chuxuan Li 07/07/2022
# plot heatmap for gene expression and gene activity score for 18line with
#genes separated into categories according to trend in 1v0 6v1 changes

# init ####
library(readr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

# load data ####
npglut_1v0 <- read_csv("./pseudobulk_DE/res_use_combatseq_mat/filtered_by_basemean_only/npglut_1v0_basemean_filtered_only.csv")
npglut_6v1 <- read_csv("./pseudobulk_DE/res_use_combatseq_mat/filtered_by_basemean_only/npglut_6v1_basemean_filtered_only.csv")
npglut_6v0 <- read_csv("./pseudobulk_DE/res_use_combatseq_mat/filtered_by_basemean_only/npglut_6v0_basemean_filtered_only.csv")
#load("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/multiomic_obj.RData")

# define categories ####

# lists of genes from one of the time periods
ux <- npglut_1v0[npglut_1v0$log2FoldChange > 0 & npglut_1v0$pvalue < 0.1, ]
xu <- npglut_6v0[npglut_6v0$log2FoldChange > 0 & npglut_6v0$pvalue < 0.1, ]
dx <- npglut_1v0[npglut_1v0$log2FoldChange < 0 & npglut_1v0$pvalue < 0.1, ]
xd <- npglut_6v0[npglut_6v0$log2FoldChange < 0 & npglut_6v0$pvalue < 0.1, ]
fx <- npglut_1v0[npglut_1v0$pvalue > 0.1, ]
xf <- npglut_6v0[npglut_6v0$pvalue > 0.1, ]

# find intersections between the above lists to find the 8 specific categories
uu_gene <- intersect(ux$gene, xu$gene) # up up
ud_gene <- intersect(ux$gene, xd$gene) # up down
uf_gene <- intersect(ux$gene, xf$gene) # up flat (not significant in 6v1)
du_gene <- intersect(dx$gene, xu$gene) # down up
dd_gene <- intersect(dx$gene, xd$gene) # down down
df_gene <- intersect(dx$gene, xf$gene) # down flat
fu_gene <- intersect(fx$gene, xu$gene) # flat up
fd_gene <- intersect(fx$gene, xd$gene) # flat down
#ff_gene <- intersect(fx$gene, xf$gene)
total_gene <- intersect(npglut_1v0$gene, npglut_6v1$gene)
write.table(x = uu_gene, 
            file = "./pseudobulk_DE/res_use_combatseq_mat/0-1-6_dynamic_analysis/uugene.txt",
            quote = F, sep = "\t", row.names = F, col.names = F)
write.table(x = total_gene, 
            file = "./pseudobulk_DE/res_use_combatseq_mat/0-1-6_dynamic_analysis/bgdgene.txt",
            quote = F, sep = "\t", row.names = F, col.names = F)

# make data frame with 0hr at baseline, fold changes at 1hr, 6hr ####
genelist <- c(uu_gene, ud_gene, uf_gene,
              du_gene, dd_gene, df_gene,
              fu_gene, fd_gene)
genecategories <- data.frame(gene = genelist, 
                             category = c(rep_len("uu", length(uu_gene)),
                                          rep_len("ud", length(ud_gene)),
                                          rep_len("uf", length(uf_gene)),
                                          rep_len("du", length(du_gene)),
                                          rep_len("dd", length(dd_gene)),
                                          rep_len("df", length(df_gene)),
                                          rep_len("fu", length(fu_gene)),
                                          rep_len("fd", length(fd_gene))))
logFCdf_0 <- data.frame(gene = genelist, exp = as.numeric(rep_len(0, length(genelist))),
                        time = rep_len("0hr", length(genelist)), category = rep_len("ff", length(genelist)))
logFCdf_1 <- data.frame(gene = genelist, exp = as.numeric(rep_len(0, length(genelist))),
                        time = rep_len("1hr", length(genelist)), category = rep_len("ff", length(genelist)))
logFCdf_6 <- data.frame(gene = genelist, exp = as.numeric(rep_len(0, length(genelist))), 
                        time = rep_len("6hr", length(genelist)), category = rep_len("ff", length(genelist)))

fillLogFCdf <- function(gene, resdf) {
  fc <- resdf$log2FoldChange[resdf$gene == gene]
  return(fc)
}
fillCategory <- function(gene) {
  cat <- genecategories$category[genecategories$gene == gene]
  return(cat)
}
logFCdf_1$exp <- as.numeric(lapply(X = genelist, FUN = fillLogFCdf, resdf = npglut_1v0))
logFCdf_6$exp <- as.numeric(lapply(X = genelist, FUN = fillLogFCdf, resdf = npglut_6v0))
logFCdf_0$category <- as.character(lapply(X = genelist, FUN = fillCategory))
logFCdf_1$category <- as.character(lapply(X = genelist, FUN = fillCategory))
logFCdf_6$category <- as.character(lapply(X = genelist, FUN = fillCategory))

logFCdf <- rbind(logFCdf_0, logFCdf_1, logFCdf_6)

# plot line plots with log2FoldChange ####
ggplot(logFCdf, aes(x = time, y = exp, group = gene, color = category)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = brewer.pal(6, "Set3")) +
  theme_minimal()

ggplot(logFCdf, aes(x = time, y = exp, group = gene, color = category)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("#8DD3C7", "grey", "grey", "grey", "grey", "grey")) +
  theme_minimal()

# compute average gex and gact matrices ####
unique(multiomic_obj$cell.type)
nmglut <- subset(multiomic_obj, cell.type == "NEFM_neg_glut")
npglut <- subset(multiomic_obj, cell.type == "NEFM_pos_glut")
GABA <- subset(multiomic_obj, cell.type == "GABA")


times <- sort(unique(npglut$time.ident))
gex_df <- data.frame(npglut_zero = genelist, npglut_one = genelist, npglut_six = genelist)
gact_df <- data.frame(npglut_zero = genelist, npglut_one = genelist, npglut_six = genelist)

makeGexDf <- function(object, genelist) {
  for (j in 1:length(times)) {
    obj <- subset(object, time.ident == times[j])
    for (i in 1:length(genelist)) {
      g <- genelist[i]
      if (i == 1) {
        DefaultAssay(obj) <- "SCT"
        gex_col <- sum(exp(obj@assays$SCT@data[rownames(obj) == g, ])) /
          ncol(obj)
      } else {
        DefaultAssay(obj) <- "SCT"
        col_to_append <- sum(exp(obj@assays$SCT@data[rownames(obj) == g, ])) /
          ncol(obj)
        gex_col <- rbind(gex_col, col_to_append)
      }
    }
    if (j == 1) {
      gex_mat <- gex_col
    } else {
      gex_mat <- rbind(gex_mat, gex_col)
    }
    rownames(gex_mat) <- genelist
  }
  return(gex_mat)
}
npglut_gex_df <- makeGexDf(GABA, genelist)
gex_clust <- hclust(dist(gex_df, method = "maximum"))
heatmap.2(x = as.matrix(gex_df),
          scale = "row", 
          trace = "none", 
          Rowv = gex_clust$order, 
          Colv = T, 
          dendrogram = "none",
          main = "gene expression")

