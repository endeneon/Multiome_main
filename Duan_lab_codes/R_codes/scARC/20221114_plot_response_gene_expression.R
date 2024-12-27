# Chuxuan Li 11/14/2022
# plot response gene expression

library(ggplot2)
library(Seurat)
library(cowplot)

load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/integrated_labeled.RData")

RNAobj <- subset(integrated_labeled, cell.type != "unidentified")
RNAobj$typextime <- ""
types <- unique(RNAobj$cell.type)
times <- unique(RNAobj$time.ident)
for (i in 1:length(types)) {
  for (j in 1:length(times)) {
    RNAobj$typextime[RNAobj$cell.type == types[i] &
                       RNAobj$time.ident == times[j]] <- paste(types[i], times[j], sep = " ")
  }
}
unique(RNAobj$typextime)
resp_genes <- c("FOS", "NPAS4", "BDNF", "VGF")
names(resp_genes) <- c("darkred", "darkred", "deepskyblue3", "deepskyblue3")

# feature plot ####
for (i in 1:length(resp_genes)) {
  g <- resp_genes[i]
  c <- names(resp_genes)[i]
  cat(g, c)
  p1 <- FeaturePlot(GABA, features = g, cols = c("grey", c), split.by = "time.ident") +
    theme(legend.position = "right")
  p2 <- FeaturePlot(nmglut, features = g, cols = c("grey", c), split.by = "time.ident") +
    theme(legend.position = "right")
  p3 <- FeaturePlot(npglut, features = g, cols = c("grey", c), split.by = "time.ident") +
    theme(legend.position = "right")
  png(filename = paste0("./Analysis_part2_GRCh38/featplot_for_gex/", g, "_split_by_time.png"), width = 1000, height = 1000)
  p <- plot_grid(p1, p2, p3, align = "v", nrow = 3, label_size = 10, labels = c("GABA", "NEFM- glut", "NEFM+ glut"))
  print(p)
  dev.off()
}


