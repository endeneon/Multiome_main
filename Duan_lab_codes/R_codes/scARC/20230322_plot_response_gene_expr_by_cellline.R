# Chuxuan Li 02/01/2023
# Look at expression of response genes for 029 RNA data

# init ####
library(Seurat)
library(ggplot2)
library(RColorBrewer)

times <- unique(integrated_labeled$time.ident)

# feature plots ####
for (t in times) {
  tobj <- subset(integrated_labeled, time.ident == t)
  png(paste0("./Dimplots/response_gene_expression_featplot_", t, ".png"), 
      width = 900, height = 900)
  p <- FeaturePlot(tobj, c("FOS", "NPAS4", "BDNF", "VGF"), ncol = 2)
  print(p)
  dev.off()
  pdf(paste0("./Dimplots/response_gene_expression_featplot_", t, ".pdf"), 
      width = 900, height = 900)
  p <- FeaturePlot(tobj, c("FOS", "NPAS4", "BDNF", "VGF"), ncol = 2)
  print(p)
  dev.off()
}

# dotplot ####
lines <- sort(unique(integrated_labeled$cell.line.ident))
integrated_labeled$linextime.ident <- ""
for (l in lines) {
  for (t in times) {
    integrated_labeled$linextime.ident[integrated_labeled$cell.line.ident == l &
                                         integrated_labeled$time.ident == t] <- 
      paste(l, t, sep = " ")
  }
}
unique(integrated_labeled$linextime.ident)
DotPlot(integrated_labeled, assay = "SCT", cols = c("white", "red3"),
        features = c("EGR1", "NR4A1", "NPAS4", "FOS"), group.by = "linextime.ident") +
  coord_flip() +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10))

DotPlot(integrated_labeled, assay = "SCT", cols = c("white", "red3"),
        features = c("IGF1", "VGF", "BDNF"), group.by = "linextime.ident") +
  coord_flip() +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10))

save(integrated_labeled, file = "029_RNA_integrated_labeled.RData")

# use the following RData for 74 line plotting
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/018-029_RNA_integrated_labeled_with_harmony.RData")

table(integrated_labeled$fine.cell.type)
table(integrated_labeled$cell.type.forplot)
table(integrated_labeled$fine.cell.type)
DefaultAssay(integrated_labeled) <- "integrated"

DimPlot(integrated_labeled,
        reduction = "umap",
        group.by  = "cell.type")
