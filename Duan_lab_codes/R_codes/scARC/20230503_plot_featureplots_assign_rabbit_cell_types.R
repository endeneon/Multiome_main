# Chuxuan Li 2023/05/03
# Plot cell type specific markers for rabbit cells and assign cell types 

# init ####
library(Seurat)
library(ggplot2)
library(pals)

load("./oc2_mm10_analysis/integrated_rabbit_only_object_filtered2.RData")

# plot feature plots ####
DefaultAssay(integrated_filtered) <- "SCT"
p <- FeaturePlot(integrated_filtered, features = c("ENSOCUG00000012644", 
                                          "ENSOCUG00000017480",
                                          "ENSOCUG00000016389"#,
                                          #"ENSOCUG00000012799"
                                          )) 
p$patches$plots[[1]] <- p$patches$plots[[1]] + ggtitle("GAD1")
p$patches$plots[[2]] <- p$patches$plots[[2]] + ggtitle("GAD2")
p$labels$title <- "SLC17A6"            
p

# NPC, neuron
genes <- c("VIM", "MAP2", "PAX6") #no SOX2 or SOX3
ensids <- c("ENSOCUG00000009222", "ENSOCUG00000005163", "ENSOCUG00000011075")
p <- FeaturePlot(integrated_filtered, features = ensids)
p$patches$plots[[1]] <- p$patches$plots[[1]] + ggtitle("VIM")
p$patches$plots[[2]] <- p$patches$plots[[2]] + ggtitle("MAP2")
p$labels$title <- "PAX6"
print(p)

# astrocyte
genes <- c("GFAP", "AQP4", "SLC1A2", "S100B")
ensids <- c("ENSOCUG00000006895", "ENSOCUG00000008833",
            "ENSOCUG00000000718", "ENSOCUG00000001473")
p <- FeaturePlot(integrated_filtered, features = ensids)
p$patches$plots[[1]] <- p$patches$plots[[1]] + ggtitle("GFAP")
p$patches$plots[[2]] <- p$patches$plots[[2]] + ggtitle("AQP4")
p$patches$plots[[3]] <- p$patches$plots[[3]] + ggtitle("SLC1A2")
p$labels$title <- "S100B"
print(p)


# endothelial
p <- FeaturePlot(integrated_filtered, features = "ENSOCUG00000002835")
p$labels$title <- "CCN2"
p

# monocytes
p <- FeaturePlot(integrated_filtered, features = c("ENSOCUG00000004218",
                                          "ENSOCUG00000015493",
                                          #"ENSOCUG00000006382",
                                          "ENSOCUG00000034782"))
p$patches$plots[[1]] <- p$patches$plots[[1]] + ggtitle("CD14")
p$patches$plots[[2]] <- p$patches$plots[[2]] + ggtitle("ANXA1")
#p$patches$plots[[3]] <- p$patches$plots[[3]] + ggtitle("IFI44")
p$labels$title <- "OAS1"
p

# oligodendrocytes
p <- FeaturePlot(integrated_filtered, features = c(#"MOBP",
                                          "ENSOCUG00000027389",
                                          "ENSOCUG00000002880"
                                          ))
p$patches$plots[[1]] <- p$patches$plots[[1]] + ggtitle("OLIG1")
p$labels$title <- "OPALIN"
p

# microglia
p <- FeaturePlot(integrated_filtered, features = c("ENSOCUG00000021056"))
p$labels$title <- "CSF1R"
p

# plot stacked vlnplot ####
p <- StackedVlnPlot(integrated_filtered, features = c("ENSOCUG00000012644", 
                                                 "ENSOCUG00000017480",
                                                 "ENSOCUG00000016389",
                                                 "ENSOCUG00000009222", "ENSOCUG00000005163", "ENSOCUG00000011075",
                                                 "ENSOCUG00000006895", "ENSOCUG00000008833",
                                                 "ENSOCUG00000000718", "ENSOCUG00000001473",
                                                 "ENSOCUG00000002835",
                                                 "ENSOCUG00000004218",
                                                 "ENSOCUG00000015493",
                                                 #"ENSOCUG00000006382",
                                                 "ENSOCUG00000034782",
                                                 "ENSOCUG00000027389",
                                                 "ENSOCUG00000002880",
                                                 "ENSOCUG00000021056")) +
  coord_flip()
length(p$patches$plots)
genelist <- c("GAD1", "GAD2", "SLC17A6", 
              "VIM", "MAP2", "PAX6", 
              "GFAP", "AQP4", "SLC1A2", "S100B", 
              "CCN2", 
              "CD14", "ANXA1", "OAS1", 
              "OLIG1", "OPALIN", 
              "CSF1R")
for (i in 1:length(p$patches$plots)) {
  p$patches$plots[[i]] <- p$patches$plots[[i]] + ylab(genelist[i])
}
p

VlnPlot(integrated_filtered, features = "ENSOCUG00000012644", assay = "integrated") + ggtitle("GAD1") 
VlnPlot(integrated_filtered, features = "ENSOCUG00000017480", assay = "integrated") + ggtitle("GAD2") 
VlnPlot(integrated_filtered, features = "ENSOCUG00000016389", assay = "integrated") + ggtitle("SLC17A6")

# assign cell types ####
new.cluster.ids <-
  c("GABA", "GABA", "glut", "GABA", "glut",
    "astrocyte", "microglia", "unknown neuron", "oligodendrocyte", "monocyte",
    "unknown neuron", "unknown neuron", "immature neuron", "unknown")
unique(new.cluster.ids)
length(new.cluster.ids)
length(unique(integrated_filtered$seurat_clusters))
names(new.cluster.ids) <- levels(integrated_filtered)
integrated_labeled <- RenameIdents(integrated_filtered, new.cluster.ids)
integrated_labeled$cell.type <- Idents(integrated_labeled)

DimPlot(integrated_labeled, reduction = "umap", label = TRUE, repel = T,
        pt.size = 0.3, cols = rev(as.character(kelly(12)))) +
  ggtitle("Clustering of rabbit cells (integrated all libraries)\nLabeled by cell type") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

integrated_labeled$cell.type.counts <- as.character(integrated_labeled$cell.type)
types <- unique(integrated_labeled$cell.type.counts)
for (i in 1:length(types)) {
  count <- sum(integrated_labeled$cell.type == types[i])
  print(count)
  integrated_labeled$cell.type.counts[integrated_labeled$cell.type.counts == types[i]] <-
    paste0(types[i], "\n", count)
}
unique(integrated_labeled$cell.type.counts)
DimPlot(integrated_labeled, reduction = "umap", group.by = "cell.type.counts", 
        label = TRUE, repel = T, pt.size = 0.3, 
        cols = c("#008856", #astro
        "#0067A5", #GABA
        "#E68FAC", #glut
        "#F38400", #immature neu 
        "#848482", #microglia
        "#A1CAF1", #monocyte
        "#BE0032", #oligo
        "#875692", #unknown 
        "#C2B280" #unknown neu
)) +
  ggtitle("Clustering of rabbit cells (integrated all libraries)\nLabeled by cell type") +
  theme(text = element_text(size = 10), axis.text = element_text(size = 10)) +
  NoLegend()
 

save(integrated_labeled, file = "./oc2_mm10_analysis/integrated_labeled_after_filtering2.RData")
