# Siwei 25 Oct 2023
# make some plots for ASHG 2023

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(future)
library(parallel)
library(grDevices)


plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)
options(future.seed = T)

load("018-029_RNA_integrated_labeled.RData")
Assays(integrated_labeled)
DefaultAssay(integrated_labeled)
Idents(integrated_labeled) <- "cell.type"
# integrated_labeled$

DefaultAssay(integrated_labeled) <-
  "integrated"

Idents(integrated_labeled) <- "fine.cell.type"

integrated_labeled$cell.type.new <-
  integrated_labeled$cell.type.forplot
integrated_labeled$cell.type.new[integrated_labeled$cell.type.forplot %in% "unknown"] <-
  "Not used for analysis"
integrated_labeled$cell.type.new[integrated_labeled$cell.type.forplot %in% '?glut'] <-
  "Not used for analysis"



DimPlot(integrated_labeled, 
        # reduction = "umap", 
        label = TRUE,
        repel = T,
        # pt.size = 0.3,
        # # alpha = 1,
        cols = brewer.pal(n = 11,
                          name = "Paired")
        ) +
  # geom_point(alpha = 1) +
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 10))

# need to increase resolution
system.time({
  integrated_labeled <- 
    FindClusters(integrated_labeled, 
                 resolution = 5, 
                 random.seed = 42)
})
DimPlot(integrated_labeled, 
        group.by = "seurat_clusters",
        # reduction = "umap", 
        label = TRUE,
        repel = T,
        # pt.size = 0.3,
        # # alpha = 1,
        # cols = brewer.pal(n = 6,
        #                   name = "Dark2")
) +
  # geom_point(alpha = 1) +
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 10))

Idents(integrated_labeled) <- 
  "seurat_clusters"


cell_type_to_export <-
  vector(mode = "list",
         length = 3L)
names(cell_type_to_export) <-
  c("GABA",
    "nmglut",
    "npglut")

cell_type_to_export[["GABA"]] <-
  subset(x = integrated_labeled,
         idents = 0, 
         downsample = 2000)
cell_type_to_export[["nmglut"]] <-
  subset(x = integrated_labeled,
         idents = 8, 
         downsample = 2000)
cell_type_to_export[["npglut"]] <-
  subset(x = integrated_labeled,
         idents = 7, 
         downsample = 2000)

DefaultAssay(cell_type_to_export[[1]]) <- "RNA"
DefaultAssay(cell_type_to_export[[2]]) <- "RNA"
DefaultAssay(cell_type_to_export[[3]]) <- "RNA"

cell_type_to_export[[1]]@assays$integrated <- NULL
cell_type_to_export[[2]]@assays$integrated <- NULL
cell_type_to_export[[3]]@assays$integrated <- NULL

cell_type_to_export[[1]]@assays$SCT <- NULL
cell_type_to_export[[2]]@assays$SCT <- NULL
cell_type_to_export[[3]]@assays$SCT <- NULL

saveRDS(object = cell_type_to_export,
        file = "cell_3_types_from_76_lines_2K_each_4_Velmeshev.RDs")

{
  integrated_labeled$cell.type.new <-
    "Not used in analysis"
  integrated_labeled$cell.type.new[integrated_labeled$seurat_clusters %in% 
                                     c(76, 93, 60, 70, 69, 
                                       71, 81, 91)] <-
    "NPC"
  # integrated_labeled$cell.type.new[integrated_labeled$seurat_clusters %in% 
  #                                    c(65, 94, 88, 80,
  #                                      74, 87, 39, 84,
  #                                      48, 75, 36, 55, 42, 64, 51,
  #                                      57, 66, 52, 36,
  #                                      64, 45, 27, 51)] <-
  #   "Not used in analysis"
  integrated_labeled$cell.type.new[integrated_labeled$seurat_clusters %in% 
                                     c(34, 43, 44, 35, 58,
                                       25, 89, 1, 11, 23,
                                       28, 29, 2, 79,
                                       16, 38, 42,
                                       0, 14, 73, 3,
                                       15,
                                       22, 30, 26, 5, 9)] <-
    "GABA"
  integrated_labeled$cell.type.new[integrated_labeled$seurat_clusters %in% 
                                     c(63, 
                                       67, 
                                       31, 
                                       50,
                                       18,
                                       # 62, 
                                       10,
                                       # 59, 
                                       25,
                                       40,
                                       82,
                                       47,
                                       # 41, 
                                       24, 
                                       19,
                                       17,
                                       46,
                                       # 56, 
                                       13
                                       , 21
                                       )] <-
    'NEFM+ Glut'
  integrated_labeled$cell.type.new[integrated_labeled$seurat_clusters %in% 
                                     c(7, 8, 6, 33, 85, 
                                       54, 32, 61, 92, 37,
                                       20, 12, 4, 78)] <-
    'NEFM- Glut'
  # [1] "27" "36" "39" "41" "45" "48" "49" "51" "55" "56" "57" "59" "62" "64" "65" "66"
  # [17] "68" "75" "77" "80" "83" "84" "86" "87" "88" "90" "94" "95"
  integrated_labeled$cell.type.new[integrated_labeled$seurat_clusters %in% 
                                     c(82, 
                                       # 17, 
                                       53, 54, 52, 72, 74, 
                                       45)] <-
    ' '
  
  integrated_labeled$cell.type.new <-
    factor(integrated_labeled$cell.type.new,
           levels = c('NEFM- Glut',
                      'NEFM+ Glut',
                      "GABA",
                      "Not used in analysis",
                      "NPC",
                      ' '))
  
  DimPlot(integrated_labeled, 
          group.by = "cell.type.new",
          # reduction = "umap", 
          label = TRUE,
          repel = T,        
          # pt.size = ifelse(test = integrated_labeled$cell.type.new %in% 
          #                 "not_to_plot",
          #               yes = 0,
          #               no = 1),
          cols = c("orange", #nmglut
                   "#CCAA7AFF", #npglut
                   "#B33E52FF", #GABA
                   # "#f29116", #?glut
                   "darkseagreen", #NPC
                   "#4293DBFF",#unknown
                   "#FFFFFF00"#not plotted
          )
          # pt.size = 0.3,
          # # alpha = 1,
          # cols = brewer.pal(n = 6,
          #                   name = "Dark2")
  ) +
    # geom_point(alpha = 1) +
    ggtitle("labeled by cell type") +
    theme(text = element_text(size = 12), 
          axis.text = element_text(size = 10))
  
}

sort(unique(as.character(integrated_labeled$seurat_clusters[integrated_labeled$cell.type.new %in%
                                                              "Not used in analysis"])))




Idents(integrated_labeled) <- "fine.cell.type"
DimPlot(integrated_labeled, 
        reduction = "umap", 
        label = TRUE, 
        repel = T,
        pt.size = 0.3, 
        cols = c("#E6D2B8", #nmglut
                 "#CCAA7A", #npglut
                 "#B33E52", #GABA
                 "#f29116", #?glut
                 "#347545", #NPC
                 "#4293db"#unknown
        )) +
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 10))


Assays(multiomic_obj)
DefaultAssay(integrated_labeled)
Idents(integrated_labeled)

DefaultAssay(integrated_labeled) <-
  "integrated"

DimPlot(integrated_labeled, 
        reduction = "umap", 
        label = TRUE, 
        repel = T,
        pt.size = 0.3, 
        cols = c("#E6D2B8", #nmglut
                 "#CCAA7A", #npglut
                 "#B33E52", #GABA
                 "#f29116", #?glut
                 "#347545", #NPC
                 "#4293db", #not in use
                 "#00000000"#not plotted
        )) +
  ggtitle("labeled by cell type") +
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 10))
