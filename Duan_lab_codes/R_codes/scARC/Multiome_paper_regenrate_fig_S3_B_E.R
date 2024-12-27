# Siwei Zhang 23 Aug 2024
# Re-generate RNA umaps pre and post Harmony
# RScript for 20230515_QC_plots_for_018-030_data.R part 2 

# init ####
{
  library(Seurat)
  library(Signac)
  library(SeuratObject)
  
  library(RColorBrewer)
  library(viridis)
  library(ggplot2)
  
  library(future)
  
  library(purrr)
  library(patchwork)
  library(scales)
  
  library(stringr)
  library(readr)
}

{
  plan("multisession", workers = 2)
  options(expressions = 20000)
  options(future.globals.maxSize = 207374182400)
  
}

source("../Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/stacked_vln_copy.R")




setwd("~/NVME/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/")

# load("~/NVME/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/integrated_018-030_RNAseq_obj_test_QC.RData")

system.time(load("018-030_RNA_integrated_labeled_with_harmony.RData"))
# user  system elapsed 
# 793.036  69.511 910.933 

# load("anchors_after_cleaned_lst_QC.RData")
# integrated <- IntegrateData(anchorset = anchors,
#                             verbose = T)
#save(cleaned_lst, file = "final_lst_prepped_4_SCT_integration.RData")
#save(anchors, file = "anchors_4_integration.RData")
# save(integrated, file = "integrated_018-030_RNAseq_obj_test_remove_10lines.RData")
# load("integrated_018-030_RNAseq_obj_test_remove_10lines.RData")

# load("integrated_018-030_RNAseq_obj_test_QC.RData")
# integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, assay = "RNA",
#                                                    pattern = "^MT-")
# integrated$time.ident <-
#   paste0(str_sub(integrated$orig.ident,
#                  start = -1L), "hr")
# 
# DefaultAssay(integrated) <- "RNA"
# 
# 
DefaultAssay(integrated_labeled) <- "integrated" 

Idents(integrated_labeled) <- "orig.ident"
Idents(integrated_labeled) <- "seq.batch.ident"
Idents(integrated_labeled) <- "time.ident"
Idents(integrated_labeled) <- "seurat_clusters"

unique(integrated_labeled$cell.type.forplot)

show_col(hue_pal()(23))

DimPlot(integrated_labeled, # 9x4.5 in or 5x4.5 in
        reduction = "umap", 
        label = T,
        repel = T)# + ggtitle("Test")

DimPlot(integrated_labeled,
        label = F,
        alpha = 0.2,
        group.by = "time.ident")# +
DimPlot(integrated_labeled,
        label = F,
        alpha = 0.2,
        group.by = "seq.batch.ident")# +
DimPlot(integrated_labeled,
        label = F,
        alpha = 0.2,
        group.by = "orig.ident") +
  NoLegend()
  # ggtitle("by time point") +
  # theme(text = element_text(size = 12))

integrated_labeled$cell.type.forplot <-
  factor(integrated_labeled$cell.type.forplot,
         levels = c("GABA",
                    'NEFM+ glut',
                    'NEFM- glut',
                    'glut?',
                    'unidentified'))
DimPlot(integrated_labeled,
        label = F,
        alpha = 0.2,
        cols = c("#B33E52", 
                 "darkorange",
                 
                 "#CCAA7A",
                 "#E6D2B8", 
                 "darkgreen"),
        # cols = brewer.pal(n = 5,
        #                   name = "Dark2"),
        group.by = "cell.type.forplot") +
  NoLegend()# +
# integrated_labeled_21clusters[, integrated_labeled_21clusters$seurat_clusters %in% 21] <- 20

integrated_labeled$seurat_clusters

Idents(integrated_labeled) <- "seurat_clusters"
typeof(Idents(integrated_labeled))


integrated_labeled_21clusters <-
  subset(x = integrated_labeled,
         cells = (integrated_labeled$seurat_clusters %in% c(21, 22)),
         invert = T)

cell_lines_all <-
  unique(integrated_labeled@meta.data$cell.line.ident)

DimPlot(integrated_labeled_21clusters, # 9x4.5 in or 5x4.5 in
        reduction = "umap", 
        label = T,
        repel = T,
        cols = hue_pal()(21),
        alpha = 0.2)# + ggtitle("Test")
# DimPlot(integrated_labeled, 
#         reduction = "harmony")

Idents(integrated_labeled_21clusters)
DefaultAssay(integrated_labeled_21clusters) <- "SCT"
FeaturePlot(object = integrated_labeled_21clusters,
            features = c("POU5F1",
                         "NANOG",
                         "NEFM",
                         "MAP2"),
            ncol = 2)

trimmed_markers <- c("GAD1", "GAD2", "SLC17A6", "SLC17A7",
                     "EBF1", # striatal
                     "SEMA3E", # subcerebral
                     "BCL11B",  # cortical
                     "SST", # inhibitory
                     "SATB2",  "NEFM", # excitatory
                     "VIM", "SOX2",  #NPC
                     "SERTAD4", "FOXG1",  # forebrain
                     "POU3F2", "LHX2", # general cortex
                     "ADCYAP1", "CUX1", "CUX2", 
                     "MAP2", "DCX") # pan-neuron

StackedVlnPlot(obj = integrated_labeled_21clusters, 
               features = trimmed_markers) +
  coord_flip()



image_4_pdf <-
  vector(mode = "list",
         length = length(cell_lines_all))

for (i in 1:length(cell_lines_all)) {
  print(paste("i =",
              i))
  sub_2_plot <-
    subset(integrated_labeled,
           subset = (cell.line.ident == cell_lines_all[i]))
  Idents(sub_2_plot) <- "orig.ident"
  # pdf(file = paste0("Fig_S3E/",
  #                   "Line_",
  #                   cell_lines_all[i],
  #                   ".pdf"),
  #     width = 2.25,
  #     height = 3, 
  #     compress = T,
  #     title = paste0("Line_",
  #                    cell_lines_all[i]))
  image_4_pdf[[i]] <-
    DimPlot(sub_2_plot,
            reduction = "umap", 
            cols = brewer.pal(n = 4,
                              name = "Set1"),
            alpha = 0.2) +
    ggtitle(paste0("Line_",
                   cell_lines_all[i]))
  # print(image_4_pdf)
  # dev.off()
}

saveRDS(image_4_pdf,
        file = "all_100_lines_umap.RData")

integrated_labeled$lib.ident.4.plot <-
  str_c('Library ',
        str_split(integrated_labeled$lib.ident,
                  pattern = '-',
                  simplify = T)[, 1])

integrated_labeled$time.ident.4.plot <-
  str_c(str_split(integrated_labeled$lib.ident,
                  pattern = '-',
                  simplify = T)[, 2], 
        ' hr')

libraries_all <-
  unique(integrated_labeled$lib.ident.4.plot)
image_4_pdf <-
  vector(mode = "list",
         length = length(libraries_all))

for (i in 1:length(libraries_all)) {
  print(paste("i =",
              i))
  sub_2_plot <-
    subset(integrated_labeled,
           subset = (lib.ident.4.plot == libraries_all[i]))
  Idents(sub_2_plot) <- "time.ident.4.plot"
  # pdf(file = paste0("Fig_S3E/",
  #                   "Line_",
  #                   cell_lines_all[i],
  #                   ".pdf"),
  #     width = 2.25,
  #     height = 3, 
  #     compress = T,
  #     title = paste0("Line_",
  #                    cell_lines_all[i]))
  image_4_pdf[[i]] <-
    DimPlot(sub_2_plot,
            reduction = "umap", 
            cols = brewer.pal(n = 4,
                              name = "Set1"),
            alpha = 0.2) +
    ggtitle(libraries_all[i])
  # print(image_4_pdf)
  # dev.off()
}

assembled_plot <-
  wrap_plots(image_4_pdf,
             ncol = 6,
             byrow = T,
             guides = "keep",
             axes = "keep")
  
print(assembled_plot)


#####
dev.off()


DefaultAssay(integrated_labeled) <- "RNA"
DefaultAssay(integrated_labeled) <- "SCT"
FeaturePlot(integrated_labeled,
            features = c("SLC17A6",
                         "SLC17A7",
                         "NEFM",
                         "GAD1",
                         "GAD2",
                         "FOS"),
            ncol = 3)
Idents(integrated_labeled) <- "cell.type.forplot"
unique(integrated_labeled$cell.type.forplot)

integrated_labeled$cell.type.forplot <-
  factor(integrated_labeled$cell.type.forplot,
         levels = c("GABA",
                    'NEFM+ glut',
                    'NEFM- glut',
                    'glut?',
                    'unidentified'))
DimPlot(integrated_labeled,
        alpha = 0.2,
        cols = brewer.pal(n = 5,
                          name = "Dark2"))

Idents(integrated_labeled) <- "lib.ident"
integrated_labeled$lib.time.4.plot <-
  integrated_labeled$lib.ident
unique(integrated_labeled$lib.time.4.plot)

all_lib_times <-
  unique(integrated_labeled$lib.time.4.plot)

for (i in 1:54) {
  print(i)
  integrated_labeled$lib.time.4.plot[integrated_labeled$lib.time.4.plot == all_lib_times[i]] <-
    str_c("CD",
          all_lib_times[i],
          sep = "")
}

integrated_labeled$lib.time.4.plot <-
  str_c(integrated_labeled$lib.time.4.plot,
        "h",
        sep = "")
# integrated_labeled$lib.time.4.plot <-
#   str_c("CD",
#         '-',
#         integrated_labeled$lib.ident,
#         "hr")
# unique(integrated_labeled$orig.ident)
# unique(integrated_labeled$lib.ident.4.plot)

# integrated_labeled$lib.time.4.plot <-
#   factor(integrated_labeled$lib.time.4.plot)
# levels(integrated_labeled$lib.time.4.plot)
# integrated_labeled$lib.time.4.plot <-
#   str_replace_all(string = integrated_labeled$lib.time.4.plot,
#                  pattern = "^CD-20087",
# #                  replacement = "20087")
# unique(integrated_labeled$lib.time.4.plot)
# 
# lib_time_4_plot <-
#   unique(integrated_labeled$lib.time.4.plot)

# for (i in 55:length(lib_time_4_plot)) {
#   integrated_labeled$lib.time.4.plot[integrated_labeled$lib.time.4.plot == lib_time_4_plot[i]] <-
#     str_remove_all(string = lib_time_4_plot[i],
#                    pattern = "^CD-")
# }



Idents(integrated_labeled) <- "lib.time.4.plot"
DotPlot(integrated_labeled,
        features = rev(c("NPAS4", "FOS", 
                         "VGF", "BDNF")),
        # idents = "lib.ident",
        cols = c("lightgrey", "brown2")) +
  theme(axis.text.x = element_text(angle = 315,
                                   vjust = 0,
                                   hjust = 0),
        legend.position = "bottom") +
  coord_flip() 


ncol(integrated_labeled)
unique(integrated_labeled$seurat_clusters)

integrated_labeled_21clusters <-
  subset(x = integrated_labeled,
         cells = (integrated_labeled$seurat_clusters %in% c(21, 22)),
         invert = T)
ncol(integrated_labeled_21clusters)
unique(integrated_labeled_21clusters$seurat_clusters)
# user  system elapsed 
# 21.569   0.866  22.489
system.time({
  sub_2_plot <-
    subset(x = integrated_labeled,
           subset = (cell.line.ident == cell_lines_all[1]))
})
DimPlot(sub_2_plot,
        reduction = "umap", 
        cols = brewer.pal(n = 4,
                          name = "Set1"),
        alpha = 0.2) +
  ggtitle(paste0("Line_",
                 cell_lines_all[1]))

DefaultAssay(integrated_labeled)
Idents(integrated_labeled) <- "orig.ident"
# ! if use subset(cells = ...), invert= must be T (no idea why) #####
system.time({
  sub_2_plot <-
    subset(x = integrated_labeled,
           cells = (integrated_labeled$seurat_clusters %in% c(1, 2)),
           invert = F)
})

sub_2_plot <-
  subset(x = integrated_labeled,
         cells = !(integrated_labeled$orig.ident %in% cell_lines_all[2]),
         invert = T)
ncol(sub_2_plot)


sub_2_plot <-
  subset(x = integrated_labeled,
         cells = (integrated_labeled$seurat_clusters %in% c(21, 22)),
         invert = F)



Idents(sub_2_plot) <- "orig.ident"
# pdf(file = paste0("Fig_S3E/",
#                   "Line_",
#                   cell_lines_all[i],
#                   ".pdf"),
#     width = 2.25,
#     height = 3, 
#     compress = T,
#     title = paste0("Line_",
#                    cell_lines_all[i]))
# image_4_pdf[[i]] <-
DimPlot(sub_2_plot,
        reduction = "umap", 
        cols = brewer.pal(n = 4,
                          name = "Set1"),
        alpha = 0.2) +
  ggtitle(paste0("Line_",
                 cell_lines_all[1]))
# print(image_4_pdf)
# dev.off()


system.time(
  integrated_labeled <-
    readRDS("all_100_lines_umap.RData")
)

unique(integrated_labeled_21clusters$lib.ident)
unique(integrated_labeled_21clusters$cell.line.ident)

df_100_metadata <-
  read.csv("100line_metadata_df.csv")
df_100_metadata <-
  df_100_metadata[df_100_metadata$time_point == "0hr", ]
df_100_aff_case_metadata <-
  df_100_metadata[df_100_metadata$aff == "case", ]

# integrated_labeled$aff <- NA

integrated_labeled_21clusters$aff <- "control"

unique(df_100_aff_case_metadata$cell_line)
integrated_labeled_21clusters$aff[integrated_labeled_21clusters$cell.line.ident %in%
                                    unique(df_100_aff_case_metadata$cell_line)] <-
  "case"

integrated_labeled_21clusters <-
  readRDS("integrated_018_030_21_clusters_w_aff.RData")
# 
# for (i in 1:nrow(df_100_metadata)) {
#   integrated_labeled_21clusters[integrated_labeled_21clusters$cell.line.ident]
# }
# for (i in 1:nrow(df_100_metadata)) {
#   integrated_labeled_21clusters[integrated_labeled_21clusters$cell.line.ident]
# }

Idents(integrated_labeled_21clusters) <-
  "aff"
DefaultAssay(integrated_labeled_21clusters) <- "RNA"

integrated_labeled_21clusters$aff <-
  factor(integrated_labeled_21clusters$aff,
         levels = c("case",
                    "control"))
DimPlot(integrated_labeled_21clusters,
        reduction = "umap",
        cells = WhichCells(integrated_labeled_21clusters, 
                           idents = "case"),
        cols = c("darkred"),
        alpha = 0.5) +
  NoLegend()

DimPlot(integrated_labeled_21clusters,
        reduction = "umap",
        # cells = WhichCells(integrated_labeled_21clusters, 
        #                    idents = "case"),
        cols = brewer.pal(n = 3,
                          name = "Set1"),
        # cols = c("darkred",
        #          "darkblue"),
        alpha = 0.2) 


table(integrated_labeled_21clusters$aff)
DimPlot(integrated_labeled_21clusters)
DimPlot(integrated_labeled)
saveRDS(integrated_labeled_21clusters,
        file = "integrated_018_030_21_clusters_w_aff.RData")


DefaultAssay(integrated_labeled) <- "SCT"
FeaturePlot(integrated_labeled,
            features = "FGF5")
