# Chuxuan Li 02/21/2022
# Redo differential experession analysis on 5-line data using MAST and check GO
#term enrichment in the DEGs

# init ####
library(MAST)
library(Seurat)
library(ggplot2)
library(stringr)
library(future)

plan("multisession", workers = 2)
options(future.globals.maxSize = 207374182400)

setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata")
load("renormalized_4celltype_obj.RData")

# function from 20-line analysis to calculate time-specific DEG ####
find_time_DEG <- function(obj, metadata, metadata_name, test2use, min.pct){
  
  # Generic function to calculate time-specific DEG for an obj within subgroups
  #diveded according to another ident, e.g. cell type. 
  # Inputs:
  # obj: a Seurat object, containing counts and metadata for time and another identity
  # metadata: a named vector, one metadata column from obj
  # metadata_name: a string, the name of the metadata column to be used in FindMarkers()
  # test2use: a string, indicates which test to send to FindMarkers() as input 
  # min.pct: a numerical, minimum percentage of cells expressing the gene, send to FindMarkers()
  # Returns: 
  # df_list_2_return: a list, contains two data frames df_1vs0 and df_6vs0, which
  #contains DEG comparing 1hr vs 0hr and DEG comparing 6hr vs 0hr
  
  Idents(obj) <- metadata_name
  j <- 0
  # calc hour 1 vs 0
  for (i in sort(unique(metadata))) {
    print(i)
    j = j + 1
    if (j == 1) {
      print("first cell type")
      df_1vs0 <- FindMarkers(object = obj,
                             slot = "scale.data",
                             ident.1 = "1hr",
                             group.by = 'time.ident', 
                             subset.ident = i,
                             ident.2 = "0hr",
                             min.pct = min.pct,
                             logfc.threshold = 0.0,
                             test.use = test2use)
      df_1vs0$gene_symbol <- rownames(df_1vs0)
      df_1vs0$cell_type <- i
    } else {
      print("others")
      df_to_append <- FindMarkers(object = obj,
                                  slot = "scale.data",
                                  ident.1 = "1hr",
                                  group.by = 'time.ident',
                                  subset.ident = i,
                                  ident.2 = "0hr",
                                  min.pct = min.pct,
                                  logfc.threshold = 0.0,
                                  test.use = test2use)
      df_to_append$gene_symbol <- rownames(df_to_append)
      df_to_append$cell_type <- i
      df_1vs0 <- rbind(df_1vs0,
                       df_to_append,
                       make.row.names = F,
                       stringsAsFactors = F)
    }
  }
  
  j <- 0
  # calculate 6hr vs 0hr
  for (i in sort(unique(metadata))) {
    print(i)
    j = j + 1
    if (j == 1) {
      print("first cell type")
      df_6vs0 <- FindMarkers(object = obj,
                             slot = "scale.data",
                             ident.1 = "6hr",
                             group.by = 'time.ident',
                             subset.ident = i,
                             ident.2 = "0hr",
                             min.pct = min.pct,
                             logfc.threshold = 0.0,
                             test.use = test2use)
      df_6vs0$gene_symbol <- rownames(df_6vs0)
      df_6vs0$cell_type <- i
    } else {
      print("others")
      df_to_append <- FindMarkers(object = obj,
                                  slot = "scale.data",
                                  ident.1 = "6hr",
                                  group.by = 'time.ident',
                                  subset.ident = i,
                                  ident.2 = "0hr",
                                  min.pct = min.pct,
                                  logfc.threshold = 0.0,
                                  test.use = test2use)
      df_to_append$gene_symbol <- rownames(df_to_append)
      df_to_append$cell_type <- i
      df_6vs0 <- rbind(df_6vs0,
                       df_to_append,
                       make.row.names = F,
                       stringsAsFactors = F)
    }
  }
  df_list_2_return <- vector(mode = "list", length = 2L)
  df_list_2_return[[1]] <- df_1vs0
  df_list_2_return[[2]] <- df_6vs0
  
  return(df_list_2_return)
}

# Test for DE ####
unique(human_only$cell.line.ident)
# sum(RNAseq_integrated_labeled$cell.line.ident == "unmatched") # 25173
# 
# pure_human <- subset(RNAseq_integrated_labeled, cell.line.ident %in% c("CD_54", "CD_27",
#                                                                        "CD_26", "CD_08",
#                                                                        "CD_25"))
# pure_human$spec.cell.type <- "others"
# pure_human$spec.cell.type[pure_human$broad.cell.type == "NEFM_pos_glut"] <- "NEFM_pos_glut"
# pure_human$spec.cell.type[pure_human$broad.cell.type == "NEFM_neg_glut"] <- "NEFM_neg_glut"
# pure_human$spec.cell.type[pure_human$broad.cell.type == "NPC"] <- "NPC"
# pure_human$spec.cell.type[pure_human$broad.cell.type == "GABA"] <- "GABA"
# unique(pure_human$spec.cell.type)

DefaultAssay(filtered_obj) <- "RNA"
unique(filtered_obj$cell.type)
filtered_obj <- ScaleData(filtered_obj)
DEG_5line_0.04 <- find_time_DEG(filtered_obj, 
                                metadata = filtered_obj$cell.type, 
                                metadata_name = "cell.type", 
                                test2use = "MAST", 
                                min.pct = 0.04)
DEG_5line_0.05 <- find_time_DEG(filtered_obj, 
                                metadata = filtered_obj$cell.type, 
                                metadata_name = "cell.type", 
                                test2use = "MAST", 
                                min.pct = 0.05)
DEG_5line_0.06 <- find_time_DEG(filtered_obj, 
                                metadata = filtered_obj$cell.type, 
                                metadata_name = "cell.type", 
                                test2use = "MAST", 
                                min.pct = 0.06)
DEG_5line_0.07 <- find_time_DEG(filtered_obj, 
                                metadata = filtered_obj$cell.type, 
                                metadata_name = "cell.type", 
                                test2use = "MAST", 
                                min.pct = 0.07) 


DEG_5line_0.1 <- find_time_DEG(filtered_obj, 
                                metadata = filtered_obj$cell.type, 
                                metadata_name = "cell.type", 
                                test2use = "MAST", 
                                min.pct = 0.1)
DEG_5line_0.15 <- find_time_DEG(filtered_obj, 
                               metadata = filtered_obj$cell.type, 
                               metadata_name = "cell.type", 
                               test2use = "MAST", 
                               min.pct = 0.15)
save("DEG_5line_0.06", 
     file = "renormalized_DEG_5line_0.06.RData")
save("DEG_5line_0.05", "DEG_5line_0.06", "DEG_5line_0.07",
     file = "renormalized_DEG_5line_0.05-7.RData")
save("DEG_5line_0.08", "DEG_5line_0.09", "DEG_5line_0.1",
     file = "renormalized_DEG_5line_0.08-10.RData")

#save("pure_human", file = "mapped_to_demuxed_barcodes_pure_human.RData")
DefaultAssay(integrated) <- "integrated"
#integrated <- ScaleData(integrated)

for (i in unique(integrated$lib.ident)){
  temp <- subset(integrated, lib.ident == i)
  for (j in unique(temp$cell.line.ident)){
    print(i)
    print(j)
    line <- subset(temp, cell.line.ident == j)
    print(sum(line[rownames(line) == "BDNF", ]@assays$RNA@counts))
  }
}
