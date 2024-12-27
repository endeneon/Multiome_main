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
  library(stringr)
  library(readr)
}


{
  plan("multisession", workers = 6)
  options(expressions = 20000)
  options(future.globals.maxSize = 207374182400)
  
}
# load("~/NVME/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/integrated_018-030_RNAseq_obj_test_QC.RData")

# system.time(load("018-030_RNA_integrated_labeled_with_harmony.RData"))
# user  system elapsed 
# 793.036  69.511 910.933 
# setwd("~/NVME/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/")
load("anchors_after_cleaned_lst_QC.RData")
integrated <- IntegrateData(anchorset = anchors,
                            verbose = T)
#save(cleaned_lst, file = "final_lst_prepped_4_SCT_integration.RData")
#save(anchors, file = "anchors_4_integration.RData")
save(integrated, file = "integrated_018-030_RNAseq_obj_test_remove_10lines_no_harmony.RData")


load("integrated_018-030_RNAseq_obj_test_remove_10lines.RData")

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
# DefaultAssay(integrated_labeled) <- "integrated"
# 
# DimPlot(integrated_labeled, 
#         reduction = "umap")
# DimPlot(integrated_labeled, 
#         reduction = "harmony")
