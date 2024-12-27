# run under Rscript
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
  plan("multisession", workers = 1)
  options(expressions = 20000)
  options(future.globals.maxSize = 207374182400)
  
}

setwd("~/NVME/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/")

# load("~/NVME/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/integrated_018-030_RNAseq_obj_test_QC.RData")

system.time(load("018-030_RNA_integrated_labeled_with_harmony.RData"))

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

saveRDS(integrated_labeled,
        file = "all_integrated_data_4_dotplot_27Aug2024.RDs")