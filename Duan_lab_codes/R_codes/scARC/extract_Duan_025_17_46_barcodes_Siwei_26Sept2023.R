# Siwei Zhang 26 Sept 2023
# Extract barcodes of Duan-025 17/46
# make a subset object

# Get all libraries from 018 to 030, normalize them separately, then integrate

# init ####
{
  
  library(Seurat)
  library(EnsDb.Hsapiens.v86)
  library(patchwork)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(RColorBrewer)
  library(dplyr)
  library(viridis)
  library(graphics)
  library(readxl)
  library(readr)
  
  library(future)
}

plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 161061273600)

load("~/NVME/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/018-030_RNA_integrated_labeled_with_harmony.RData")

unique(integrated_labeled$lib.ident)

integrate_17_46_subset <-
  subset(integrated_labeled,
         subset = (lib.ident %in% c("17-0",
                                    "17-1",
                                    "17-6",
                                    "46-0",
                                    "46-1",
                                    "46-6")))

unique(integrate_17_46_subset$cell.type)
unique(integrate_17_46_subset$time.ident)
unique(integrate_17_46_subset$cell.line.ident)

# extract 17
for (time_ident in unique(integrate_17_46_subset$time.ident)) {
  for (cell_line_ident in c("CD_17",
                            "CD_18",
                            "CD_45",
                            "CD_47")) {
    for (cell_type in c("npglut", 
                        "GABA", 
                        "nmglut")) {
      output_file_name <-
        str_c(time_ident,
              cell_line_ident,
              cell_type,
              "barcodes.txt",
              sep = "_")
      print(output_file_name)
      
      current_barcode_list <-
        colnames(integrate_17_46_subset)[(integrate_17_46_subset$time.ident %in% time_ident) &
                                       (integrate_17_46_subset$cell.line.ident %in% cell_line_ident) &
                                       (integrate_17_46_subset$cell.type %in% cell_type)]
      # remove the "_1" suffix
      # print(current_barcode_list)
      try({ # in case some combinations have very low cell content, 
            # vector of zero length may be produced
        current_barcode_list <-
          str_split(string = current_barcode_list,
                    pattern = "_",
                    simplify = T)[, 1]
        
        # write out the barcode list
        if (!dir.exists("barcodes_by_time_line_type_lib17")) {
          dir.create("barcodes_by_time_line_type_lib17")
        }
        write.table(current_barcode_list,
                    file = str_c("barcodes_by_time_line_type_lib17",
                                 output_file_name,
                                 sep = "/"),
                    quote = F,
                    sep = "\t",
                    row.names = F,
                    col.names = F)
      })
    }
  } 
}

# extract 46
for (time_ident in unique(integrate_17_46_subset$time.ident)) {
  for (cell_line_ident in c("CD_44",
                            "CD_46",
                            "CD_58",
                            "CD_61")) {
    for (cell_type in c("npglut", 
                        "GABA", 
                        "nmglut")) {
      output_file_name <-
        str_c(time_ident,
              cell_line_ident,
              cell_type,
              "barcodes.txt",
              sep = "_")
      print(output_file_name)
      
      current_barcode_list <-
        colnames(integrate_17_46_subset)[(integrate_17_46_subset$time.ident %in% time_ident) &
                                           (integrate_17_46_subset$cell.line.ident %in% cell_line_ident) &
                                           (integrate_17_46_subset$cell.type %in% cell_type)]
      # remove the "_1" suffix
      # print(current_barcode_list)
      try({ # in case some combinations have very low cell content, 
        # vector of zero length may be produced
        current_barcode_list <-
          str_split(string = current_barcode_list,
                    pattern = "_",
                    simplify = T)[, 1]
        
        # write out the barcode list
        if (!dir.exists("barcodes_by_time_line_type_lib46")) {
          dir.create("barcodes_by_time_line_type_lib46")
        }
        write.table(current_barcode_list,
                    file = str_c("barcodes_by_time_line_type_lib46",
                                 output_file_name,
                                 sep = "/"),
                    quote = F,
                    sep = "\t",
                    row.names = F,
                    col.names = F)
      })
    }
  } 
}
