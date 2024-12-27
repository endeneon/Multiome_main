# Siwei 29 Jun 2023
# Extract cell barcodes by line:time:type (nmglut/npglut/gaba)

# init ####
library(Seurat)
library(ggplot2)
library(pals)
library(stringr)
library(future)

plan("multisession", workers = 4)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

# load data from Lexi's RData #####
load("Analysis_part2_mapped_to_Hg38_only/integrated_labeled.Rdata")

Idents(integrated_labeled)
unique(integrated_labeled$cell.line.ident)
# [1] "CW20112" "CW20079" "CW20137" "CD_38"   "CW20063" "CW20111" "CW50099" "CW70030"
# [9] "CW50148" "CW50059" "CW50021" "CW50150" "CW50094" "CW70282" "CW50023" "CW70344"
# [17] "CW70305" "CD_39"   "CW70079" "CW50037"
unique(integrated_labeled$time.ident) # "0hr" "1hr" "6hr"
unique(integrated_labeled$cell.type) # "npglut"       "GABA"         "unidentified" "nmglut"

## extract CW20063 group
for (time_ident in unique(integrated_labeled$time.ident)) {
  for (cell_line_ident in c("CD_38", "CW20063", "CW20079", "CW20112", "CW20137")) {
    for (cell_type in c("npglut", "GABA", "nmglut")) {
      output_file_name <-
        str_c("Duan030",
              time_ident,
              cell_line_ident,
              cell_type,
              "barcodes.txt",
              sep = "_")
      print(output_file_name)

      current_barcode_list <-
        colnames(integrated_labeled)[(integrated_labeled$time.ident %in% time_ident) &
                                       (integrated_labeled$cell.line.ident %in% cell_line_ident) &
                                       (integrated_labeled$cell.type %in% cell_type)]
      # remove the "_1" suffix
      # print(current_barcode_list)
      try({ # line CW30154 has very low GABA content, vector of zero length may be produced
        current_barcode_list <-
          str_split(string = current_barcode_list,
                    pattern = "_",
                    simplify = T)[, 1]

        if (!dir.exists("barcodes_by_time_line_type_CW20063")) {
          dir.create("barcodes_by_time_line_type_CW20063")
        }

        # write out the barcode list
        write.table(current_barcode_list,
                    file = str_c("barcodes_by_time_line_type_CW20063",
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


## extract CW20111 group
for (time_ident in unique(integrated_labeled$time.ident)) {
  for (cell_line_ident in c("CW20111", "CW50059", "CW50099", "CW50148", "CW70030")) {
    for (cell_type in c("npglut", "GABA", "nmglut")) {
      output_file_name <-
        str_c("Duan030",
              time_ident,
              cell_line_ident,
              cell_type,
              "barcodes.txt",
              sep = "_")
      print(output_file_name)

      current_barcode_list <-
        colnames(integrated_labeled)[(integrated_labeled$time.ident %in% time_ident) &
                                       (integrated_labeled$cell.line.ident %in% cell_line_ident) &
                                       (integrated_labeled$cell.type %in% cell_type)]
      # remove the "_1" suffix
      # print(current_barcode_list)
      try({ # line CW30154 has very low GABA content, vector of zero length may be produced
        current_barcode_list <-
          str_split(string = current_barcode_list,
                    pattern = "_",
                    simplify = T)[, 1]

        if (!dir.exists("barcodes_by_time_line_type_CW20111")) {
          dir.create("barcodes_by_time_line_type_CW20111")
        }

        # write out the barcode list
        write.table(current_barcode_list,
                    file = str_c("barcodes_by_time_line_type_CW20111",
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

## extract CW50094 group
for (time_ident in unique(integrated_labeled$time.ident)) {
  for (cell_line_ident in c("CW50021", "CW50023", "CW50094", "CW50150", "CW70282")) {
    for (cell_type in c("npglut", "GABA", "nmglut")) {
      output_file_name <-
        str_c("Duan030",
              time_ident,
              cell_line_ident,
              cell_type,
              "barcodes.txt",
              sep = "_")
      print(output_file_name)

      current_barcode_list <-
        colnames(integrated_labeled)[(integrated_labeled$time.ident %in% time_ident) &
                                       (integrated_labeled$cell.line.ident %in% cell_line_ident) &
                                       (integrated_labeled$cell.type %in% cell_type)]
      # remove the "_1" suffix
      # print(current_barcode_list)
      try({ # line CW30154 has very low GABA content, vector of zero length may be produced
        current_barcode_list <-
          str_split(string = current_barcode_list,
                    pattern = "_",
                    simplify = T)[, 1]

        if (!dir.exists("barcodes_by_time_line_type_CW50094")) {
          dir.create("barcodes_by_time_line_type_CW50094")
        }

        # write out the barcode list
        write.table(current_barcode_list,
                    file = str_c("barcodes_by_time_line_type_CW50094",
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

## extract CW70305 group
for (time_ident in unique(integrated_labeled$time.ident)) {
  for (cell_line_ident in c("CD_39", "CW50037", "CW70079", "CW70305", "CW70344")) {
    for (cell_type in c("npglut", "GABA", "nmglut")) {
      output_file_name <-
        str_c("Duan030",
              time_ident,
              cell_line_ident,
              cell_type,
              "barcodes.txt",
              sep = "_")
      print(output_file_name)

      current_barcode_list <-
        colnames(integrated_labeled)[(integrated_labeled$time.ident %in% time_ident) &
                                       (integrated_labeled$cell.line.ident %in% cell_line_ident) &
                                       (integrated_labeled$cell.type %in% cell_type)]
      # remove the "_1" suffix
      # print(current_barcode_list)
      try({ # line CW30154 has very low GABA content, vector of zero length may be produced
        current_barcode_list <-
          str_split(string = current_barcode_list,
                    pattern = "_",
                    simplify = T)[, 1]

        if (!dir.exists("barcodes_by_time_line_type_CW70305")) {
          dir.create("barcodes_by_time_line_type_CW70305")
        }

        # write out the barcode list
        write.table(current_barcode_list,
                    file = str_c("barcodes_by_time_line_type_CW70305",
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
