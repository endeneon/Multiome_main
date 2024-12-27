# Chuxuan Li 03/30/2023
# Check NPC percentage for control cell lines for Hanwen

# init ####
library(Seurat)
library(readr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(readxl)

load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/integrated_labeled.RData")
obj_025 <- integrated_labeled
rm(integrated_labeled)
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000_with_covars_from_30Sep2022_single_cell_DE.RData")
obj_024 <- obj
rm(obj)
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/renormalized_4celltype_obj.RData")
obj_018 <- filtered_obj
rm(filtered_obj)
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/20line_codes_and_rdata/demux_20line_integrated_labeled_obj.RData")
obj_022 <- integrated_labeled
rm(integrated_labeled)

data_excel <- read_excel("iPSC lines_selection-MINND.xlsx", sheet = "EA")

# find cell lines in objects ####
lines <- str_replace(data_excel$`new iPSC line ID`, "CD00000", "CD_")
obj_lst <- list(obj_018, obj_022, obj_024, obj_025)
types_lst <- vector("list", length(lines))
subobj_lst <- vector("list", length(lines))
for (j in 1:length(lines)) {
  l <- lines[j]
  for (i in 1:length(obj_lst)) {
    obj <- obj_lst[[i]]
    obj_lines <- unique(obj$cell.line.ident)
    if (l %in% obj_lines) {
      subobj <- subset(obj, cell.line.ident == l)
      subobj_lst[[j]] <- subobj
      types <- unique(subobj$cell.type)
      types_lst[[j]] <- types
    }
  }
}
names(types_lst) <- lines

# calculate NPC percentage ####
data_excel$cleaned_cell_line_ID <- str_replace(data_excel$`new iPSC line ID`, "CD00000", "CD_")

for (i in 1:length(types_lst)) {
  l <- names(types_lst)[i]
  types <- types_lst[[i]]
  print(types)
  if ("NPC" %in% types) {
    data_excel$`neuron diff`[data_excel$cleaned_cell_line_ID == l] <- sum(subobj_lst[[i]]$cell.type == "NPC") / length(subobj_lst[[i]]$cell.type)
  } else {
    data_excel$`neuron diff`[data_excel$cleaned_cell_line_ID == l] <- "NA"   
  }
}

write_csv(data_excel, file = "iPSC lines_selection_with_NPC_proportion.csv", quote = "none")
