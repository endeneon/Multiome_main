# Chuxuan Li 03/09/2023
# Generate covariate table for 018-029 data for SZ PRS linear regression 

# init ####
library(Seurat)
library(readr)
library(stringr)
library(readxl)

load("../case_control_DE_analysis/018-029_covar_table_final.RData")
MGS_w_PRS_clozapine <- read_excel("~/NVME/scARC_Duan_018/SZ_PRS_DE_analysis/MGS_SZ case control iPSC lines for neuron stimulation_w_PRS_clozapine.xlsx")
MGS_w_PRS_clozapine$`new iPSC line ID` <- str_replace(MGS_w_PRS_clozapine$`new iPSC line ID`,
                                                      "CD00000", "CD_")
#MGS_w_PRS_clozapine$SZ_PRS[is.na(MGS_w_PRS_clozapine$SZ_PRS)] <- "NA"

covar_table_final$PRS <- NA
for (i in 1:length(covar_table_final$cell_line)) {
  line <- covar_table_final$cell_line[i]
  print(line)
  if (line %in% MGS_w_PRS_clozapine$`new iPSC line ID`) {
    covar_table_final$PRS[i] <-
      MGS_w_PRS_clozapine$SZ_PRS[MGS_w_PRS_clozapine$`new iPSC line ID` == line]
  }
}

covar_table_PRS <- covar_table_final[!is.na(covar_table_final$PRS), ]
sum(covar_table_PRS$aff == "case" & covar_table_PRS$time == "0hr") #24
sum(covar_table_PRS$aff == "control" & covar_table_PRS$time == "0hr") #24

covar_table_PRS$time <- as.factor(covar_table_PRS$time)
covar_table_PRS$aff <- as.factor(covar_table_PRS$aff)
covar_table_PRS$sex <- as.factor(covar_table_PRS$sex)
covar_table_PRS$batch <- as.factor(covar_table_PRS$batch)

save(covar_table_PRS, file = "covar_table_only_lines_w_PRS.RData")
write.table(covar_table_PRS, file = "covar_table_only_lines_w_PRS.csv", 
            quote = F, sep = ",", row.names = F, col.names = T)

# subset count matrix ####
load("../case_control_DE_analysis/018-029_by_type_linextime_combat_adj_mat_lst.RData")

for (i in 1:length(adj_mat_lst)) {
  mat <- adj_mat_lst[[i]]
  adj_mat_lst[[i]] <- mat[, str_extract(colnames(mat), "CD_[0-9][0-9]") %in% covar_table_PRS$cell_line]
}

save(adj_mat_lst, file = "sep_by_type_col_by_linextime_adj_mat_lst_w_only_lines_w_PRS.RData")


load("/nvmefs/scARC_Duan_018/case_control_DE_analysis/018-029_by_type_and_time_combat_adj_mat_lst.RData")
for (i in 1:length(adj_mat_lst)) {
  mat <- adj_mat_lst[[i]]
  adj_mat_lst[[i]] <- mat[, str_extract(colnames(mat), "CD_[0-9][0-9]") %in% covar_table_PRS$cell_line]
}

save(adj_mat_lst, file = "sep_by_type_time_col_by_line_adj_mat_lst_w_only_lines_w_PRS.RData")
