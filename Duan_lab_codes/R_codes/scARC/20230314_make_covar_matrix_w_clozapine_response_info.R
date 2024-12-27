# Chuxuan Li 03/14/2023
# Generate covariate table for 018-029 data for SZ clozapine linear regression 

# init ####
library(Seurat)
library(readr)
library(stringr)
library(readxl)

load("../case_control_DE_analysis/018-029_covar_table_final.RData")
MGS_w_PRS_clozapine <- read_excel("~/NVME/scARC_Duan_018/SZ_PRS_DE_analysis/MGS_SZ case control iPSC lines for neuron stimulation_w_PRS_clozapine.xlsx")
MGS_w_PRS_clozapine$`new iPSC line ID` <- str_replace(MGS_w_PRS_clozapine$`new iPSC line ID`,
                                                      "CD00000", "CD_")
covar_table_final$clozapine <- NA
for (i in 1:length(covar_table_final$cell_line)) {
  line <- covar_table_final$cell_line[i]
  print(line)
  if (line %in% MGS_w_PRS_clozapine$`new iPSC line ID`) {
    covar_table_final$clozapine[i] <-
      MGS_w_PRS_clozapine$`Taking clozapine when interviewed (DIGS for MGS)`[MGS_w_PRS_clozapine$`new iPSC line ID` == line]
  }
}

covar_table_clozapine <- covar_table_final[covar_table_final$clozapine %in% c("yes", "no"), ]
sum(covar_table_clozapine$aff == "case" & covar_table_clozapine$time == "0hr") #26
sum(covar_table_clozapine$aff == "control" & covar_table_clozapine$time == "0hr") #0

covar_table_clozapine$time <- as.factor(covar_table_clozapine$time)
covar_table_clozapine$aff <- as.factor(covar_table_clozapine$aff)
covar_table_clozapine$sex <- as.factor(covar_table_clozapine$sex)
covar_table_clozapine$batch <- as.factor(covar_table_clozapine$batch)

save(covar_table_clozapine, file = "covar_table_only_lines_w_clozapine_response.RData")
write.table(covar_table_clozapine, file = "covar_table_only_lines_w_clozapine.csv", 
            quote = F, sep = ",", row.names = F, col.names = T)

# subset count matrix ####
load("../case_control_DE_analysis/018-029_by_type_linextime_combat_adj_mat_lst.RData")

for (i in 1:length(adj_mat_lst)) {
  mat <- adj_mat_lst[[i]]
  adj_mat_lst[[i]] <- mat[, str_extract(colnames(mat), "CD_[0-9][0-9]") %in% covar_table_clozapine$cell_line]
}

save(adj_mat_lst, file = "sep_by_type_col_by_linextime_adj_mat_lst_w_only_lines_w_clozapine.RData")


load("/nvmefs/scARC_Duan_018/case_control_DE_analysis/018-029_by_type_and_time_combat_adj_mat_lst.RData")
for (i in 1:length(adj_mat_lst)) {
  mat <- adj_mat_lst[[i]]
  adj_mat_lst[[i]] <- mat[, str_extract(colnames(mat), "CD_[0-9][0-9]") %in% covar_table_clozapine$cell_line]
}

save(adj_mat_lst, file = "sep_by_type_time_col_by_line_adj_mat_lst_w_only_lines_w_clozapine.RData")
