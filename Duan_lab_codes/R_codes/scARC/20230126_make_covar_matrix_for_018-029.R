# Chuxuan Li 01/25/2023
# Generate covariate table for 018-029 data for case v control DE analysis

# init ####
library(Seurat)
library(readr)
library(stringr)
library(readxl)

MGS_table <- read_excel("~/NVME/scARC_Duan_018/MGS_iPSC lines_60samples_scRNA_ATAC-seq_bulkATAC-seq status_AK_HZ.xlsx")
CIRM_table <- read_excel("~/NVME/scARC_Duan_018/CIRM control iPSC lines_40_Duan (003).xlsx")

# clean and make metadata table ####
MGS_table <- MGS_table[, c(9, 1, 5, 4, 3)]
colnames(MGS_table) <- c("group", "cell_line", "age", "sex", "aff")
MGS_table$sex[MGS_table$sex == 2] <- "F"
MGS_table$sex[MGS_table$sex == 1] <- "M"
MGS_table$cell_line <- str_replace(MGS_table$cell_line, "CD00000", "CD_")
MGS_table <- MGS_table[!is.na(MGS_table$group), ]
MGS_table$group[str_detect(MGS_table$group, "& [0-9]+$")] <- 
  str_replace(str_extract(MGS_table$group[str_detect(MGS_table$group, "& [0-9]+$")], "& [0-9]+$"), "& ", "")

CIRM_table <- CIRM_table[, c(4, 6, 13, 16, 21)]
CIRM_table <- CIRM_table[!is.na(CIRM_table$`Catalog ID`), ]
colnames(CIRM_table) <- c("group", "cell_line", "age", "sex", "aff")
CIRM_table$sex[CIRM_table$sex == "Female"] <- "F"
CIRM_table$sex[CIRM_table$sex == "Male"] <- "M"
CIRM_table$aff <- "control"

covar_table <- rbind(MGS_table, CIRM_table)
covar_table$time <- "0hr"
covar_table_full <- rbind(covar_table, covar_table, covar_table)
covar_table_full$time[111:220] <- "1hr"
covar_table_full$time[221:330] <- "6hr"


# add cellular composition data ####
load("~/NVME/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/025_covariates_full_df.RData")
load("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/029_RNA_integrated_labeled.RData")

integrated_labeled$linextime.ident <- NA
times <- sort(unique(integrated_labeled$time.ident))
lines <- sort(unique(integrated_labeled$cell.line.ident))
for (t in times) {
  for (l in lines) {
    integrated_labeled$linextime.ident[integrated_labeled$time.ident == t &
                                         integrated_labeled$cell.line.ident == l] <-
      paste(l, t, sep = "_")
  }
}
lts <- unique(integrated_labeled$linextime.ident)
med_df <- data.frame(cell_line = rep_len("", length(lts)),
                     time = rep_len("", length(lts)),
                     GABA_counts = rep_len("", length(lts)),
                     GABA_fraction = rep_len("", length(lts)),
                     nmglut_counts = rep_len("", length(lts)),
                     nmglut_fraction = rep_len("", length(lts)),
                     npglut_counts = rep_len("", length(lts)),
                     npglut_fraction = rep_len("", length(lts)),
                     total_counts = rep_len("", length(lts)))
unique(integrated_labeled$cell.type)
integrated_labeled <- subset(integrated_labeled, cell.type != "unidentified")
for (i in 1:length(lts)) {
  obj <- subset(integrated_labeled, linextime.ident == lts[[i]])
  line <- unique(obj$cell.line.ident)
  time <- unique(obj$time.ident)
  types <- unique(obj$cell.type)
  print(line)
  print(time)
  print(types)
  med_df$cell_line[i] <- line
  med_df$time[i] <- time
  med_df$GABA_counts[i] <- sum(obj$cell.type == "GABA")
  med_df$nmglut_counts[i] <- sum(obj$cell.type == "nmglut")
  med_df$npglut_counts[i] <- sum(obj$cell.type == "npglut")
  med_df$total_counts[i] <- ncol(obj)
  med_df$GABA_fraction[i] <- sum(obj$cell.type == "GABA") / ncol(obj)
  med_df$nmglut_fraction[i] <- sum(obj$cell.type == "nmglut") / ncol(obj)
  med_df$npglut_fraction[i] <- sum(obj$cell.type == "npglut") / ncol(obj)
}
med_df$seq_batch <- "029"

# combine into covariate table
covar_add <- med_df
covar_add$aff <- rep_len("", nrow(covar_add))
covar_add$age <- rep_len("", nrow(covar_add))
covar_add$sex <- rep_len("", nrow(covar_add))
covar_add$batch <- rep_len("", nrow(covar_add))
for (i in 1:length(covar_add$cell_line)) {
  line <- covar_add$cell_line[i]
  print(line)
  covar_add$aff[i] <- covar_table$aff[covar_table$cell_line == line]
  covar_add$age[i] <- covar_table$age[covar_table$cell_line == line]
  covar_add$sex[i] <- covar_table$sex[covar_table$cell_line == line]
  covar_add$batch[i] <- covar_table$group[covar_table$cell_line == line]
}
covar_df$seq_batch <- "025"
covar_table_final <- rbind(covar_df, covar_add[, match(colnames(covar_df), colnames(covar_add))])

load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/covariates_pooled_together_for_all_linextime_5+18+20.RData")
colnames(covar_table) <- str_replace(colnames(covar_table), "disease", "aff")
colnames(covar_table) <- str_replace(colnames(covar_table), "group", "batch")
# remove cell lines in covar_table (contains 022/bad 20 line data) that are in covar_df (025 redid)
covar_table_51820_clean <- covar_table[!covar_table$cell_line %in% covar_df$cell_line, ]
redone_lines <- unique(covar_table$cell_line[covar_table$cell_line %in% covar_df$cell_line])
save(redone_lines, file = "redone_lines_char_vector.RData")

covar_table_51820_clean$seq_batch <- "022"
covar_table_51820_clean$seq_batch[covar_table_51820_clean$cell_line %in% c("CD_27", "CD_54", "CD_08", "CD_25", "CD_26")] <- "018"
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000_with_covars_from_30Sep2022_single_cell_DE.RData")
obj_024 <- obj
rm(obj)
lines_024 <- unique(obj_024$cell.line.ident)
covar_table_51820_clean$seq_batch[covar_table_51820_clean$cell_line %in% lines_024] <- "024"

covar_table_final <- rbind(covar_table_51820_clean, covar_table_final[, match(colnames(covar_table_51820_clean), colnames(covar_table_final))])

covar_table_final <- covar_table_final[order(covar_table_final$cell_line), ]
covar_table_final <- covar_table_final[order(covar_table_final$time), ]
covar_table_final$GABA_fraction <- as.numeric(covar_table_final$GABA_fraction)
covar_table_final$nmglut_fraction <- as.numeric(covar_table_final$nmglut_fraction)
covar_table_final$npglut_fraction <- as.numeric(covar_table_final$npglut_fraction)
covar_table_final$age <- as.numeric(covar_table_final$age)
save(covar_table_final, file = "018-029_covar_table_final.RData")

sum(covar_table_final$aff == "case" & covar_table_final$time == "0hr") #26
sum(covar_table_final$aff == "control" & covar_table_final$time == "0hr") #50
