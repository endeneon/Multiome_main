# Chuxuan Li 05/16/2023
# Generate metadata including sex, age, cellular composition, aff for integrated
# data frame for Lifan

# init ####
library(stringr)

load("./018-030_RNA_integrated_labeled_with_harmony.RData")

# clean original metadata files ####
MGS_df <- readxl::read_xlsx("../MGS_iPSC lines_60samples_scRNA_ATAC-seq_bulkATAC-seq status_AK_HZ.xlsx")
MGS_df <- MGS_df[, c(1, 3:5)]
colnames(MGS_df) <- c("cell_line", "aff", "sex", "age")
MGS_df$cell_line <- str_replace(MGS_df$cell_line, "CD00000", "CD_")
MGS_df$sex <- str_replace(MGS_df$sex, "1", "M")
MGS_df$sex <- str_replace(MGS_df$sex, "2", "F")

CIRM_df <- readxl::read_xlsx("../CIRM control iPSC lines_40_Duan (003).xlsx")
CIRM_df <- CIRM_df[, c(6, 21, 16, 13)]
colnames(CIRM_df) <- c("cell_line", "aff", "sex", "age")
CIRM_df <- CIRM_df[!is.na(CIRM_df$cell_line), ]
CIRM_df$sex <- str_replace(CIRM_df$sex, "Male", "M")
CIRM_df$sex <- str_replace(CIRM_df$sex, "Female", "F")
CIRM_df$aff <- str_replace(CIRM_df$aff, "No", "control")

# clean and add to combined metadata ####
meta_df <- rbind(MGS_df, CIRM_df)

used_lines <- unique(integrated_labeled$cell.line.ident)
meta_df <- meta_df[meta_df$cell_line %in% used_lines, ]
meta_df$aff[is.na(meta_df$aff)] <- "control"

meta_df$coculture_batch <- ""
for (l in used_lines) {
  print(l)
  lib <- unique(integrated_labeled$lib.ident[integrated_labeled$cell.line.ident == l])
  print(lib)
  meta_df$coculture_batch[meta_df$cell_line == l] <- str_remove(lib, "-[0|1|6]$")
}

meta_df$seq_batch <- ""
for (l in used_lines) {
  print(l)
  sbatch <- unique(integrated_labeled$seq.batch.ident[integrated_labeled$cell.line.ident == l])
  print(sbatch)
  meta_df$seq_batch[meta_df$cell_line == l] <- sbatch
}

# add time point information #####
meta_df_0hr <- meta_df
meta_df_1hr <- meta_df
meta_df_6hr <- meta_df
meta_df_0hr$time_point <- "0hr"
meta_df_1hr$time_point <- "1hr"
meta_df_6hr$time_point <- "6hr"
meta_df <- rbind(meta_df_0hr, meta_df_1hr, meta_df_6hr)

# add cellular composition ####
meta_df$GABA_count <- 0
meta_df$nmglut_count <- 0
meta_df$npglut_count <- 0
meta_df$intermediateGlut_count <- 0
meta_df$unidentified_count <- 0

meta_df$GABA_porportion <- 0
meta_df$nmglut_porportion <- 0
meta_df$npglut_porportion <- 0
meta_df$intermediateGlut_porportion <- 0
meta_df$unidentified_porportion <- 0

meta_df$total_number_of_cells <- 0

lines <- sort(unique(meta_df$cell_line))
times <- sort(unique(meta_df$time_point))

for (l in lines) {
  for (t in times) {
    selectInd <- meta_df$cell_line == l & meta_df$time_point == t
    meta_df$GABA_count[selectInd] <- 
      sum(integrated_labeled$cell.line.ident == l & integrated_labeled$time.ident == t & 
            integrated_labeled$cell.type == "GABA")
    meta_df$nmglut_count[selectInd] <- 
      sum(integrated_labeled$cell.line.ident == l & integrated_labeled$time.ident == t & 
            integrated_labeled$cell.type == "nmglut")
    meta_df$npglut_count[selectInd] <- 
      sum(integrated_labeled$cell.line.ident == l & integrated_labeled$time.ident == t & 
            integrated_labeled$cell.type == "npglut")
    meta_df$intermediateGlut_count[selectInd] <- 
      sum(integrated_labeled$cell.line.ident == l & integrated_labeled$time.ident == t & 
            integrated_labeled$cell.type == "glut?")
    meta_df$unidentified_count[selectInd] <- 
      sum(integrated_labeled$cell.line.ident == l & integrated_labeled$time.ident == t & 
            integrated_labeled$cell.type == "unidentified")
    
    meta_df$total_number_of_cells[selectInd] <- 
      sum(integrated_labeled$cell.line.ident == l & integrated_labeled$time.ident == t)
    
    meta_df$GABA_porportion[selectInd] <- 
      meta_df$GABA_count[selectInd] / meta_df$total_number_of_cells[selectInd]
    meta_df$nmglut_porportion[selectInd] <- 
      meta_df$nmglut_count[selectInd] / meta_df$total_number_of_cells[selectInd]
    meta_df$npglut_porportion[selectInd] <- 
      meta_df$npglut_count[selectInd] / meta_df$total_number_of_cells[selectInd]
    meta_df$intermediateGlut_porportion[selectInd] <- 
      meta_df$intermediateGlut_count[selectInd] / meta_df$total_number_of_cells[selectInd]
    meta_df$unidentified_porportion[selectInd] <- 
      meta_df$unidentified_count[selectInd] / meta_df$total_number_of_cells[selectInd]
  }
}

write.table(meta_df, file = "100line_metadata_df.csv", quote = F, sep = ",", row.names = F)
