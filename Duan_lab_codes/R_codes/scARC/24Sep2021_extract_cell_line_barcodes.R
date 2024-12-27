# Chuxuan Li Sept 24, 2021
# with merged group 2 and group 8 data already assigned cell types in 14Sep2021_g2_g8.R,
# extract barcodes from the Seurat object, used for subsetting BAM files

# init
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)
library(readr)
library(ggplot2)
library(future)
library(stringr)
# set threads and parallelization
plan("multisession", workers = 6)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

combined <- read_rds("//nvmefs/scARC_Duan_018/R_scARC_Duan_018/DESeq2/combined_seurat_object.rds")

DimPlot(combined, label = T) + NoLegend()
## separate cell types
Glut_barcodes <- 
  combined@assays$SCT@counts@Dimnames[[2]][combined@meta.data$seurat_clusters %in% 
                                             c("2", "3", "5", "19", 
                                               "16", "17", "8", "7",
                                               "9", "12", "13", "18", "22")]
GABA_barcodes <- 
  combined@assays$SCT@counts@Dimnames[[2]][combined@meta.data$seurat_clusters %in% 
                                             c("1", "11", "6", "4")]
Glut_barcodes <- str_sub(Glut_barcodes, end = -5)
GABA_barcodes <- str_sub(GABA_barcodes, end = -5)


## import cell ident barcodes
g_2_0_CD_27_barcodes <- 
  read_csv("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_0_common_CD27_CD54.best.CD_27", 
           col_names = FALSE)
g_2_0_CD_27_barcodes <- unlist(g_2_0_CD_27_barcodes)


CD_27_Glut_0hr_barcodes <- g_2_0_CD_27_barcodes[g_2_0_CD_27_barcodes %in% Glut_barcodes]
write.table(CD_27_Glut_0hr_barcodes,
            file = "line_type_barcodes/CD_27_Glut_0hr_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)
CD_27_Glut_1hr_barcodes <- g_2_1_CD_27_barcodes[g_2_1_CD_27_barcodes %in% Glut_barcodes]
write.table(CD_27_Glut_1hr_barcodes,
            file = "line_type_barcodes/CD_27_Glut_1hr_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)
CD_27_Glut_6hr_barcodes <- g_2_6_CD_27_barcodes[g_2_6_CD_27_barcodes %in% Glut_barcodes]
write.table(CD_27_Glut_6hr_barcodes,
            file = "line_type_barcodes/CD_27_Glut_6hr_barcodes.txt", 
            quote = F, sep = "\t", row.names = F, col.names = F)



file_list <- list.files(path = "cell_line_barcodes_copy", 
                        recursive = F, 
                        include.dirs = F)
for (f in file_list){
  bc <- read_csv(paste("cell_line_barcodes_copy/", f, sep = ""), col_names = FALSE)
  bc <- unlist(bc)
  subsetted_bc_glut <- bc[bc %in% Glut_barcodes]
  filename_glut <- paste("line_type_barcodes/", str_sub(f, start = -5), "_", str_sub(f, 5, 5), "hr_glut_barcodes.txt", # format: line_xhr_glut/GABA_barcodes.txt
                        sep = "")
  write.table(subsetted_bc_glut,
              file = filename_glut,
              quote = F, sep = "\t", row.names = F, col.names = F)
  
  subsetted_bc_GABA <- bc[bc %in% GABA_barcodes]
  filename_GABA <- paste("line_type_barcodes/", str_sub(f, start = -5), "_", str_sub(f, 5, 5), "hr_GABA_barcodes.txt", # format: line_xhr_glut/GABA_barcodes.txt
                         sep = "")
  write.table(subsetted_bc_GABA,
              file = filename_GABA,
              quote = F, sep = "\t", row.names = F, col.names = F)
  # print(filename_glut)
}
