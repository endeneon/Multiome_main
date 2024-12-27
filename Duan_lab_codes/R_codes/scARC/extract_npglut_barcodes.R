# Siwei 12 Feb 2024
# extract npglut 016 cell barcodes

load("~/NVME/scARC_Duan_018/018-029_combined_analysis/df_Lexi_cells_470K.RData")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/metadata_multiomic_obj_470K_cells.RData")

npglut_barcodes <-
  rownames(metadata_multiomic_obj[metadata_multiomic_obj$RNA.cell.type == "npglut", ])
save(npglut_barcodes,
     file = "df_npglut_cells_62K.RData")


## extract npglut & GABA 0hr barcodes
npglut_GABA_0hr_barcodes <-
  rownames(metadata_multiomic_obj[(metadata_multiomic_obj$RNA.cell.type %in% c("npglut", "GABA")) &
                                    (metadata_multiomic_obj$time.ident == "0hr"), ])
save(npglut_GABA_0hr_barcodes,
     file = "df_0hr_npglut_GABA_cells_50K.RData")
