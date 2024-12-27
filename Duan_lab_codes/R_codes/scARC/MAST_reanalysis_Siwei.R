# Siwei 22 Oct 2024
# Use MAST to do single cell DE analysis with 018-029 data combined

# init ####
{
  library(ggplot2)
  library(ggrepel)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
  library(Seurat)
}

# make an RData for 100 line metadata csv ####
raw_100line_meta <-
  readr::read_csv("100line_metadata_df.csv")
load("018-030_RNA_integrated_labeled_with_harmony.RData")
# load("/nvmefs/scARC_Duan_018/case_control_DE_analysis/018-029_covar_table_final.RData")

#options(mc.cores = 30)
