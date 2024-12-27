# Siwei 2 Dec 2024
# Run MAST use different cell+time combinations

# init ####
{
  library(Seurat)
  library(Signac)
  
  library(ggplot2)
  library(ggrepel)
  library(GGally)
  # library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  # library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(stringr)
  # library(NMF)
  # library(rsvd)
  library(RColorBrewer)
  library(MAST)
  
}

{
  plan("multisession", workers = 4)
  options(expressions = 20000)
  options(future.globals.maxSize = 207374182400)
  
}

# load data and assign aff #####
integrated_56_matched_samples <-
  readRDS("integrated_56_samples_case_ctrl_21Dec2024.RDs")

aff_ident <-
  read.csv("100line_metadata_df.csv")
case_samples <-
  aff_ident$cell_line[aff_ident$aff == "case"]

integrated_56_matched_samples$aff <-
  "ctrl"
integrated_56_matched_samples$aff[integrated_56_matched_samples$cell.line.ident %in% case_samples] <-
  "case"


integrated_56_matched_samples@meta.data

data_re_integration <-
  CreateSeuratObject(counts = integrated_56_matched_samples@assays$RNA$counts)


data_re_integration@meta.data <-
  