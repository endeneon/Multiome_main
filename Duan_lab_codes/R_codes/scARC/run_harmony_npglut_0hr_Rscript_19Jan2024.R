# Siwei 18 Jan 2024
# simply run an rscript
# try to run DA analysis and link Lexi's 76 line RNASeq data to ATACSeq data

# init ####
library(Seurat)
library(Signac)

# library(EnsDb.Hsapiens.v86)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)

# library(stringr)
# library(plyranges)

###
load("DA_seurat_peaks_19Jan2023_w_results_harmony.RData")

seurat_npglut_0hr_Harmony_backup <-
  seurat_npglut_0hr_Harmony

## restart from here
# seurat_npglut_0hr_Harmony <-
#   ScaleData(seurat_npglut_0hr)

# rm(data_matrix_raw)

seurat_npglut_0hr_Harmony <-
  RunPCA(seurat_npglut_0hr_Harmony,
         seed.use = 42,
         verbose = T)
# > colnames(seurat_npglut_0hr@meta.data)
# [1] "orig.ident"            "nCount_RNA"            "nFeature_RNA"
# [4] "nCount_ATAC"           "nFeature_ATAC"         "group.ident"
# [7] "cell.line.ident"       "nucleosome_signal"     "nucleosome_percentile"
# [10] "TSS.enrichment"        "TSS.percentile"        "ATAC_snn_res.0.5"
# [13] "seurat_clusters"       "RNA.cell.type"         "RNA.fine.cell.type"
# [16] "time.ident"            "timextype.ident"       "nCount_peaks"
# [19] "nFeature_peaks"        "aff"

library(future)
library(harmony)


{
  plan("multisession", workers = 8)
  set.seed(42)
  options(future.globals.maxSize = 229496729600)
}


seurat_npglut_0hr_Harmony <-
  RunHarmony(seurat_npglut_0hr_Harmony,
             group.by.vars = c("group.ident",
                               "TSS.percentile",
                               "nucleosome_percentile"))

save.image("DA_seurat_peaks_19Jan2023_w_results_harmony.RData")
