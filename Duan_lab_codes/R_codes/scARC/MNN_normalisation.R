# Siwei 04 March 2022
# Use MNN method to normalise 6 libraries of RNA-seq


# init
library(SingleCellExperiment)
library(scater)
library(scMerge)
library(BiocParallel)
library(BiocSingular)

library(Seurat)
library(Signac)

library(batchelor)
library(scran)


library(stringr)
library(parallel)

options(mc.cores = 16)
### load data

# setwd("/nvmefs/scARC_Duan_018/R_pseudo_bulk_RNA_seq")
### load data
# load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/raw_data_seurat_list.RData")
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/raw_data_seurat_list.RData")

# data("segList", package = "scMerge")

set.seed(42)
i <- 1L
sce_list <- vector(mode = "list", length = 6L)
names(sce_list) <- c("g_2_0", "g_2_1", "g_2_6",
                    "g_8_0", "g_8_1", "g_8_6")

# for (i in i:length(sce_list)) {
#   print(i)
#   sce_list[[i]][["sce_raw"]] <- 
#     as.SingleCellExperiment(obj_lst[[i]])
#   # logNormCount all.sce
#   sce_list[[i]][["all.sce"]] <-
#     lapply(sce_list[[i]][["sce_raw"]], logNormCounts)
#   sce_list[[i]][["all.dec"]] <- 
#     lapply(sce_list[[i]][["all.sce"]], modelGeneVar)
#   sce_list[[i]][["all.hvgs"]] <-
#     lapply(sce_list[[i]][["all.dec"]], getTopHVGs, prop = 0.1)
# }
name_list <- c("g_2_0", "g_2_1", "g_2_6",
               "g_8_0", "g_8_1", "g_8_6")

for (i in 1:length(sce_list)) {
  print(i)
  sce_list[[i]] <- 
    as.SingleCellExperiment(obj_lst[[i]])
  sce_list[[i]]$batch <- name_list[i]
  # # logNormCount all.sce
  # sce_list[[i]][["all.sce"]] <-
  #   lapply(sce_list[[i]][["sce_raw"]], logNormCounts)
  # sce_list[[i]][["all.dec"]] <- 
  #   lapply(sce_list[[i]][["all.sce"]], modelGeneVar)
  # sce_list[[i]][["all.hvgs"]] <-
  #   lapply(sce_list[[i]][["all.dec"]], getTopHVGs, prop = 0.1)
}

# for (i in i:length(sce_list)) {
#   sce_list[[i]]$batch <- name_list[i]
# }
# library size normalisation
# sce_list <-
#   librarySizeFactors(sce_list,
#                      BPPARAM = MulticoreParam(16))

# normalisation
sce_list <- 
  mclapply(sce_list, logNormCounts,
           mc.cores = 16)

# dec_list <-
#   mclapply(sce_list, modelGeneVar,
#            mc.cores = 16)
dec_list <-
  lapply(sce_list, modelGeneVar)

hvgs_list <-
  lapply(dec_list, getTopHVGs, prop = 0.1)

# sce_list_normalise <-
#   do.call(multiBatchNorm, sce_list,
#           list(preserve.single = T,
#                normalize.all = T,
#                BPPARAM = MulticoreParam(16)))
sce_list_normalise <-
  multiBatchNorm(sce_list,
                 normalize.all = T,
                 assay.type = "counts",
                 BPPARAM = MulticoreParam(16))


# combined_dec <-
#   do.call(combineVar, sce_list_normalise)
combined_dec <-
  combineVar(dec_list)
# combined_dec <-
#   combineVar(dec_list[[1]], dec_list[[2]], dec_list[[3]],
#              dec_list[[4]], dec_list[[5]], dec_list[[6]],
#              BPPARAM = MulticoreParam(16))

combined_hvg <-
  getTopHVGs(combined_dec,
             n = 20000)

## perform the merging+normalisation here ##
mnn_out <- 
  fastMNN(sce_list_normalise,
          subset.row = combined_hvg,
          BSPARAM = RandomParam(),
          BPPARAM = MulticoreParam(16))
# mnn_out <- do.call(fastMNN,
#                    c(sce_list_normalise,
#                      list(subset.row = combined_hvg,
#                           BSPARAM = RandomParam())))
# i <- 1L
# for (i in 1:length(mnn_out)) {
#   mnn_out[[i]] <-
#     fastMNN(sce_list_normalise[[i]],
#             subset.row = combined_hvg,
#             BSPARAM = BiocSingular::RandomParam(),
#             BPPARAM = MulticoreParam(16))
# }
save(mnn_out,
     file = "mnn_normalised_20K_output.RData",
     compress = T)

save.image(file = "MNN_normalisation.RData")


# cannot use MultiCoreParam here since this will deconvolute the list
# hence found the "no genes found" error
quick_corrected <-
  quickCorrect(sce_list$g_2_0, sce_list$g_2_1, sce_list$g_2_6,
               sce_list$g_8_0, sce_list$g_8_1, sce_list$g_8_6,
               precomputed = list(dec_list$g_2_0, dec_list$g_2_1, dec_list$g_2_6,
                                  dec_list$g_8_0, dec_list$g_8_1, dec_list$g_8_6),
               # BPPARAM = MulticoreParam(16),
               PARAM = FastMnnParam(BSPARAM = BiocSingular::RandomParam()))

quick_sce <- quick_corrected$corrected

set.seed(42)
quick_sce <- runTSNE(quick_sce, 
                     BPPARAM = MulticoreParam(16),
                     dimred = "corrected")
quick_sce$batch <- factor(quick_sce$batch)
plotTSNE(quick_sce, 
         colour_by = "batch")

### MNN merge and normalisation
combined_dec <-
  combineVar(dec_list$g_2_0, dec_list$g_2_1, dec_list$g_2_6,
             dec_list$g_8_0, dec_list$g_8_1, dec_list$g_8_6)
chosen_hvgs <-
  (combined_dec$bio > 0)

mnn_out <- fastMNN(sce_list$g_2_0, sce_list$g_2_1, sce_list$g_2_6,
                   sce_list$g_8_0, sce_list$g_8_1, sce_list$g_8_6,
                   d = 50, k = 20,
                   subset.row = chosen_hvgs,
                   # BPPARAM = MulticoreParam(16),
                   BSPARAM = BiocSingular::RandomParam(deferred = T))
metadata(mnn_out)$merge.info$lost.var
mnn_out@assays@data$reconstructed




mnn_out <- runTSNE(mnn_out,
                   dimred = "corrected",
                   BPPARAM = MulticoreParam(8))
mnn_out$batch <- factor(mnn_out$batch)
plotTSNE(mnn_out, 
         colour_by = "batch") 
plotTSNE(mnn_out, 
         colour_by = "batch") 
## convert back to Seurat object
mnn_seurat_out <- as.Seurat(mnn_out)


# check the shape of uncorrected matrix
uncorrected <- cbind(sce_list$g_2_0, sce_list$g_2_1, sce_list$g_2_6,
                     sce_list$g_8_0, sce_list$g_8_1, sce_list$g_8_6)

combined_dec <-
  combineVar(dec_list$g_2_0, dec_list$g_2_1, dec_list$g_2_6,
             dec_list$g_8_0, dec_list$g_8_1, dec_list$g_8_6)
chosen_hvgs <-
  (combined_dec$bio > 0)
set.seed(42)
uncorrected <- runPCA(uncorrected, subset_row = chosen_hvgs,
                      BPPARAM = MulticoreParam(16),
                      BSPARAM = BiocSingular::RandomParam())
uncorrected <- runTSNE(uncorrected,
                       dimred = "PCA",
                       BPPARAM = MulticoreParam(8))

plotTSNE(uncorrected, colour_by = "batch")
# rm(list = c("sce_merged",
#             "sce_merged_backup",
#             "scMerge_semi_supervised_merged",
#             "scMerge_semi_supervised_merged_fast",
#             "seurat_merged",
#             "temp_seurat"))

