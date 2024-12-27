# recall a collection of peaks from all npglut cells
# Siwei 16 Feb 2024

# init ####
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
# library(GenomeInfoDb)
# library(GenomicFeatures)
# library(AnnotationDbi)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

library(RColorBrewer)

library(stringr)
library(future)

plan("multisession", workers = 8)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# load data ####
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/018-029_multiomic_obj_clustered_added_metadata.RData")

DefaultAssay(multiomic_obj)
npglut_cells <-
  colnames(multiomic_obj[, multiomic_obj$RNA.cell.type == "npglut"])
fragments_raw <-
  CreateFragmentObject(path = "./018_to_029_combined/outs/atac_fragments.tsv.gz",
                       cells = npglut_cells,
                       validate.fragments = F,
                       verbose = T)
save(fragments_raw,
     file = "fragments_all_npglut.RData")

peaks_npglut_combined <-
  CallPeaks(object = fragments_raw,
            # assay = "ATAC",
            # group.by = "time.ident",
            macs2.path = "~/Data/Anaconda3-envs/aligners/bin/macs2",
            combine.peaks = F,
            cleanup = T,
            verbose = T,
            additional.args = '-q 0.05 --seed 42')
save(peaks_npglut_combined,
     file = "MACS2_peaks_all_npglut.RData")
