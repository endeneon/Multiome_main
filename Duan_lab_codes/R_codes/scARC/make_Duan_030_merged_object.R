# Siwei 05 Sept 2024
# Count TSS and NS of Duan-030

# init ####
{
  library(Seurat)
  library(Signac)
  library(sctransform)
  library(glmGamPoi)
  library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
  # library(GenomicFeatures)
  library(AnnotationDbi)
  # library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)

  library(patchwork)
  library(readr)
  library(ggplot2)
  library(RColorBrewer)
  library(dplyr)
  library(viridis)
  library(graphics)
  library(stringr)
  library(future)

}

plan("multisession", workers = 4)
set.seed(42)
options(future.globals.maxSize = 229496729600)

# assemble file list ####
# h5list <-
#   list.files(path = "~/Data/FASTQ/Duan_Project_030",
#              pattern = "*filtered_feature_bc_matrix.h5",
#              full.names = T,
#              recursive = T)
per_barcode_metrics <-
  list.files(path = "~/Data/FASTQ/Duan_Project_030",
             pattern = "per_barcode_metrics.csv$",
             full.names = T,
             recursive = T)
atac_fraglist <-
  list.files(path = "~/Data/FASTQ/Duan_Project_030",
             pattern = "*atac_fragments.tsv.gz$",
             full.names = T,
             recursive = T)
gex_molecule_info <-
  list.files(path = "~/Data/FASTQ/Duan_Project_030",
             pattern = "gex_molecule_info.h5$",
             full.names = T,
             recursive = T)
lib_names <-
  str_split(string = atac_fraglist,
            pattern = "\\/",
            simplify = T)[, 9]
lib_names <-
  str_split(string = lib_names,
            pattern = '_',
            simplify = T)[, 3]

df_4_cellranger_arc <-
  data.frame(library_id = lib_names,
             atac_fragments = atac_fraglist,
             per_barcode_metrics = per_barcode_metrics,
             gex_molecule_info = gex_molecule_info)

write.table(df_4_cellranger_arc,
            file = "Duan_030_merge_all_libs.csv",
            quote = F, sep = ",",
            row.names = F, col.names = T)



# get gene annotations #####
annotation_hg38 <-
  GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation_hg38) <-
  paste0("chr",
         seqlevels(annotation_hg38))

# make ATAC project for merging #####
# master_atac_obj_list <-
#   vector(mode = "list",
#          length = length(fraglist))
# names(master_atac_obj_list) <-
#   lib_names
#
# for (i in 1:length(fraglist)) {
#   print(i)
#   raw_counts <-
#     Read10X_h5(h5list[i])
#   master_atac_obj_list[[i]] <-
#     CreateChromatinAssay(counts = raw_counts$Peaks,
#                          sep = c(":",
#                                  "-"),
#                          fragments = fraglist[i],
#                          # genome = "hg38",
#                          annotation = annotation_hg38,
#                          verbose = T)
# }
