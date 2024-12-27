# Siwei Zhang 22 Nov 2024

# plot the peaks called using single cell method, then also for peaks called by
# Siwei using pseudobulk method
# plot HDAC9

{
  library(Signac)
  library(Seurat)
  library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  
  library(ggplot2)
  library(stringr)
}

plan("multisession", workers = 4)
set.seed(42)
options(future.globals.maxSize = 229496729600)



load("~/NVME/scARC_Duan_018/018-029_combined_analysis/multiomic_obj_with_new_peaks_labeled.RData")
load("./EnsDb_UCSC_hg38_annotation_new.RData")

# Annotation(multiomic_obj_new) <- ens_use

gene_ranges <- EnsDb.Hsapiens.v86
ensembldb::seqlevelsStyle(gene_ranges) <- "UCSC"
annotation_ranges <-
  GetGRangesFromEnsDb(ensdb = gene_ranges,
                      verbose = T, 
                      standard.chromosomes = T)

Annotation(multiomic_obj_new) <-
  annotation_ranges
DefaultAssay(multiomic_obj_new) <- "RNA"
multiomic_obj_new <-
  SCTransform(multiomic_obj_new,
              assay = "RNA", 
              new.assay.name = "SCT", 
              seed.use = 42,
              verbose = T)


group_78_idents <-
  sort(unique(multiomic_obj_new$timextype.ident))[1:9]


# HDAC9 ####
Idents(multiomic_obj_new) <- "timextype.ident"

CoveragePlot(
  object = multiomic_obj_new,
  # region = "HDAC9",
  region = "chr7-17600000-18600000",
  features = "HDAC9",
  expression.assay = "RNA", 
  # idents = sort(unique(multiomic_obj_new$timextype.ident))[1:9], 
  # region.highlight = StringToGRanges(c("chr21-46509746-46510304",
  #                                      "chr21-46410983-46411607",
  #                                      "chr21-46412422-46413653")), 
  # max.downsample = ncol(obj_with_link_list[[i]]), 
  # ymax = 1000, 
  extend.upstream = 0,
  extend.downstream = 0
)

CoveragePlot(
  object = multiomic_obj_new,
  # region = "chr7-18086825-19002416", # HDAC9
  region = "HDAC9", # HDAC9
  features = "HDAC9",
  expression.assay = "RNA",
  idents = sort(unique(multiomic_obj_new$timextype.ident))[1:9],
  extend.upstream = 1000,
  extend.downstream = 1000,
  
  group.by = "timextype.ident")

sort(unique(multiomic_obj_new$timextype.ident))[1:9]

# FOS ####
npglut <- subset(multiomic_obj_new, RNA.cell.type == "npglut")
CoveragePlot(
  object = npglut,
  region = "chr14-75221000-75290000", # FOS
  # region.highlight = GRanges(seqnames = "chr11",
  #                            ranges = IRanges(start = 66418000, end = 66418700)),
  group.by = "time.ident")







DefaultAssay(multiomic_obj_new)
multiomic_obj_new[["ATAC"]]@annotation
Annotation(multiomic_obj_new) <-
  annotation_genome

CoveragePlot(
  object = multiomic_obj_new,
  # region = "chr7-18086825-19002416", # HDAC9
  region = "MS4A1", # HDAC9
  features = "MS4A1",
  expression.assay = "RNA",
  idents = sort(unique(multiomic_obj_new$timextype.ident))[1:9],
  extend.upstream = 1000,
  extend.downstream = 1000,
  group.by = "timextype.ident",
  sep = c('-', '-'),
  annotation = T)


#####
load("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/multiomic_obj.RData")
gene_ranges <- EnsDb.Hsapiens.v86
ensembldb::seqlevelsStyle(gene_ranges) <- "UCSC" 



annotation_ranges <-
  GetGRangesFromEnsDb(ensdb = gene_ranges,
                      verbose = T, 
                      standard.chromosomes = T)

ident_list <-
  sort(unique(multiomic_obj$timextype.ident))
Annotation(multiomic_obj) <-
  annotation_ranges

Idents(multiomic_obj) <- "timextype.ident"
# multiomic_obj$timextype.ident <-
Idents(multiomic_obj) <-
  factor(Idents(multiomic_obj),
         levels = sort(levels(multiomic_obj)))

# levels(multiomic_obj) <-
#   sort(unique(multiomic_obj$timextype.ident))



CoveragePlot(
  object = multiomic_obj,
  # region = "HDAC9",
  # region = "chr7-17600000-18600000",
  region = "chr7-17586948-18586949",
  features = "HDAC9",
  expression.assay = "SCT", 
  idents = sort(unique(multiomic_obj$timextype.ident)), 
  # region.highlight = StringToGRanges(c("chr21-46509746-46510304",
  #                                      "chr21-46410983-46411607",
  #                                      "chr21-46412422-46413653")), 
  # max.downsample = ncol(obj_with_link_list[[i]]), 
  # ymax = 1000, 
  extend.upstream = 0,
  extend.downstream = 0
)
