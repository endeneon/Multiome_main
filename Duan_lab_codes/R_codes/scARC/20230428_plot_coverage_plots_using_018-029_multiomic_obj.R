# Chuxuan Li 08/05/2022
# plot the peaks called using single cell method, then also for peaks called by
#Siwei using pseudobulk method
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
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/multiomic_obj_with_new_peaks_labeled.RData")

GABA <- subset(multiomic_obj_new, RNA.cell.type == "GABA")
nmglut <- subset(multiomic_obj_new, RNA.cell.type == "nmglut")
npglut <- subset(multiomic_obj_new, RNA.cell.type == "npglut")

#load("../R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")
load("./EnsDb_UCSC_hg38_annotation_new.RData")

Annotation(npglut) <- ens_use


# FOS ####
CoveragePlot(
  object = npglut,
  region = "chr14-75221000-75290000", # FOS
  # region.highlight = GRanges(seqnames = "chr11",
  #                            ranges = IRanges(start = 66418000, end = 66418700)),
  group.by = "time.ident")

CoveragePlot(
  object = npglut,
  region = "chr11-66414647-66440411", # NPAS4
  group.by = "timextype.ident")

CoveragePlot(
  object = npglut,
  region = "chr5-96310000-96520000", # PCSK1
  group.by = "timextype.ident")

CoveragePlot(
  object = npglut,
  region = "chr7-21100000-21530000", # SP4
  group.by = "timextype.ident")

CoveragePlot(
  object = npglut,
  region = "chr11-27630000-27775000", # BDNF
  group.by = "timextype.ident")

# plot cell types ceparately, time points separately ####
load("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/multiomic_obj.RData")
gene_ranges <- EnsDb.Hsapiens.v86
ensembldb::seqlevelsStyle(gene_ranges) <- "UCSC"
annotation_ranges <-
  GetGRangesFromEnsDb(ensdb = gene_ranges,
                      verbose = T, 
                      standard.chromosomes = T)

idents_list <- unique(multiomic_obj$timextype.ident)
idents_list <- sort(idents_list)
names(obj_with_link_list) <- idents_list
DefaultAssay(obj_with_link_list[[1]])
unique(Idents(obj_with_link_list[[1]]))

for (i in 1:length(obj_with_link_list)) {
  Idents(obj_with_link_list[[i]]) <- str_replace_all(Idents(obj_with_link_list[[i]]),
                                                     "_", " ")
  filename <- paste0("./link_peak_to_gene_jubao_genes/updated_plots_08312022/", 
                     names(obj_with_link_list)[i], "_DIP2A_coverage_plot.pdf")
  pdf(file = filename, width = 8.98, height = 5.08)
  p <- CoveragePlot(
    object = obj_with_link_list[[i]],
    region = "DIP2A",
    features = "DIP2A",
    expression.assay = "SCT", 
    region.highlight = StringToGRanges(c("chr21-46509746-46510304",
                                         "chr21-46410983-46411607",
                                         "chr21-46412422-46413653")), 
    max.downsample = ncol(obj_with_link_list[[i]]), ymax = 1000, 
    extend.upstream = 62470,
    extend.downstream = 21000
  )
  print(p)
  dev.off()
}

for (i in 1:length(obj_with_link_list)) {
  Idents(obj_with_link_list[[i]]) <- str_replace_all(Idents(obj_with_link_list[[i]]),
                                                     "_", " ")
  filename <- paste0("./link_peak_to_gene_jubao_genes/updated_plots_08312022/", 
                     names(obj_with_link_list)[i], "_GRIN2A_coverage_plot.pdf")
  pdf(file = filename, width = 8.98, height = 5.08)
  p <- CoveragePlot(
    object = obj_with_link_list[[i]],
    region = "GRIN2A",
    features = "GRIN2A",
    expression.assay = "SCT", 
    region.highlight = StringToGRanges(c("chr16-9735960-9737200",
                                         "chr16-10079470-10081122")), 
    max.downsample = ncol(obj_with_link_list[[i]]), ymax = 400, 
    extend.upstream = 100000,
    extend.downstream = 100000
  )
  
  print(p)
  dev.off()
}

CoveragePlot(
  object = obj_with_link_list[[6]],
  region = "GRIN2A",
  features = "GRIN2A",
  expression.assay = "SCT", 
  region.highlight = StringToGRanges(c("chr16-9735960-9737200",
                                       "chr16-10079470-10081122")), 
  max.downsample = ncol(obj_with_link_list[[9]]), #ymax = 500, 
  extend.upstream = 100000,
  extend.downstream = 100000
)

CoveragePlot(
  object = obj_with_link_list[[6]],
  region = "DIP2A",
  features = "DIP2A",
  expression.assay = "SCT", 
  region.highlight = StringToGRanges(c("chr21-46509746-46510304",
                                       "chr21-46410983-46411607",
                                       "chr21-46412422-46413653")), 
  max.downsample = ncol(obj_with_link_list[[6]]), ymax = 1000, 
  extend.upstream = 62470,
  extend.downstream = 21000
)
names(obj_with_link_list)
l <- Links(obj_with_link_list[[7]])
l[l$gene == "GRIN2A",]

# plot (normalize) by cell type ####

GABA <- subset(multiomic_obj, subset = cell.type == "GABA")
npglut <- subset(multiomic_obj, subset = cell.type == "NEFM_pos_glut")
nmglut <- subset(multiomic_obj, subset = cell.type == "NEFM_neg_glut")

CoveragePlot(
  object = GABA,
  region = "GRIN2A",
  features = "GRIN2A", annotation = T,
  expression.assay = "SCT", 
  region.highlight = StringToGRanges(c("chr16-9735961-9737093",
                                       "chr16-9787990-9788950",
                                       "chr16-10079470-10081122")),
  assay.scale = "common", group.by = "timextype.ident", 
  ymax = 700, heights = c(4, 1, 1),
  extend.upstream = 100000,
  extend.downstream = 100000
)
Links(obj_with_link_list[[1]])

Annotation(nmglut)
CoveragePlot(
  object = nmglut,
  region = "chr16-9750000-9800000",
  #features = "GRIN2A", annotation = T,
  expression.assay = "SCT", assay = "peaks",
  #region.highlight = StringToGRanges(c(#"chr16-9735961-9737093",
  #"chr16-9787990-9788950",
  #                                     "chr16-10079470-10081122")),
  assay.scale = "common", group.by = "timextype.ident", 
  ymax = 700, #downsample.rate = 1#, 
  #extend.upstream = 100000,
  #extend.downstream = 100000
)
#rs16966337: chr16-9788096
CoveragePlot(
  object = npglut,
  region = "GRIN2A",
  features = "GRIN2A", annotation = T, assay = "peaks",
  expression.assay = "SCT", 
  region.highlight = StringToGRanges(c("chr16-9735961-9737093",
                                       "chr16-9787990-9788950",
                                       "chr16-10079470-10081122")),
  assay.scale = "common", group.by = "timextype.ident", 
  ymax = 700, heights = c(4, 1, 1),
  extend.upstream = 100000,
  extend.downstream = 100000
)
CoveragePlot(
  object = GABA,
  region = "DIP2A",
  features = "DIP2A",
  expression.assay = "SCT", assay = "peaks",
  region.highlight = StringToGRanges(c("chr21-46509746-46510304",
                                       "chr21-46410983-46411607",
                                       "chr21-46412422-46413653")), 
  max.downsample = ncol(npglut), group.by = "timextype.ident", ymax = 1000, 
  extend.upstream = 62470,
  extend.downstream = 21000
)
