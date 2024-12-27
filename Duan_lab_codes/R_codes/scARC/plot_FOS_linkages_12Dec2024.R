# remove lines 02, 05, 06 and calculate again ####

# init ####
{
  library(Signac)
  library(Seurat)
  library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  
  library(stringr)
  library(readr)
  library(dplyr)
  
  library(viridis)
  library(graphics)
  library(ggplot2)
  library(RColorBrewer)
  
  library(ggforce)
}


plan("multisession", workers = 4)
set.seed(42)
options(future.globals.maxSize = 229496729600)



gene_ranges <- EnsDb.Hsapiens.v86
ensembldb::seqlevelsStyle(gene_ranges) <- "UCSC"
annotation_ranges <-
  GetGRangesFromEnsDb(ensdb = gene_ranges,
                      verbose = T, 
                      standard.chromosomes = T)


unique(multiomic_obj$cell.line.ident)
removed_lines <- subset(multiomic_obj, 
                        cell.line.ident %in% c("CD_62", "CD_50", "CD_56", "CD_34",
                                               "CD_35", "CD_64", "CD_51", "CD_33",
                                               "CD_66", "CD_63", "CD_65", "CD_60",
                                               "CD_55", "CD_04", "CD_52"))
idents_list <- unique(removed_lines$timextype.ident)
idents_list <- sort(idents_list)
idents_list
celltype_time_sep_list <- vector(mode = "list", 
                                 length = length(idents_list))
obj_with_link_list <- vector(mode = "list", 
                             length = length(idents_list))

# USE FOS
for (i in 1:length(celltype_time_sep_list)){
  print(idents_list[i])
  temp_obj <- subset(removed_lines, subset = timextype.ident == idents_list[i])
  # first compute the GC content for each peak
  DefaultAssay(temp_obj) <- "peaks"
  temp_obj <- RegionStats(temp_obj, 
                          genome = BSgenome.Hsapiens.UCSC.hg38)
  print("LinkPeaks() running")
  # link peaks to genes
  obj_with_link_list[[i]] <- LinkPeaks(
    object = temp_obj,
    peak.assay = "peaks",
    expression.assay = "SCT",
    score_cutoff = 0, 
    pvalue_cutoff = 0.05,
    genes.use = "FOS")
  print("LinkPeaks() finished")
}


# run from here !!!!! ######
# install ggforce
# save(list = c("obj_with_link_list",
#               "annotation_ranges"),
#      file = "plot_FOS_RData_12Dec2024.RData")
idents_list <-
  c("GABA_0hr",
    "GABA_1hr",
    "GABA_6hr")

for (i in 1:3) {
  Annotation(obj_with_link_list[[i]]) <- annotation_ranges # this step is run
  #because in the initial pass, the annotation was not working, and was only 
  #detected when plotting coverage plots
  filename <- paste0("./link_peak_to_gene_plots/", idents_list[i], 
                     "_FOS_linkPeak_plot.pdf")
  print(filename)
  Idents(obj_with_link_list[[i]]) <- "cell.type"
  pdf(file = filename, height = 4, width = 4)
  p <- CoveragePlot(
    object = obj_with_link_list[[i]],
    region = "FOS",
    features = "FOS",
    expression.assay = "SCT",
    extend.upstream = 100000,
    extend.downstream = 100000
  )
  print(p)
  dev.off()
}

idents_list <-
  c("npglut_0hr",
    "npglut_1hr",
    "npglut_6hr")

for (i in 7:9) {
  Annotation(obj_with_link_list[[i]]) <- annotation_ranges # this step is run
  #because in the initial pass, the annotation was not working, and was only 
  #detected when plotting coverage plots
  filename <- paste0("./link_peak_to_gene_plots/", idents_list[i-6], 
                     "_FOS_linkPeak_plot.pdf")
  print(filename)
  Idents(obj_with_link_list[[i]]) <- "cell.type"
  pdf(file = filename, height = 4, width = 4)
  p <- CoveragePlot(
    object = obj_with_link_list[[i]],
    region = "FOS",
    features = "FOS",
    expression.assay = "SCT",
    extend.upstream = 1000,
    extend.downstream = 2000
  )
  print(p)
  dev.off()
}

# use different number of sample that makes up the background peaks ####
idents_list <- unique(multiomic_obj$timextype.ident)
idents_list <- sort(idents_list)
idents_list
celltype_time_sep_list <- vector(mode = "list", 
                                 length = length(idents_list))
obj_with_link_list <- vector(mode = "list", 
                             length = length(idents_list))
for (i in 1:length(celltype_time_sep_list)){
  print(idents_list[i])
  temp_obj <- subset(multiomic_obj, subset = timextype.ident == idents_list[i])
  # first compute the GC content for each peak
  DefaultAssay(temp_obj) <- "peaks"
  temp_obj <- RegionStats(temp_obj, 
                          genome = BSgenome.Hsapiens.UCSC.hg38)
  print("LinkPeaks() running")
  # link peaks to genes
  obj_with_link_list[[i]] <- LinkPeaks(
    object = temp_obj,
    peak.assay = "peaks",
    expression.assay = "SCT",
    genes.use = c("BDNF", "VGF"),
    n_sample = 100)
  print("LinkPeaks() finished")
}

for (i in 1:length(obj_with_link_list)){
  Annotation(obj_with_link_list[[i]]) <- annotation_ranges # this step is run
  #because in the initial pass, the annotation was not working, and was only 
  #detected when plotting coverage plots
  filename <- paste0("./link_peak_to_gene_plots/nsample_100/", idents_list[i], 
                     "_linkPeak_plot.pdf")
  print(filename)
  Idents(obj_with_link_list[[i]]) <- "cell.type"
  pdf(file = filename, height = 4, width = 8)
  p <- CoveragePlot(
    object = obj_with_link_list[[i]],
    region = "BDNF",
    features = "BDNF",
    expression.assay = "SCT",
    extend.upstream = 100000,
    extend.downstream = 100000
  )
  print(p)
  dev.off()
}
View(Links(obj_with_link_list[[5]])[str_detect(Links(obj_with_link_list[[5]])$peak, "2777")])

