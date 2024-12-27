# Chuxuan Li 04/15/2022
# Use the 18-line 10x-aggred ATACseq object with new MACS2 called peaks,
#perform link peak to genes analysis and plot resulting links in Coverage plots

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
}


plan("multisession", workers = 4)
set.seed(42)
options(future.globals.maxSize = 229496729600)



load("./ATAC_new_obj_with_MACS2_new_peaks_and_RNA_labels.RData")
load("../Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")


# combine RNA and ATAC objects ####
unique(ATAC_new$timextype.ident)
ATAC_for_links <- subset(ATAC_new, timextype.ident %in% c("GABA_0hr", "GABA_6hr",
                                                          "NEFM_neg_glut_0hr", "GABA_1hr",
                                                          "NEFM_pos_glut_1hr", "NEFM_neg_glut_1hr",
                                                          "NEFM_neg_glut_6hr", "NEFM_pos_glut_6hr",
                                                          "NEFM_pos_glut_0hr"))
unique(integrated_labeled$cell.type)
RNA_for_links <- subset(integrated_labeled, cell.type != "unknown")
colnames(ATAC_for_links)
unique(str_extract(colnames(ATAC_for_links), "-[0-9]+"))

unique(str_extract(colnames(RNA_for_links), "_[0-9]+"))
unique(str_extract(colnames(RNA_for_links), "[0-9]+_"))
RNA_for_links <- RenameCells(object = RNA_for_links, 
                             new.names = str_remove(colnames(RNA_for_links), "1_"))
unique(str_extract(colnames(RNA_for_links), "-[0-9]+"))
colnames(RNA_for_links)

# match ATAC and RNA cell barcodes
ATAC_for_links$cells.in.RNA <- "F" 
ATAC_for_links$cells.in.RNA[colnames(ATAC_for_links) %in% colnames(RNA_for_links)] <- "T"
unique(ATAC_for_links$cells.in.RNA)
ATAC_for_links <- subset(ATAC_for_links, cells.in.RNA == "T")
RNA_for_links$cells.in.ATAC <- "F" 
RNA_for_links$cells.in.ATAC[colnames(RNA_for_links) %in% colnames(ATAC_for_links)] <- "T"
unique(RNA_for_links$cells.in.ATAC)
RNA_for_links <- subset(RNA_for_links, cells.in.ATAC == "T")

# combine RNA and ATAC
multiomic_obj <- ATAC_for_links
multiomic_obj[["SCT"]] <- RNA_for_links@assays$SCT


# coverage plot ####
gene_ranges <- EnsDb.Hsapiens.v86
ensembldb::seqlevelsStyle(gene_ranges) <- "UCSC"
annotation_ranges <-
  GetGRangesFromEnsDb(ensdb = gene_ranges,
                      verbose = T, 
                      standard.chromosomes = T)
Annotation(multiomic_obj[["peaks"]]) <- annotation_ranges

idents <- unique(multiomic_obj$cell.type)
idents
DefaultAssay(multiomic_obj) <- "peaks"
Annotation(multiomic_obj[["peaks"]])

# 
# for (i in 1:length(idents)){
#   celltype <- idents[i]
#   print(celltype)
#   file_name <- paste0("./peak_track_plots/NPAS4_",
#                       celltype,
#                       "_peak_tracks.pdf")
#   pdf(file = file_name)
#   p <- CoveragePlot(
#     object = subset(ATAC_new, subset = cell.type == idents[i]),
#     region = "chr11-66414647-66427411", # NPAS4
#     region.highlight = GRanges(seqnames = "chr11",
#                                ranges = IRanges(start = 66418000, end = 66418700)),
#     group.by = "time.ident")
#   
#   print(p)
#   dev.off()
# }

multiomic_obj$timextype.ident.for.plot <- str_replace_all(multiomic_obj$timextype.ident,
                                                     "_", " ")
multiomic_obj$timextype.ident.for.plot <- str_replace(multiomic_obj$timextype.ident.for.plot,
                                                 " pos", "+")
multiomic_obj$timextype.ident.for.plot <- str_replace(multiomic_obj$timextype.ident.for.plot,
                                                 " neg", "-")
unique(multiomic_obj$timextype.ident.for.plot)
CoveragePlot(
  object = multiomic_obj,
  region = "chr11-27689113-27774038", # BDNF
  region.highlight = GRanges(seqnames = "chr11",
                             ranges = IRanges(start = 27767000, end = 27777000)),
  group.by = "timextype.ident.for.plot")
CoveragePlot(
  object = multiomic_obj,
  region = "chr11-66414647-66427411", # NPAS4
  region.highlight = GRanges(seqnames = "chr11",
                             ranges = IRanges(start = 66418000, end = 66418700)),
  group.by = "timextype.ident.for.plot")
CoveragePlot(
  object = multiomic_obj,
  region = "chr14-75286000-75298000", # FOS
  #region.highlight = GRanges(seqnames = "chr11",
  #                           ranges = IRanges(start = 66418000, end = 66418700)),
  group.by = "timextype.ident.for.plot")
CoveragePlot(
  object = multiomic_obj,
  region = "chr7-101160000-101172000", # VGF
  #region.highlight = GRanges(seqnames = "chr11",
  #                           ranges = IRanges(start = 66418000, end = 66418700)),
  group.by = "timextype.ident.for.plot")


# link peak to genes for each time x type combination ####
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
    expression.assay = "SCT")
  print("LinkPeaks() finished")
}

for (i in 1:length(obj_with_link_list)){
  Annotation(obj_with_link_list[[i]]) <- annotation_ranges # this step is run
  #because in the initial pass, the annotation was not working, and was only 
  #detected when plotting coverage plots
  celltype_time_sep_list[[i]] <- as.data.frame(Links(obj_with_link_list[[i]]))
  write.table(x = as.data.frame(celltype_time_sep_list[[i]]), 
              file = paste0("./link_peak_to_gene_results/",
                            idents_list[i], ".csv"),
              quote = F, sep = ",", col.names = T, row.names = F)
}

# plot peak links of early/late response genes ####
idents_list
# plotted VGF and BDNF using the script below, switching different gene names
for (i in 1:length(obj_with_link_list)){
  filename <- paste0("./link_peak_to_gene_plots/", "BDNF_", idents_list[i], 
                     "_linkPeak_plot.pdf")
  print(filename)
  Idents(obj_with_link_list[[i]]) <- "timextype.ident"
  pdf(file = filename, height = 6, width = 12)
  p <- CoveragePlot(
    object = obj_with_link_list[[i]],
    region = "BDNF",
    features = "BDNF",
    expression.assay = "SCT", 
    extend.upstream = 50000,
    extend.downstream = 150000
  )
  print(p)
  dev.off()
}

# compute links by cell type only ####
idents_list <- unique(multiomic_obj$cell.type)
idents_list
type_obj_with_link_list <- vector(mode = "list", 
                             length = length(idents_list))
for (i in 1:length(type_obj_with_link_list)){
  print(idents_list[i])
  temp_obj <- subset(multiomic_obj, subset = cell.type == idents_list[i])
  # first compute the GC content for each peak
  DefaultAssay(temp_obj) <- "peaks"
  temp_obj <- RegionStats(temp_obj, 
                          genome = BSgenome.Hsapiens.UCSC.hg38)
  print("LinkPeaks() running")
  # link peaks to genes
  type_obj_with_link_list[[i]] <- LinkPeaks(
    object = temp_obj,
    peak.assay = "peaks",
    expression.assay = "SCT")
  print("LinkPeaks() finished")
}

for (i in 1:length(type_obj_with_link_list)){
  Annotation(type_obj_with_link_list[[i]]) <- annotation_ranges # this step is run
  #because in the initial pass, the annotation was not working, and was only 
  #detected when plotting coverage plots
  celltype_time_sep_list[[i]] <- as.data.frame(Links(type_obj_with_link_list[[i]]))
  write.table(x = as.data.frame(celltype_time_sep_list[[i]]), 
              file = paste0("./link_peak_to_gene_results/sep_by_type_",
                            idents_list[i], ".csv"),
              quote = F, sep = ",", col.names = T, row.names = F)
}

save(type_obj_with_link_list, file = "type_obj_with_link_list.RData")

# plot cell type specific only peak links ####
idents_list
# plotted VGF and BDNF using the script below, switching different gene names
for (i in 1:length(type_obj_with_link_list)){
  filename <- paste0("./link_peak_to_gene_plots/", "VGF_", idents_list[i], 
                     "_linkPeak_plot.pdf")
  print(filename)
  Idents(type_obj_with_link_list[[i]]) <- "cell.type"
  pdf(file = filename, height = 4, width = 8)
  p <- CoveragePlot(
    object = type_obj_with_link_list[[i]],
    region = "VGF",
    features = "VGF",
    expression.assay = "SCT",
    extend.upstream = 30000,
    extend.downstream = 30000
  )
  print(p)
  dev.off()
}

# combine two glut groups and call linkPeaks() again ####
two_gluts <- subset(multiomic_obj, RNA.cell.type %in% c("NEFM_pos_glut", "NEFM_neg_glut"))
rownames(two_gluts@assays$peaks)
counts <- GetAssayData(two_gluts, assay = "peaks")
counts <- counts[str_detect(rownames(two_gluts@assays$peaks), "chr11"), ]
two_gluts <- subset(two_gluts, features = rownames(counts))
multiomic_obj$cells.in.ATAC <- "F" 
multiomic_obj$cells.in.ATAC[colnames(multiomic_obj) %in% colnames(two_gluts)] <- "T"
unique(multiomic_obj$cells.in.ATAC)
multiomic_obj <- subset(multiomic_obj, cells.in.ATAC == "T")

two_gluts <- RegionStats(two_gluts, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
two_gluts[["SCT"]] <- multiomic_obj@assays$SCT
two_gluts <- LinkPeaks(
  object = two_gluts,
  peak.assay = "peaks",
  expression.assay = "SCT")
Links(two_gluts)

# link peak to genes no score cutoff for each time x type combination ####
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
    score_cutoff = 0, 
    pvalue_cutoff = 0.5,
    genes.use = "BDNF")
  print("LinkPeaks() finished")
}

for (i in 1:length(obj_with_link_list)){
  Annotation(obj_with_link_list[[i]]) <- annotation_ranges # this step is run
  #because in the initial pass, the annotation was not working, and was only 
  #detected when plotting coverage plots
  filename <- paste0("./link_peak_to_gene_plots/", "no_score_cutoff_pval_0.5_BDNF_", idents_list[i], 
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


# remove lines 02, 05, 06 and calculate again ####
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
    genes.use = "BDNF")
  print("LinkPeaks() finished")
}

for (i in 1:length(obj_with_link_list)){
  Annotation(obj_with_link_list[[i]]) <- annotation_ranges # this step is run
  #because in the initial pass, the annotation was not working, and was only 
  #detected when plotting coverage plots
  filename <- paste0("./link_peak_to_gene_plots/BDNF_removed_3lines/", idents_list[i], 
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

