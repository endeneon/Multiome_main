# Chuxuan Li 10/27/2021
# Plot peaks for each big clusters of GABA and glut cells mimicking the plots in the grant 

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)

library(stringr)
library(readr)

unique(aggr_signac_unfiltered$RNA.cell.type.ident)
# group the subtypes into broad categories
aggr_signac_unfiltered$broad.cell.type <- NA
aggr_signac_unfiltered$broad.cell.type[aggr_signac_unfiltered$RNA.cell.type.ident %in% 
                                         c("GAD1+/GAD2+ GABA", "GAD1+ GABA (less mature)", "GAD1+ GABA")] <- "GABA"
aggr_signac_unfiltered$broad.cell.type[aggr_signac_unfiltered$RNA.cell.type.ident %in% 
                                         c("NEFM+/CUX2- glut", "NEFM-/CUX2+ glut",
                                           "NEFM+/CUX2-, SST+ glut", "NEFM-/CUX2+, ADCYAP1+ glut",
                                           "NEFM+/CUX2-, ADCYAP1+ glut", "NEFM-/CUX2+, TCERG1L+ glut",
                                           "NEFM+/CUX2+, SST+ glut")] <- "glut"
aggr_signac_unfiltered$broad.cell.type[aggr_signac_unfiltered$RNA.cell.type.ident %in% 
                                         c("NPC", "immature neuron", "MAP2+ NPC", "subcerebral immature neuron")] <- "NPC"
# divide the glut
aggr_signac_unfiltered$fine.cell.type <- NA
aggr_signac_unfiltered$fine.cell.type[aggr_signac_unfiltered$RNA.cell.type.ident %in% 
                                        c("GAD1+/GAD2+ GABA", "GAD1+ GABA (less mature)", "GAD1+ GABA")] <- "GABA"
aggr_signac_unfiltered$fine.cell.type[aggr_signac_unfiltered$RNA.cell.type.ident %in% 
                                        c("NEFM+/CUX2- glut", "NEFM+/CUX2-, SST+ glut", 
                                          "NEFM+/CUX2-, ADCYAP1+ glut")] <- "NEFM+/CUX2- glut"
aggr_signac_unfiltered$fine.cell.type[aggr_signac_unfiltered$RNA.cell.type.ident %in% 
                                        c("NEFM-/CUX2+, glut",
                                          "NEFM-/CUX2+, ADCYAP1+ glut",
                                          "NEFM-/CUX2+, TCERG1L+ glut"
                                        )] <- "NEFM-/CUX2+ glut"
aggr_signac_unfiltered$fine.cell.type[aggr_signac_unfiltered$RNA.cell.type.ident %in% 
                                        c("NPC", "immature neuron", "MAP2+ NPC", "subcerebral immature neuron")] <- "NPC"

# separate into cell lines
# import cell ident barcodes
times <- c("0", "1", "6")
lines_g2 <- c("CD_27", "CD_54")
lines_g8 <- c("CD_08", "CD_25", "CD_26")
g_2_0_bc_lst <- vector(mode = "list", length = 2)
g_2_1_bc_lst <- vector(mode = "list", length = 2)
g_2_6_bc_lst <- vector(mode = "list", length = 2)
g_8_0_bc_lst <- vector(mode = "list", length = 3)
g_8_1_bc_lst <- vector(mode = "list", length = 3)
g_8_6_bc_lst <- vector(mode = "list", length = 3)

for (t in times){
  for (i in 1:length(lines_g2)){
    l <- lines_g2[i]
    file_name <- paste0("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_2/common_barcodes/g_2_",
                        t,
                        "_common_CD27_CD54.best.",
                        l)
    print(file_name)
    bc <- read_csv(file_name, col_names = FALSE)
    bc <- unlist(bc)
    bc <- str_sub(bc, end = -3L)
    if (t %in% "0"){
      g_2_0_bc_lst[[i]] <- bc
    } else if (t %in% "1"){
      g_2_1_bc_lst[[i]] <- bc
    } else {
      g_2_6_bc_lst[[i]] <- bc
    }
  }
}

for (t in times){
  for (i in 1:length(lines_g8)){
    l <- lines_g8[i]
    file_name <- paste0("/nvmefs/scARC_Duan_018/hg38_Rnor6_mixed/group_8/common_barcodes/g_8_",
                        t,
                        "_common_CD08_CD25_CD26.best.",
                        l)
    print(file_name)
    bc <- read_csv(file_name, col_names = FALSE)
    bc <- unlist(bc)
    bc <- str_sub(bc, end = -3L)
    if (t %in% "0"){
      g_8_0_bc_lst[[i]] <- bc
    } else if (t %in% "1"){
      g_8_1_bc_lst[[i]] <- bc
    } else {
      g_8_6_bc_lst[[i]] <- bc
    }
  }
}
# assign cell line identities to the object
aggr_signac_unfiltered$cell.line.ident <- NA
aggr_signac_unfiltered$cell.line.ident[str_sub(aggr_signac_unfiltered@assays$peaks@counts@Dimnames[[2]], end = -3L) %in% 
                            c(g_2_0_bc_lst[[1]], g_2_1_bc_lst[[1]], g_2_6_bc_lst[[1]])] <- "CD_27"
aggr_signac_unfiltered$cell.line.ident[str_sub(aggr_signac_unfiltered@assays$peaks@counts@Dimnames[[2]], end = -3L) %in% 
                            c(g_2_0_bc_lst[[2]], g_2_1_bc_lst[[2]], g_2_6_bc_lst[[2]])] <- "CD_54"

aggr_signac_unfiltered$cell.line.ident[str_sub(aggr_signac_unfiltered@assays$peaks@counts@Dimnames[[2]], end = -3L) %in% 
                            c(g_8_0_bc_lst[[1]], g_8_1_bc_lst[[1]], g_8_6_bc_lst[[1]])] <- "CD_08"
aggr_signac_unfiltered$cell.line.ident[str_sub(aggr_signac_unfiltered@assays$peaks@counts@Dimnames[[2]], end = -3L) %in% 
                            c(g_8_0_bc_lst[[2]], g_8_1_bc_lst[[2]], g_8_6_bc_lst[[2]])] <- "CD_25"
aggr_signac_unfiltered$cell.line.ident[str_sub(aggr_signac_unfiltered@assays$peaks@counts@Dimnames[[2]], end = -3L) %in% 
                            c(g_8_0_bc_lst[[3]], g_8_1_bc_lst[[3]], g_8_6_bc_lst[[3]])] <- "CD_26"

# plot peak tracks for only the big cell types
aggr_glut <- subset(aggr_signac_unfiltered, subset = broad.cell.type == "glut")
aggr_GABA <- subset(aggr_signac_unfiltered, subset = broad.cell.type == "GABA")
aggr_NpCm_glut <- subset(aggr_signac_unfiltered, subset = fine.cell.type == "NEFM+/CUX2- glut")
aggr_NmCo_glut <- subset(aggr_signac_unfiltered, subset = fine.cell.type == "NEFM-/CUX2+ glut")

CoveragePlot(
  object = aggr_glut,
  region = "chr11-27689113-27774038", # BDNF
  region.highlight = GRanges(seqnames = "chr11",
                             ranges = IRanges(start = 27767000, end = 27777000)),
  group.by = "time.ident"
)

CoveragePlot(
  object = aggr_glut,
  region = "chr11-66416900-66423100", # NPAS4
  region.highlight = GRanges(seqnames = "chr11",
                             ranges = IRanges(start = 66418000, end = 66418700)),
  group.by = "time.ident"
)


# first look at response markers for glut: NPAS4 (early) and BDNF (late)
CD27 <- subset(x = aggr_glut, subset = cell.line.ident == "CD_27")
CD54 <- subset(x = aggr_glut, subset = cell.line.ident == "CD_54")
CD25 <- subset(x = aggr_glut, subset = cell.line.ident == "CD_25")
CD26 <- subset(x = aggr_glut, subset = cell.line.ident == "CD_26")
CD08 <- subset(x = aggr_glut, subset = cell.line.ident == "CD_08")

objs_by_line <- c(CD27, CD54, CD08, CD25, CD26)
obj_names <- c("CD27", "CD54", "CD08", "CD25", "CD26")
rm(file_name)

for (i in 1:length(objs_by_line)){
  obj <- objs_by_line[[i]]
  obj_name <- obj_names[i]
  file_name <- paste0(obj_name, "_glut_BDNF_peak_tracks.pdf")
  
  pdf(file = file_name)
  p <- CoveragePlot(
    object = obj,
    region = "chr11-27689113-27774038", # BDNF
    group.by = "time.ident")
  print(p)
  dev.off()
}

rm(file_name)
for (i in 1:length(objs_by_line)){
  obj <- objs_by_line[[i]]
  obj_name <- obj_names[i]
  file_name <- paste0(obj_name, "_glut_NPAS4_peak_tracks.pdf")
  
  pdf(file = file_name)
  p <- CoveragePlot(
    object = obj,
    region = "chr11-66416900-66423100", # NPAS4
    group.by = "time.ident"
  ) + ggtitle(obj_name)
  print(p)
  dev.off()
}

# now process GABA 
CoveragePlot(
  object = aggr_GABA,
  region = "chr8-28309000-28326519", # PNOC
  region.highlight = GRanges(seqnames = "chr11",
                            ranges = IRanges(start = 28312000, end = 28313500)),
  group.by = "time.ident"
)

CoveragePlot(
  object = aggr_GABA,
  region = "chr11-66416900-66423100", # NPAS4
  region.highlight = GRanges(seqnames = "chr11",
                             ranges = IRanges(start = 66418000, end = 66418700)),
  group.by = "time.ident"
)

CoveragePlot(
  object = aggr_GABA,
  region = "chr11-27689113-27774038", # BDNF
  region.highlight = GRanges(seqnames = "chr11",
                             ranges = IRanges(start = 27767000, end = 27777000)),
  group.by = "time.ident"
)

CD27 <- subset(x = aggr_GABA, subset = cell.line.ident == "CD_27")
CD54 <- subset(x = aggr_GABA, subset = cell.line.ident == "CD_54")
CD25 <- subset(x = aggr_GABA, subset = cell.line.ident == "CD_25")
CD26 <- subset(x = aggr_GABA, subset = cell.line.ident == "CD_26")
CD08 <- subset(x = aggr_GABA, subset = cell.line.ident == "CD_08")

objs_by_line <- c(CD27, CD54, CD08, CD25, CD26)
obj_names <- c("CD27", "CD54", "CD08", "CD25", "CD26")
rm(file_name)

for (i in 1:length(objs_by_line)){
  obj <- objs_by_line[[i]]
  obj_name <- obj_names[i]
  file_name <- paste0(obj_name, "_GABA_PNOC_peak_tracks.pdf")
  
  pdf(file = file_name)
  p <- CoveragePlot(
    object = obj,
    region.highlight = GRanges(seqnames = "chr11",
                               ranges = IRanges(start = 28312000, end = 28313500)),
    region = "chr8-28309000-28326519", # PNOC
    group.by = "time.ident")
  print(p)
  dev.off()
}

rm(file_name)
for (i in 1:length(objs_by_line)){
  obj <- objs_by_line[[i]]
  obj_name <- obj_names[i]
  file_name <- paste0(obj_name, "_GABA_NPAS4_peak_tracks.pdf")
  
  pdf(file = file_name)
  p <- CoveragePlot(
    object = obj,
    region = "chr11-66416900-66423100", # NPAS4
    region.highlight = GRanges(seqnames = "chr11",
                               ranges = IRanges(start = 66418000, end = 66418700)),
    group.by = "time.ident"
  ) + ggtitle(obj_name)
  print(p)
  dev.off()
}

for (i in 1:length(objs_by_line)){
  obj <- objs_by_line[[i]]
  obj_name <- obj_names[i]
  file_name <- paste0(obj_name, "_GABA_BDNF_peak_tracks.pdf")
  
  pdf(file = file_name)
  p <- CoveragePlot(
    object = obj,
    region = "chr11-27689113-27774038", # BDNF
    region.highlight = GRanges(seqnames = "chr11",
                               ranges = IRanges(start = 27767000, end = 27777000)),
    group.by = "time.ident")
  print(p)
  dev.off()
}

# plot two glut groups 
CoveragePlot(
  object = aggr_NmCo_glut,
  region = "chr11-27689113-27774038", # BDNF
  region.highlight = GRanges(seqnames = "chr11",
                             ranges = IRanges(start = 27767000, end = 27777000)),
  group.by = "time.ident"
) + ggtitle("NEFM-/CUX2+ glut")

CoveragePlot(
  object = aggr_NpCm_glut,
  region = "chr11-27689113-27774038", # BDNF
  region.highlight = GRanges(seqnames = "chr11",
                             ranges = IRanges(start = 27777000, end = 27772300)),
  group.by = "time.ident"
) + ggtitle("NEFM+/CUX2- glut")

# zoom in onto the enhancer region
CoveragePlot(
  object = aggr_NmCo_glut,
  region = "chr11-27770250-27771400", # BDNF
  group.by = "time.ident"
) + ggtitle("NEFM-/CUX2+ glut")
