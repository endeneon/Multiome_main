# Chuxuan Li 04/01/2022
# find co-accessible peaks on 5-line ATACseq data with Cicero

# init ####
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(cicero)
library(monocle3)

library(ggplot2)
library(patchwork)
library(future)

plan("multisession", workers = 1)

# load data ####
load("./ATACseq_5line_objonly_labeled.RData")
unique(obj_complete$fine.cell.type)
nmglut <- subset(obj_complete, fine.cell.type == "NEFM-/CUX2+ glut")
npglut <- subset(obj_complete, fine.cell.type == "NEFM+/CUX2- glut")
GABA <- subset(obj_complete, fine.cell.type == "GABA")

# create Cicero object ####
nmglut.cds <- SeuratWrappers::as.cell_data_set(x = nmglut)
nmglut.cicero <- make_cicero_cds(nmglut.cds, 
                                 reduced_coordinates = reducedDims(nmglut.cds)$UMAP)

# get the chromosome sizes from the Seurat object
genome <- 249250621

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = "chr1", "length" = genome)

# run cicero
conns <- run_cicero(nmglut.cicero, genomic_coords = genome.df, sample_num = 100)
head(conns)
ccans <- generate_ccans(conns, coaccess_cutoff_override = )
# add links to Seurat object
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(nmglut) <- links
# check VPS45 reagion
CoveragePlot(bone, region = "chr1-40189344-40252549")
