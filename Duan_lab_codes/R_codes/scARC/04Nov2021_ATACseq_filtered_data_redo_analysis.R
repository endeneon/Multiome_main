# Chuxuan Li 11/4/2021
# 1. data: started from 10x-aggregated filtered matrix
# 2. preprocessing: removed the peaks that are mapped onto alternative transcripts
# 3. analyses in this script: TFIDF, SVD, UMAP, clustering, mapped onto RNAseq, response marker gene
#peak tracks are examined by cell-type, link gene to peaks, find differentially accessible
#peaks
# 4. output (already saved): "chrname_subsetted_10x_aggr_labeled.RData",
#"chrname_subsetted_10x_aggr_raw.RData", "RNAseq_ATACseq_matched_files.RData"
# 5. output (not saved): dimplots mapping ATACseq on RNAseq and vice versa, by cell type coverage plots for BDNF and NPAS4,
#by cell type and cell line coverage plots for BDNF and NPAS4, 
# 6. notes: this script contains the analyses done in "29Oct2021_diff_chromatin_accessibility.R",
#"11Oct2021_ATACseq_clustering.R", "27Oct2021_ATACseq_peaks_plot.R"; all three used the
#unfiltered 10x-aggregated matrix, which included chromosome names unaccepted by the steps
#used to link peaks to genes.

# init
library(Seurat)
library(sctransform)
library(glmGamPoi)

library(readr)
library(stringr)
library(dplyr)

library(ggplot2)
library(RColorBrewer)
library(viridis)
library(graphics)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)

library(future)

set.seed(1886)
# set threads and parallelization
plan("multisession", workers = 3)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 21474836480)

# read in the combined object obtained by the -aggr method from 10x
aggr_raw <- Read10X_h5(filename = "/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/filtered_feature_bc_matrix.h5")

aggr <- CreateSeuratObject(
  counts = aggr_raw$`Gene Expression`,
  assay = "RNA"
)

peaks <- aggr_raw$Peaks

# choose the peaks that only start with chr
peaks_subsetted <- peaks[str_sub(peaks@Dimnames[[1]], end = 3L) %in% "chr", ]

# create chromatin assay without filtering
chrom_assay <- CreateChromatinAssay(
  counts = peaks_subsetted,
  sep = c(":", "-"),
  #genome = 'hg38',
  fragments = '/home/cli/NVME/scARC_Duan_018/GRCh38_aggregation/sample_2_8_aggr_11Oct2021/outs/atac_fragments.tsv.gz',
  annotation = annot_aggr_signac_ucsc
  # min.cells = 10,
  # min.features = 200
)

aggr[["ATAC"]] <- chrom_assay

#save("peaks_subsetted", file = "chrname_subsetted_10x_aggr_raw.RData")

# pks_raw <-
#   CreateSeuratObject(counts = peaks_subsetted,
#                      project = "aggregated_raw_pks")

aggr$time.group.ident <- str_sub(aggr@assays$ATAC@counts@Dimnames[[2]], 
                                       start = -1L)

aggr$time.ident <- NA
aggr$time.ident[aggr$time.group.ident %in% c("1", "4")] <- "0hr"
aggr$time.ident[aggr$time.group.ident %in% c("2", "5")] <- "1hr"  
aggr$time.ident[aggr$time.group.ident %in% c("3", "6")] <- "6hr"  
# check if anything is unlabeled
sum(is.na(aggr$time.ident))

# annot_aggr <-
#   GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86,
#                       standard.chromosomes = T)
# make a txdb from Gencode database
# txdb.v32 <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# 
# exon_txdb_v32 <- transcripts(txdb.v32)


#seqlevelsStyle(annot_aggr_signac) <- 'UCSC' # cannot  use FTP service


# QC #####
DefaultAssay(aggr) <- "ATAC"
aggr <- NucleosomeSignal(aggr)
aggr <- TSSEnrichment(aggr, 
                      fast = F # if True, cannot plot the accessibility profile at the TSS
)

# using fewer filters filtered object to compare with filtered results
aggr_filtered <- subset(x = aggr,
                        subset = nFeature_ATAC < 20000 &
                          #pct_reads_in_peaks > 15 &
                          #blacklist_ratio < 0.05 &
                          nucleosome_signal < 4 &
                          TSS.enrichment > 2)

#save("aggr_filtered", file = "chrname_subsetted_10x_aggr_filtered.RData")

# dimensional reduction ####
aggr_filtered <- RunTFIDF(aggr_filtered)
aggr_filtered <- FindTopFeatures(aggr_filtered, min.cutoff = 'q0')
aggr_filtered <- RunSVD(aggr_filtered)

save("aggr_filtered", file = "ATAC_RNA_combined_10x_aggregated_data_filtered.RData")

DefaultAssay(aggr_filtered) <- "RNA"
aggr_filtered <- PercentageFeatureSet(aggr_filtered, pattern = "^MT-", col.name = "percent.mt")


aggr_filtered <- SCTransform(aggr_filtered,
                    vars.to.regress = "percent.mt",
                    method = "glmGamPoi",
                    variable.features.n = 9000,
                    verbose = T,
                    seed.use = 42)
aggr_filtered <- RunPCA(aggr_filtered)

# build a joint neighbor graph using both assays
aggr_filtered <- FindMultiModalNeighbors(
  object = aggr_filtered,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
aggr_filtered <- RunUMAP(
  object = aggr_filtered,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)


aggr_filtered <- FindNeighbors(object = aggr_filtered, 
                                reduction = 'lsi', 
                                dims = 2:30) # remove the first component
aggr_filtered <- FindClusters(object = aggr_filtered, 
                               verbose = FALSE, 
                              graph.name = "ATAC_snn",
                               algorithm = 3, 
                               resolution = 0.5,
                               random.seed = 99)

# integrate with RNAseq dataset ####
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/codes/integrated_with_all_labels.RData")
# get the barcodes for each cluster from the RNAseq object
DefaultAssay(aggr_filtered) <- "ATAC"
aggr_filtered$RNA.cluster.ident <- NA
aggr_filtered$trimmed.barcodes <- str_sub(aggr_filtered@assays$ATAC@data@Dimnames[[2]],
                                           end = -3L)

for (i in 1:length(levels(integrated_renamed_1$seurat_clusters))){
  barcodes <- str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]][integrated_renamed_1$seurat_clusters %in% 
                                                                             as.character(i)], end = -7L)
  
  aggr_filtered$RNA.cluster.ident[aggr_filtered$trimmed.barcodes 
                                   %in% barcodes] <- as.character(i) 
}

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = aggr_filtered, 
        label = TRUE, 
        repel = T,
        reduction = "umap",
        group.by = "RNA.cluster.ident") + 
  ggtitle("clustering of ATACseq data with RNAseq projection")


# project ATACseq onto RNAseq
integrated_renamed_1$ATAC.cluster.ident <- NA
integrated_renamed_1$trimmed.barcodes <- str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]], end = -7L)

for (i in 1:length(levels(aggr_filtered$seurat_clusters))){
  barcodes <- str_sub(aggr_filtered@assays$ATAC@counts@Dimnames[[2]][aggr_filtered$seurat_clusters %in% 
                                                                         as.character(i)], end = -3L)
  
  integrated_renamed_1$ATAC.cluster.ident[integrated_renamed_1$trimmed.barcodes 
                                          %in% barcodes] <- as.character(i) 
}

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = integrated_renamed_1, 
        label = TRUE, 
        repel = T,
        group.by = "ATAC.cluster.ident") + 
  ggtitle("clustering of RNAseq data with ATACseq projection")

unique(aggr_filtered$RNA.cell.type.ident)
# group the subtypes into broad categories
aggr_filtered$broad.cell.type <- NA
aggr_filtered$broad.cell.type[aggr_filtered$RNA.cell.type.ident %in% 
                                c("GAD1+/GAD2+ GABA", "GAD1+ GABA (less mature)", "GAD1+ GABA")] <- "GABA"
aggr_filtered$broad.cell.type[aggr_filtered$RNA.cell.type.ident %in% 
                                c("NEFM+/CUX2- glut", "NEFM-/CUX2+ glut",
                                  "NEFM+/CUX2-, SST+ glut", "NEFM-/CUX2+, ADCYAP1+ glut",
                                  "NEFM+/CUX2-, ADCYAP1+ glut", "NEFM-/CUX2+, TCERG1L+ glut",
                                  "NEFM+/CUX2+, SST+ glut")] <- "glut"
aggr_filtered$broad.cell.type[aggr_filtered$RNA.cell.type.ident %in% 
                                c("NPC", "immature neuron", "MAP2+ NPC", "subcerebral immature neuron")] <- "NPC"
# divide the glut
aggr_filtered$fine.cell.type <- NA
aggr_filtered$fine.cell.type[aggr_filtered$RNA.cell.type.ident %in% 
                               c("GAD1+/GAD2+ GABA", "GAD1+ GABA (less mature)", "GAD1+ GABA")] <- "GABA"
aggr_filtered$fine.cell.type[aggr_filtered$RNA.cell.type.ident %in% 
                               c("NEFM+/CUX2- glut", "NEFM+/CUX2-, SST+ glut", 
                                 "NEFM+/CUX2-, ADCYAP1+ glut")] <- "NEFM+/CUX2- glut"
aggr_filtered$fine.cell.type[aggr_filtered$RNA.cell.type.ident %in% 
                               c("NEFM-/CUX2+, glut",
                                 "NEFM-/CUX2+, ADCYAP1+ glut",
                                 "NEFM-/CUX2+, TCERG1L+ glut"
                               )] <- "NEFM-/CUX2+ glut"

aggr_filtered$fine.cell.type[aggr_filtered$RNA.cell.type.ident %in% "SEMA3E+ glut"] <- "SEMA3E+ glut"

aggr_filtered$fine.cell.type[aggr_filtered$RNA.cell.type.ident %in% c("MAP2+ NPC", "NPC")] <- "NPC"

aggr_filtered$fine.cell.type[aggr_filtered$RNA.cell.type.ident %in% c("immature forebrain glut", 
                                                                      "forebrain NPC")] <- "forebrain"
aggr_filtered$fine.cell.type[aggr_filtered$RNA.cell.type.ident == "subcerebral immature neuron"] <- "subcerebral immature neuron"
aggr_filtered$fine.cell.type[aggr_filtered$RNA.cell.type.ident == "immature neuron"] <- "immature neuron"


# link peaks to genes ####

DefaultAssay(aggr_filtered) <- "ATAC"

# compute the GC content for each peak
aggr_filtered <- RegionStats(aggr_filtered, genome = BSgenome.Hsapiens.UCSC.hg38)

# separate time points
times <- c("0hr", "1hr", "6hr")
aggr_subsetted_by_time <- vector(mode = "list", length = 3L)
aggr_by_time_linked <- vector(mode = "list", length = 3L)
goi <- c("BDNF", "NPAS4",
         "VPS45","OTUD7B",
         "BCL11B", "SETD1A")

for (i in 1:length(times)){
  t <- times[i]
  print(t)
  aggr_subsetted_by_time[[i]] <- subset(aggr_filtered, subset = time.ident == t)
  # link peaks to genes
  aggr_by_time_linked[[i]] <- LinkPeaks(
    object = aggr_subsetted_by_time[[i]],
    peak.assay = "ATAC",
    expression.assay = "SCT",
    genes.use = goi
  )
}

# ranges_list <- c("chr11-27689113-27774038", #BDNF
#                  "chr11-66414779-66433926", #NPAS4
#                  "chr1-150028428-150184241", #VPS45
#                  "chr1-149919583-150028955", #OTUD7B
#                  "chr14-99143559-99297925", #BCL11B
#                  "chr16-30933106-31013838" #SETD1A
#                  )
# visualize with CoveragePlot or CoverageBrowser
Idents(aggr_filtered) <- "fine.cell.type"

for (i in 1:length(goi)){
  g <- goi[i] 
  # grange <- ranges_list[i]
  print(g)
  for (j in 1:length(aggr_by_time_linked)){
    t <- times[j]
    file_name <- paste0(g, "_", t, "_", "link_peaks_to_gene.pdf")
    print(file_name)
    
    pdf(file = file_name)
    p <- CoveragePlot(
      object = aggr_by_time_linked[[j]],
      region = g,
      features = g,
      expression.assay = "SCT",
      idents = idents.plot,
      extend.upstream = 50000,
      extend.downstream = 50000
    )
    print(p)
    dev.off()
  }
}

# look at the interactive version
CoverageBrowser(object = aggr_by_time_linked[[2]], 
                region = "chr11-27689113-27774038", 
                idents = idents.plot,
                sep = c("-", "-"))


for (i in 1:length(times)){
  t <- times[i]
  print(t)
  aggr_subsetted_by_time[[i]] <- subset(aggr_filtered, subset = time.ident == t)
  # link peaks to genes
  aggr_by_time_linked[[i]] <- LinkPeaks(
    object = aggr_subsetted_by_time[[i]],
    peak.assay = "ATAC",
    expression.assay = "SCT"
  )
}

# the version of plotting that works:
# Idents(aggr_by_time_linked[[1]])
# p <- CoveragePlot(
#   object = aggr_by_time_linked[[1]],
#   region = ranges_list[1],
#   expression.assay = "SCT",
#   idents = idents.plot
# )
# print(p)


# add an attribute onto the RNAseq data from ATACseq nCount and nFeature
integrated_renamed_1$ATAC.nCount <- NA
integrated_renamed_1$ATAC.nFeature <- NA
integrated_renamed_1$ATAC.nCount[integrated_renamed_1$trimmed.barcodes %in% aggr_filtered$trimmed.barcodes] <- 
  aggr_filtered$nCount_peaks[integrated_renamed_1$trimmed.barcodes %in% aggr_filtered$trimmed.barcodes]
integrated_renamed_1$ATAC.nFeature[integrated_renamed_1$trimmed.barcodes %in% aggr_filtered$trimmed.barcodes] <- 
  aggr_filtered$nFeature_peaks[integrated_renamed_1$trimmed.barcodes %in% aggr_filtered$trimmed.barcodes]


# determine the cutoff
hist(x = integrated_renamed_1$ATAC.nCount,
     breaks = 1000,
     xlim = c(0, 50000))
ncount_feat <- FeaturePlot(integrated_renamed_1, 
                           features = "ATAC.nCount",
                           max.cutoff = 25000,
                           cols = c("grey", "red3")
)
ncount_feat[[1]]$layers[[1]]$aes_params$alpha = .3
ncount_feat

hist(x = integrated_renamed_1$ATAC.nFeature,
     breaks = 1000,
     xlim = c(0, 20000))
nfeature_feat <- FeaturePlot(integrated_renamed_1, 
                             features = "ATAC.nFeature",
                             max.cutoff = 16000,
                             cols = c("grey", "red3"))
nfeature_feat[[1]]$layers[[1]]$aes_params$alpha = .3
nfeature_feat


# check how many cells in RNAseq are not labeled by ATACseq
sum(is.na(integrated_renamed_1$ATAC.nFeature))/length(integrated_renamed_1$ATAC.nFeature)
# check how many cells in ATACseq are not labeled by RNAseq
sum(is.na(aggr_filtered$RNA.cluster.ident))/length(aggr_filtered$RNA.cluster.ident)

# without QC for ATACseq
sum(str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]], end = -7L) %in% 
      str_sub(aggr_signac@assays$peaks@data@Dimnames[[2]],end = -3L)) / 
  length(integrated_renamed_1@assays$RNA@data@Dimnames[[2]])

sum(str_sub(aggr_signac@assays$peaks@data@Dimnames[[2]], end = -3L) %in% 
      str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]], end = -7L))/ 
  length(aggr_signac@assays$peaks@data@Dimnames[[2]])


# check duplicate barcodes 
length(aggr_filtered$trimmed.barcodes) #total number of barcodes (ATAC): 80970
length(unique(aggr_filtered$trimmed.barcodes)) # unique num of barcodes (ATAC): 77183
length(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]]) #total num of barcodes (RNA): 83093
length(unique(integrated_renamed_1$trimmed.barcodes)) # unique num of barcodes (RNA): 79159
77183/80970 # unique ATAC barcodes to all ATAC barcodes ratio: 0.9532296
79159/83093 # unique RNA barcodes to all RNA barcodes ratio: 0.9526555

sum(is.na(aggr_filtered$RNA.cluster.ident))

# create another attribute to project RNAseq cell type ident onto the ATACseq data
aggr_filtered$RNA.cell.type.ident <- NA
# find all barcodes corresponding to each cell type

for (i in levels(integrated_renamed_1@active.ident)){
  barcodes <- str_sub(integrated_renamed_1@assays$RNA@counts@Dimnames[[2]][integrated_renamed_1@active.ident %in% i], 
                      end = -7L)
  
  aggr_filtered$RNA.cell.type.ident[aggr_filtered$trimmed.barcodes %in% barcodes] <- i
}

save("aggr_filtered", file = "ATAC_RNA_combined_10x_aggr_labeled.RData")

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = aggr_filtered, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.cell.type.ident") + 
  ggtitle("clustering of ATACseq data with cell type projection")

# peak tracks ####
# plot peaks for exc/inh cells and separated by time and cell line (like in the grant)



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
aggr_filtered$cell.line.ident <- NA
aggr_filtered$cell.line.ident[str_sub(aggr_filtered@assays$peaks@counts@Dimnames[[2]], end = -3L) %in% 
                                 c(g_2_0_bc_lst[[1]], g_2_1_bc_lst[[1]], g_2_6_bc_lst[[1]])] <- "CD_27"
aggr_filtered$cell.line.ident[str_sub(aggr_filtered@assays$peaks@counts@Dimnames[[2]], end = -3L) %in% 
                                 c(g_2_0_bc_lst[[2]], g_2_1_bc_lst[[2]], g_2_6_bc_lst[[2]])] <- "CD_54"

aggr_filtered$cell.line.ident[str_sub(aggr_filtered@assays$peaks@counts@Dimnames[[2]], end = -3L) %in% 
                                 c(g_8_0_bc_lst[[1]], g_8_1_bc_lst[[1]], g_8_6_bc_lst[[1]])] <- "CD_08"
aggr_filtered$cell.line.ident[str_sub(aggr_filtered@assays$peaks@counts@Dimnames[[2]], end = -3L) %in% 
                                 c(g_8_0_bc_lst[[2]], g_8_1_bc_lst[[2]], g_8_6_bc_lst[[2]])] <- "CD_25"
aggr_filtered$cell.line.ident[str_sub(aggr_filtered@assays$peaks@counts@Dimnames[[2]], end = -3L) %in% 
                                 c(g_8_0_bc_lst[[3]], g_8_1_bc_lst[[3]], g_8_6_bc_lst[[3]])] <- "CD_26"

# plot peak tracks for only the big cell types
aggr_glut <- subset(aggr_filtered, subset = broad.cell.type == "glut")
aggr_GABA <- subset(aggr_filtered, subset = broad.cell.type == "GABA")
aggr_NpCm_glut <- subset(aggr_filtered, subset = fine.cell.type == "NEFM+/CUX2- glut")
aggr_NmCo_glut <- subset(aggr_filtered, subset = fine.cell.type == "NEFM-/CUX2+ glut")

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
                             ranges = IRanges(start = 27767000, end = 27777000)),
  group.by = "time.ident"
) + ggtitle("NEFM+/CUX2- glut")

# zoom in onto the enhancer region
CoveragePlot(
  object = aggr_NmCo_glut,
  region = "chr11-27770250-27771400", # BDNF
  group.by = "time.ident"
) + ggtitle("NEFM-/CUX2+ glut")



# find differentially accessible peaks ####
# list all cell types
cell_types <- c("GABA", "NEFM+/CUX2- glut", "NEFM-/CUX2+ glut", "SEMA3E+ glut", "NPC", "forebrain", "immature neuron", 
                "subcerebral immature neuron")
SEMA3Ep_glut <- subset(aggr_filtered, subset = fine.cell.type %in% "SEMA3E+ glut")
NPC <- subset(aggr_filtered, subset = fine.cell.type %in% "NPC")
forebrain <- subset(aggr_filtered, subset = fine.cell.type %in% "forebrain")
subcerebral <- subset(aggr_filtered, subset = fine.cell.type %in% "subcerebral immature neuron")
immature_neuron <- subset(aggr_filtered, subset = fine.cell.type %in% "immature neuron")



