# Chuxuan Li 04/15/2022
# MACS2 call peaks on the 10x-aggred objects with RNA cell types already added

# init####
library(Seurat)
library(Signac)
library(sctransform)
library(glmGamPoi)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg38)
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

plan("multisession", workers = 2)
set.seed(2001)
options(future.globals.maxSize = 229496729600)

# use MACS2 to call new peaks ####
peaks <- CallPeaks(
  object = ATAC,
  group.by = "timextype.ident",
  macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2"
)
save(peaks, file = "macs2_called_new_peaks.RData")
peaks

peaks_uncombined <- CallPeaks(
  object = ATAC, 
  assay = "ATAC", 
  group.by = "timextype.ident",
  macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2",
  combine.peaks = F
)
save(peaks_uncombined, file = "macs2_called_new_peak_by_celltype_separated.RData")

# save peak sets to bed file ####
df <- data.frame(seqnames = seqnames(peaks),
                   starts = start(peaks) - 1,
                   ends = end(peaks),
                   strands = strand(peaks))
df <- df[grep("chr", df$seqnames), ]
write.table(df, 
            file = "./MACS2_called_peak_bed_files/18line_union_peak_set.bed", 
            quote = F, 
            sep = "\t", 
            row.names = F, col.names = F)

for (i in 1:length(peaks_uncombined)){
  df <- data.frame(seqnames = seqnames(peaks_uncombined[[i]]),
                   starts = start(peaks_uncombined[[i]]) - 1,
                   ends = end(peaks_uncombined[[i]]),
                   name = peaks_uncombined[[i]]@elementMetadata@listData$name,
                   score = peaks_uncombined[[i]]@elementMetadata@listData$neg_log10qvalue_summit,
                   strands = strand(peaks_uncombined[[i]]))
  df <- df[grep("chr", df$seqnames), ]
  file_name <- paste0("./MACS2_called_peak_bed_files/peaks_called_in_",
                      unique(peaks_uncombined[[i]]$ident),
                      ".bed")
  write.table(df, 
              file = file_name, 
              quote = F, 
              sep = "\t", 
              row.names = F, col.names = F)
  
}

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak ####
macs2_counts <- FeatureMatrix(
  fragments = Fragments(ATAC),
  features = peaks,
  cells = colnames(ATAC)
)
fragpath <- "//data/FASTQ/Duan_Project_024/hybrid_output/hybrid_aggr_5groups_no2ndary/outs/atac_fragments.tsv.gz"

#load("../R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")
# create a new assay using the MACS2 peak set and add it to the Seurat object ####
ATAC[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath#,
  #annotation = ens_use
)

save(ATAC, file = "ATAC_obj_added_MACS2_peaks.RData")

DefaultAssay(ATAC_new) <- "peaks"

# calculate percent fragments in peaks 
total_fragments <- CountFragments(fragpath)
barcodes <- ATAC@assays$ATAC@counts@Dimnames[[2]]
ATAC$fragments <- total_fragments$frequency_count[total_fragments$CB %in% barcodes]

ATAC <- FRiP(
  object = ATAC,
  assay = 'peaks',
  total.fragments = 'fragments'
)
hist(ATAC$FRiP, breaks = 10000, xlim = c(0, 100))

# dimensional reduction and clustering of new peaks ####
ATAC_new <- FindTopFeatures(ATAC, min.cutoff = 5)
ATAC_new <- RunTFIDF(ATAC_new)
ATAC_new <- RunSVD(ATAC_new)

# check which components are associated with technical (nonbiological) variance (sequencing depth)
DepthCor(ATAC_new)

# clustering and downstream analysis
ATAC_new <- RunUMAP(object = ATAC_new, 
                    reduction = 'lsi', 
                    dims = 2:30, # remove the first component
                    seed.use = 999) 

ATAC_new <- FindNeighbors(object = ATAC_new, 
                          reduction = 'lsi', 
                          dims = 2:30) # remove the first component

ATAC_new <- FindClusters(object = ATAC_new, 
                         verbose = FALSE, 
                         algorithm = 3, 
                         resolution = 0.5,
                         random.seed = 99)

# plot dimplots
DimPlot(object = ATAC_new, label = F) +
  ggtitle("ATACseq data with new, MACS2-called peak set")


# project new ATACseq onto RNAseq and vice versa ####
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")
# apply RNAseq cell types ####
types <- as.vector(unique(integrated_labeled$cell.type))
ATAC_new$RNA.cell.type <- "unknown"
for (i in 1:length(types)){
  print(types[i])
  subRNA <- subset(integrated_labeled, cell.type == types[i])
  ATAC_new$RNA.cell.type[str_extract(ATAC_new@assays$ATAC@counts@Dimnames[[2]], "^[A-Z]+-") %in% 
                           str_extract(subRNA@assays$RNA@counts@Dimnames[[2]], "^[A-Z]+-")] <- types[i]
}
unique(ATAC_new$RNA.cell.type)
sum(ATAC_new$RNA.cell.type == "unknown") #8899
DimPlot(object = ATAC_new, label = TRUE, group.by = "RNA.cell.type") + 
  NoLegend() +
  ggtitle("ATACseq data projected by RNAseq cell types")


integrated_labeled$ATAC.cluster.ident <- NA
clusters <- unique(ATAC_new$seurat_clusters)
for (i in 1:length(clusters)){
  print(clusters[i])
  subATAC <- subset(ATAC_new, seurat_clusters == clusters[i])
  integrated_labeled$ATAC.cluster.ident[str_extract(integrated_labeled@assays$RNA@counts@Dimnames[[2]], "^[A-Z]+-") %in% 
                                          str_extract(subATAC@assays$ATAC@counts@Dimnames[[2]], "^[A-Z]+-")] <-
    clusters[i]
}

ATAC_new$RNA.cluster.ident <- NA
clusters <- as.vector(unique(integrated_labeled$seurat_clusters))
for (i in 1:length(clusters)){
  print(clusters[i])
  subRNA <- subset(integrated_labeled, seurat_clusters == clusters[i])
  ATAC_new$RNA.cluster.ident[str_extract(ATAC_new@assays$ATAC@counts@Dimnames[[2]], "^[A-Z]+-") %in% 
                               str_extract(subRNA@assays$RNA@counts@Dimnames[[2]], "^[A-Z]+-")] <- clusters[i]
}

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = integrated_labeled, 
        label = TRUE, 
        repel = T,
        group.by = "ATAC.cluster.ident") + 
  ggtitle("RNAseq data \n projected by ATACseq cluster numbers")
DimPlot(object = integrated_labeled, 
        label = TRUE, 
        repel = T)

DimPlot(object = ATAC_new, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.cluster.ident") + 
  ggtitle("ATACseq data \n projected by RNAseq cluster numbers") +
  NoLegend()

DimPlot(object = ATAC_new, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.cell.type") + 
  ggtitle("ATACseq data \n projected by RNAseq cell types") +
  NoLegend()

save(ATAC_new, file = "ATAC_new_obj_with_MACS2_new_peaks_and_RNA_labels.RData")

