# Chuxuan Li 02/07/2023
# Use MACS2 to recall the peaks from 029, then re-count the fragments, and repeat
#normalizing, clustering, plotting

# init ####
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

load("./029_multiomic_obj_normalized_and_clustered_w_RNA_celltypes.RData")

# use MACS2 to call new peaks ####
peaks <- CallPeaks(
  object = multiomic_obj, assay = "ATAC",  
  group.by = "timextype.ident",
  macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2"
)
save(peaks, file = "029_macs2_called_new_peaks_combined.RData")

overlap_boo <- c()
for(i in 1:length(peaks)) {
  overlap_boo <- c(overlapsAny(peaks[i], peaks[i+1]))
}
peaks_uncombined <- CallPeaks(
  object = multiomic_obj, assay = "ATAC", 
  group.by = "timextype.ident",
  macs2.path = "~/Data/cli_anaconda3/cli_conda_env/rstudio/bin/macs2",
  combine.peaks = F
)
save(peaks_uncombined, file = "029_macs2_called_new_peaks_sep_by_timextype.RData")

# save peak sets to bed file ####
df <- data.frame(seqnames = seqnames(peaks),
                 starts = start(peaks) - 1,
                 ends = end(peaks),
                 strands = strand(peaks))
df <- df[grep("chr[0-9]+", df$seqnames), ]
write.table(df, 
            file = "./MACS2_called_new_peaks/029_union_peak_set.bed", 
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
  file_name <- paste0("./MACS2_called_new_peaks/029_sep_peaks_called_in_",
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
unique(peaks@seqnames)
peaks <- peaks[str_detect(peaks@seqnames, "chr[0-9]+"), ] # remove chrX and chrY

# quantify counts in each peak ####
macs2_counts <- FeatureMatrix(
  fragments = Fragments(multiomic_obj),
  features = peaks,
  cells = colnames(multiomic_obj)
)
save(macs2_counts, file = "MACS2_recalled_029_peaks_FeatureMatrix_output.RData")

fragpath <- "/data/FASTQ/Duan_Project_029/GRCh38_only/Duan_Project_029_aggr/outs/atac_fragments.tsv.gz"

load("../R_scARC_Duan_018/codes/EnsDb_UCSC_hg38_annotation.RData")
# create a new assay using the MACS2 peak set and add it to the Seurat object ####
multiomic_obj[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = ens_use
)
save(multiomic_obj, file = "multiomic_obj_added_MACS2_new_peaks.RData")

# calculate percent fragments in peaks ####
DefaultAssay(multiomic_obj) <- "peaks"

total_fragments <- CountFragments(fragpath)
bc <- multiomic_obj@assays$peaks@counts@Dimnames[[2]]
multiomic_obj$fragment_count <- total_fragments$reads_count[total_fragments$CB %in% bc]
multiomic_obj$fragment_count[1:10]
colSums(multiomic_obj@assays$peaks@counts)[1:10]
frip <- colSums(multiomic_obj@assays$peaks@counts)/multiomic_obj$fragment_count
hist(colSums(multiomic_obj@assays$peaks@counts)[frip < 1])
hist(ATAC_clean$fragment_count[frip < 1], xlim = c(0, 2e4), breaks = 1000)
hist(ATAC_clean$fragment_count[frip > 1], xlim = c(0, 2e3), breaks = 1000)
hist(multiomic_obj$fragment_count, xlim = c(0, 2e4), breaks = 1000)
hist(frip)
multiomic_obj <- FRiP(
  object = multiomic_obj,
  assay = 'peaks',
  total.fragments = 'fragment_count'
)
hist(frip, breaks = 10000, xlim = c(0, 1), 
     main = "Fraction of fragments in peaks - 025", xlab = "FRIP")
max(frip)

# dimensional reduction and clustering of new peaks ####
multiomic_obj <- FindTopFeatures(multiomic_obj, min.cutoff = 5)
multiomic_obj <- RunTFIDF(multiomic_obj)
multiomic_obj <- RunSVD(multiomic_obj)

# check which components are associated with technical (nonbiological) variance (sequencing depth)
DepthCor(multiomic_obj)

# clustering and downstream analysis
multiomic_obj_new <- RunUMAP(object = multiomic_obj, 
                    reduction = 'lsi', 
                    dims = 2:30, # remove the first component
                    seed.use = 999) 

multiomic_obj_new <- FindNeighbors(object = multiomic_obj_new, 
                          reduction = 'lsi', 
                          dims = 2:30) # remove the first component

multiomic_obj_new <- FindClusters(object = multiomic_obj_new, 
                         verbose = FALSE, 
                         algorithm = 3, 
                         resolution = 0.5,
                         random.seed = 99)

# plot dimplots
DimPlot(object = multiomic_obj_new, label = F) +
  ggtitle("ATACseq data clustering result with new, MACS2-called peak set") +
  theme(text = element_text(size = 10))


# project new ATACseq onto RNAseq and vice versa ####
load("../Duan_project_029_RNA/029_RNA_integrated_labeled.RData")
# apply RNAseq cell types ####
types <- as.vector(unique(integrated_labeled$cell.type))
unique(multiomic_obj_new$RNA.cell.type)
sum(multiomic_obj_new$RNA.cell.type == "unidentified") #16697
DimPlot(object = multiomic_obj_new, label = TRUE, group.by = "RNA.cell.type") + 
  NoLegend() +
  ggtitle("ATACseq data w/ new peak set projected by RNAseq cell types") +
  theme(text = element_text(size = 10))

# map RNA/ATAC clusters onto each other
integrated_labeled$ATAC.cluster.ident <- NA
clusters <- unique(multiomic_obj_new$seurat_clusters)
for (c in clusters){
  print(c)
  subATAC <- subset(multiomic_obj_new, seurat_clusters == c)
  integrated_labeled$ATAC.cluster.ident[str_extract(integrated_labeled@assays$RNA@counts@Dimnames[[2]], "^[A-Z]+-") %in% 
                                          str_extract(subATAC@assays$ATAC@counts@Dimnames[[2]], "^[A-Z]+-")] <- c
}
integrated_labeled$ATAC.cluster.ident <- as.numeric(integrated_labeled$ATAC.cluster.ident)
multiomic_obj_new$RNA.cluster.ident <- NA
clusters <- as.vector(unique(integrated_labeled$seurat_clusters))
for (c in clusters){
  print(c)
  subRNA <- subset(integrated_labeled, seurat_clusters == c)
  multiomic_obj_new$RNA.cluster.ident[str_extract(multiomic_obj_new@assays$ATAC@counts@Dimnames[[2]], "^[A-Z]+-") %in% 
                               str_extract(subRNA@assays$RNA@counts@Dimnames[[2]], "^[A-Z]+-")] <- c
}
multiomic_obj_new$RNA.cluster.ident <- as.numeric(multiomic_obj_new$RNA.cluster.ident)

# plot projections from RNAseq clusters onto the dimplot
DimPlot(object = integrated_labeled, 
        label = TRUE, 
        repel = T,
        group.by = "ATAC.cluster.ident") + 
  ggtitle("RNAseq data \n projected by ATACseq cluster numbers") +
  theme(text = element_text(size = 10))

DimPlot(object = multiomic_obj_new, 
        label = TRUE, 
        repel = T,
        group.by = "RNA.cluster.ident") + 
  ggtitle("ATACseq data \n projected by RNAseq cluster numbers") +
  theme(text = element_text(size = 10))

save(multiomic_obj_new, file = "multiomic_obj_with_new_peaks_labeled.RData")

