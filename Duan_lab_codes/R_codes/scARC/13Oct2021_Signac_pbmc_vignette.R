# Chuxuan Li 13Oct2021
# Signac PBMC vignette
# Goal is to look at the pbmc object and confirm which attribute name is which

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)

set.seed(1234)

# set threads and parallelization
# plan("multisession", workers = 1)
# # plan("sequential")
# # plan()
# options(expressions = 20000)
# options(future.globals.maxSize = 21474836480)

# read in count matrix 
pbmc_counts <- Read10X_h5(filename = "/home/zhangs3/NVME/scARC_Duan_018/R_scARC_Duan_018/Signac/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
pbmc_metadata <- read.csv(
  file = "/home/zhangs3/NVME/scARC_Duan_018/R_scARC_Duan_018/Signac/atac_v1_pbmc_10k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

pbmc_chrom_assay <- CreateChromatinAssay(
  counts = pbmc_counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = '/home/zhangs3/NVME/scARC_Duan_018/R_scARC_Duan_018/Signac/atac_v1_pbmc_10k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = pbmc_chrom_assay,
  assay = "peaks",
  meta.data = pbmc_metadata
)
pbmc
pbmc[['peaks']]
granges(pbmc)


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# genome(annotations) <- "hg19"

# change to UCSC style since the data was mapped to hg19
annotations@seqnames@values <- paste0('chr', annotations@seqnames@values)
annotations@seqnames@values[annotations@seqnames@values == "chrMT"] <- "chrM"

annotations@seqinfo@seqnames <- paste0('chr', annotations@seqinfo@seqnames)
annotations@seqinfo@seqnames[annotations@seqinfo@seqnames == "chrMT"] <- "chrM"

annotations@seqinfo@genome <- rep_len("hg19", 25)

# add the gene information to the object
Annotation(pbmc) <- annotations

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
pbmc

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()

gene.activities <- GeneActivity(pbmc)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("pbmc_10k_v3.rds")

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2

pbmc <- subset(pbmc, idents = 14, invert = TRUE)
pbmc <- RenameIdents(
  object = pbmc,
  '0' = 'CD14 Mono',
  '1' = 'CD4 Memory',
  '2' = 'CD8 Effector',
  '3' = 'CD4 Naive',
  '4' = 'CD14 Mono',
  '5' = 'DN T',
  '6' = 'CD8 Naive',
  '7' = 'NK CD56Dim',
  '8' = 'pre-B',
  '9' = 'CD16 Mono',
  '10' = 'pro-B',
  '11' = 'DC',
  '12' = 'NK CD56bright',
  '13' = 'pDC'
)



# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'peaks'

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "NK CD56Dim",
  ident.2 = "CD14 Mono",
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)


plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14 Mono")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2


hist(da_peaks$avg_log2FC, breaks = 1000)
max(da_peaks$avg_log2FC)
#CD4 Naive: 0.9; CD4 memory: 0.77; CD8 Effector: 0.94; DN T: 0.85
# CD4 memory vs pre-B: 1.33
# NK CD56Dim vs DC: 1.13
# NK CD56Dim vs CD14 mono: 1.15

df <- da_peaks[, c("p_val", "avg_log2FC", "p_val_adj")]
df$significance <- "nonsig"
df$significance[df$p_val_adj < 0.05 & df$avg_log2FC > 0] <- "pos"
df$significance[df$p_val_adj < 0.05 & df$avg_log2FC < 0] <- "neg"
df$neg_log_p_val <- (0 - log10(df$p_val))

colors <- setNames(c("steelblue3", "grey", "red3"), c("neg", "nonsig", "pos"))

ggplot(data = df,
       aes(x = avg_log2FC, 
           y = neg_log_p_val, 
           color = significance)) + 
  geom_point() +
  #xlim(-0.5, 0.5) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  ggtitle("CD4 Naive vs. CD14 Mono \n (red points are CD4 Naive specific peaks)")
