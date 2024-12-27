# Siwei 07 Apr 2022
# use 18 lines, Glut cells at 0 hr to check 
# rs2027349-c1orf54 connections

# need chr1 only

#init

library(Seurat)
library(SeuratWrappers)
library(Signac)
library(future)

library(remotes)
library(monocle3)
# library(monocle)
library(cicero)

library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(stringr)

library(ggplot2)
library(patchwork)

plan("multisession", workers = 4)
options(future.globals.maxSize = 100 * 1024^3)

# load data
# load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/18line_10xaggred_labeled_by_RNA_obj.RData")

DefaultAssay(ATAC)

ATAC_chr1 <-
  ATAC[rownames(ATAC) %in% "chr1-", ]

ATAC_chr1 <-
  ATAC[str_detect(string = rownames(ATAC),
                  pattern = "chr1\\-"), ]
rm(ATAC)

# keep glut only, need to compare
# all gluts
# 2 glut types
# all 0 hr gluts
# glut 0 hr types
ATAC_chr1_glut <-
  ATAC_chr1[, ATAC_chr1$cell.type %in% 
              c("NEFM_neg_glut", "NEFM_pos_glut")]
# ATAC_chr1_glut <-
#   ATAC_chr1[(ATAC_chr1$cell.type == "NEFM_neg_glut") |
#               (ATAC_chr1$cell.type == "NEFM_pos_glut"), ]
Idents(ATAC_chr1_glut) <- "cell.type"
CoveragePlot(ATAC_chr1_glut, 
             region = "chr1-150059035-150269039") 

# check the origin of peaks in the large peak list,
# first make a peakset from all glut cells to test
# take chr1 only
names(peaks_uncombined)
# [1] "GABA_0hr"          "GABA_6hr"          "NEFM_neg_glut_0hr"
# [4] "NEFM_neg_glut_1hr" "GABA_1hr"          "NEFM_pos_glut_1hr"
# [7] "NEFM_neg_glut_6hr" "NEFM_pos_glut_6hr" "unknown_6hr"      
# [10] "NEFM_pos_glut_0hr" "unknown_1hr"       "unknown_0hr"

# 0 hr only
combined_peaks_df <-
  reduce(x = c(peaks_uncombined[[3]],
               peaks_uncombined[[10]]))
combined_peaks_df <-
  combined_peaks_df[(width(combined_peaks_df) > 20) & 
                      (width(combined_peaks_df) < 10000)]
combined_peaks_df <- 
  combined_peaks_df[combined_peaks_df@seqnames@values == "chr1"]

# test to change the range of ATAC_chr1_glut to combined_peaks_df
# ATAC_chr1_glut_backup <- ATAC_chr1_glut
ATAC_chr1_glut <- ATAC_chr1_glut_backup
# ATAC_chr1_glut@assays$ATAC@ranges <- combined_peaks_df
# need to extract information from the original ATAC_chr1_glut object
# and re-quantify use the combined_peaks_df peakset

chr1_counts <-
  FeatureMatrix(fragments = ATAC_chr1_glut@assays$ATAC@fragments[[1]],
                features = combined_peaks_df,
                cells = ATAC_chr1_glut@assays$ATAC@fragments[[1]]@cells)
gene_ranges <- 
  GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86,
                      verbose = T)
seqlevelsStyle(gene_ranges) <- "UCSC"
chr1_glut_ATAC_new <-
  CreateChromatinAssay(counts = chr1_counts,
                       fragments = ATAC_chr1_glut@assays$ATAC@fragments[[1]],
                       # genome = "GRCh38",
                       annotation = gene_ranges)
genome(chr1_glut_ATAC_new) <-
  seqinfo(BSgenome.Hsapiens.UCSC.hg38)
# need to extract seqinfo from a BSgenome object and assign

chr1_glut_ATAC_new <-
  CreateSeuratObject(chr1_glut_ATAC_new,
                     assay = "ATAC")
CoveragePlot(chr1_glut_ATAC_new, 
             region = "chr1-150059035-150269039") 

# clearly, need to re-call peaks or try 10x peaks
## try 10x peaks first
library(readr)
atac_peaks_10x <- read_table2("/data/FASTQ/Duan_Project_024/hybrid_output/hybrid_aggr_5groups_no2ndary/outs/atac_peaks.bed", 
                          col_names = FALSE, skip = 58)
atac_peaks_10x$X4 <- "."
colnames(atac_peaks_10x) <- c("chr", "start", "end", "strand")
atac_peaks_10x <-
  makeGRangesFromDataFrame(atac_peaks_10x,
                           ignore.strand = T,
                           starts.in.df.are.0based = T)
# 10x peaks are too many
CoveragePlot(chr1_glut_ATAC_new, 
             region = "chr1-150059035-150269039",
             ranges = atac_peaks_10x,
             ranges.title = "atac_peaks_10x") 
DefaultAssay(chr1_glut_ATAC_new)

# call peaks use MACS2
MACS2_peaks <-
  CallPeaks(object = chr1_glut_ATAC_new,
            macs2.path = "~/Data/Anaconda3-envs/signac/bin/macs2",
            extsize = 150,
            shift = -73,
            additional.args = '-q 0.01 -g hs')
save.image("chr1_18_samples.RData")

CoveragePlot(chr1_glut_ATAC_new, 
             region = "chr1-150059035-150269039",
             ranges = MACS2_peaks,
             ranges.title = "MACS2_peaks") 


#### calculate links
# set a backup
# chr1_glut_ATAC_new_backup <- chr1_glut_ATAC_new

### remake Seurat object using new GRanges peakset (MACS2_peaks)

MACS2_peaks_chr1 <- MACS2_peaks[MACS2_peaks@seqnames == "chr1", ]
chr1_counts <-
  FeatureMatrix(fragments = ATAC_chr1_glut@assays$ATAC@fragments[[1]],
                features = MACS2_peaks_chr1,
                cells = ATAC_chr1_glut@assays$ATAC@fragments[[1]]@cells)

# !! important, extract annotations use this way
# !! then there will be no connect attempt to ftp.ncbi.nlm.nih.gov
gene_ranges <- EnsDb.Hsapiens.v86
ensembldb::seqlevelsStyle(gene_ranges) <- "UCSC"

annotation_ranges <-
  GetGRangesFromEnsDb(ensdb = gene_ranges,
                      verbose = T, 
                      standard.chromosomes = T)


chr1_glut_ATAC_new <-
  CreateChromatinAssay(counts = chr1_counts,
                       fragments = ATAC_chr1_glut@assays$ATAC@fragments[[1]],
                       annotation = annotation_ranges)


# need to extract seqinfo from a BSgenome object and assign

chr1_glut_ATAC_new <-
  CreateSeuratObject(chr1_glut_ATAC_new,
                     assay = "ATAC")
# chr1_glut_ATAC_new@assays$ATAC

CoveragePlot(chr1_glut_ATAC_new, 
             region = "chr1-150059035-150269039") 






# assign chr1_glut_ATAC_new to ATAC_cds 

# obj_chr16 <- obj_chr16@assays$ATAC
# obj_chr16@assays$ATAC@
# ATAC_cds <- as.cell_data_set(obj_chr1)
ATAC_cds <- as.cell_data_set(chr1_glut_ATAC_new)
# ATAC_cds <-
#   make_atac_cds(as.matrix(ATAC_cds@assayData$exprs),
#                 binarize = T)
ATAC_cds <- detect_genes(ATAC_cds)
ATAC_cds <- ATAC_cds[,Matrix::colSums(exprs(ATAC_cds)) != 0]
ATAC_cds <- estimate_size_factors(ATAC_cds)
ATAC_cds <- preprocess_cds(ATAC_cds, method = "LSI")
ATAC_cds <- 
  reduce_dimension(ATAC_cds, reduction_method = 'UMAP',
                   preprocess_method = "LSI")
ATAC_reduced_coord <- 
  reducedDims(ATAC_cds)
# ATAC_reduced_coord <- 
#   reducedDim

ATAC_cicero <- 
  make_cicero_cds(ATAC_cds,
                  reduced_coordinates = ATAC_reduced_coord$UMAP)

# find connections
# test_genomeinfo <- atac_small
# test_genomeinfo@assays$peaks@seqinfo@seqlengths[1] # use chr1
# seqlengths(test_genomeinfo)

genome_len <- 249250621
genome.df <- data.frame("chr" = "chr1",
                        "length" = genome_len)
conns <- 
  run_cicero(cds = ATAC_cicero,
             genomic_coords = genome.df,
             sample_num = 100)
ccans <- generate_ccans(conns, 
                        coaccess_cutoff_override = 0.05
)
ccans <- generate_ccans(conns)
# ccans <- generate_ccans(conns, 
#                         tolerance_digits = 3)


ATAC_links <- 
  ConnectionsToLinks(conns = conns,
                     ccans = ccans)
###
Links(chr1_glut_ATAC_new) <- ATAC_links
chr1_glut_ATAC_new$plot_ident <- "glut"

Idents(chr1_glut_ATAC_new) <- "plot_ident"
DefaultAssay(chr1_glut_ATAC_new) <- "ATAC"

CoveragePlot(chr1_glut_ATAC_new, 
             region = "chr1-150059035-150269039") 

# Links(obj_chr1) <- ATAC_links
# 
# Idents(obj_chr1) <- "time.ident"
# DefaultAssay(obj_chr1)
# 
# CoveragePlot(obj_chr1, 
#              region = "chr1-150059035-150269039") 
# use this way
gene_ranges <- EnsDb.Hsapiens.v86
ensembldb::seqlevelsStyle(gene_ranges) <- "UCSC"

annotation_ranges <-
  GetGRangesFromEnsDb(ensdb = gene_ranges,
                      verbose = T, 
                      standard.chromosomes = T)



link_rs2027349 <-
  ATAC_links[(ATAC_links@ranges@start > 150050000) &
               (ATAC_links@ranges@start < 150100000), ]
# link_df <-
#   cbind("start" = data.frame(link_rs2027349@ranges@start),
#         "end" = data.frame(link_rs2027349@ranges@start +
#                            link_rs2027349@ranges@width),
#         "strand" = data.frame(link_rs2027349@strand@values,
#                             stringsAsFactors = F),
#         "score" = data.frame(link_rs2027349@elementMetadata@listData$score))

link_df <-
  cbind("start" = link_rs2027349@ranges@start,
        "end" = link_rs2027349@ranges@start +
          link_rs2027349@ranges@width,
        "strand" = link_rs2027349@strand@values,
        "score" = link_rs2027349@elementMetadata@listData$score)

