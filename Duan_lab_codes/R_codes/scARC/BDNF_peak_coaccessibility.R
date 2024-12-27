# Siwei 18 Apw 2022
# Check the correlation of BDNF peak


# need chr11 only

#init

{
  library(Seurat)
  library(SeuratWrappers)
  library(Signac)
  library(future)
  
  library(remotes)
  # library(monocle3)
  # library(monocle)
  # library(cicero)
  
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(EnsDb.Hsapiens.v86)
  library(stringr)
  
  library(ggplot2)
  library(patchwork)
}

plan("multisession", workers = 4)
options(future.globals.maxSize = 100 * 1024^3)

# load data
load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/18line_10xaggred_labeled_by_RNA_obj.RData")

ATAC_chr11 <-
  ATAC[str_detect(string = rownames(ATAC),
                  pattern = "chr11\\-"), ]
rm(ATAC)

# keep glut only, need to compare
# all gluts
# 2 glut types
# all 0 hr gluts
# glut 0 hr types
ATAC_chr11_glut <-
  ATAC_chr11[, ATAC_chr11$cell.type %in% 
              c("NEFM_neg_glut", "NEFM_pos_glut")]

ATAC_chr11$glut.ident <- "glut"
Idents(ATAC_chr11_glut) <- "glut.ident"

# call peaks use MACS2
MACS2_peaks <-
  CallPeaks(object = ATAC_chr11_glut,
            macs2.path = "~/Data/Anaconda3-envs/signac/bin/macs2",
            extsize = 150,
            shift = -73,
            additional.args = '-q 0.001 -g hs')

CoveragePlot(ATAC_chr11_glut, 
             region = "chr11-27700000-27780000",
             # ranges = MACS2_peaks,
             ranges.title = "MACS2_peaks") 
# !! important, extract annotations use this way
# !! then there will be no connect attempt to ftp.ncbi.nlm.nih.gov
gene_ranges <- EnsDb.Hsapiens.v86
ensembldb::seqlevelsStyle(gene_ranges) <- "UCSC"
annotation_ranges <-
  GetGRangesFromEnsDb(ensdb = gene_ranges,
                      verbose = T, 
                      standard.chromosomes = T)

Annotation(ATAC_chr11_glut@assays$ATAC) <- annotation_ranges
rm(ATAC_chr11)

# ATAC_chr11_glut_backup <- ATAC_chr11_glut

ATAC_cds <- as.cell_data_set(ATAC_chr11_glut)
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
# test_genomeinfo@assays$peaks@seqinfo@seqlengths[1] # use chr11
# seqlengths(test_genomeinfo)

genome_len <- 135086622 # chr11
genome.df <- data.frame("chr" = "chr11",
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
Links(ATAC_chr11_glut) <- ATAC_links
ATAC_chr11_glut$plot_ident <- "glut"

Idents(ATAC_chr11_glut) <- "plot_ident"
DefaultAssay(ATAC_chr11_glut) <- "ATAC"

CoveragePlot(ATAC_chr11_glut, 
             region = "chr11-27700000-27780000") 
