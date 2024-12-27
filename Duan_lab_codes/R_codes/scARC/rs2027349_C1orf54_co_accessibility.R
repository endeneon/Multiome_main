# 31 March 2022
# check co-accessible sites between rs2027349 and C1orf54

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
library(stringr)

library(ggplot2)
library(patchwork)

## init
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/ATACseq_5line_objonly_labeled.RData")
obj_complete@assays$RNA <- NULL
obj_complete@assays$SCT <- NULL

Idents(obj_complete)
DefaultAssay(obj_complete)
# obj_complete@assays$ATAC@data@Dim



# obj_chr1 <-
#   obj_complete[str_detect(rownames(obj_complete),
#                           pattern = "^chr1\\-"), ]
# 
# # obj_chr1$broad.cell.type
# obj_chr1 <- obj_chr1[, !is.na(obj_chr1$fine.cell.type)]
# obj_chr1 <- obj_chr1[, obj_chr1$fine.cell.type %in%
#                          c("NEFM+/CUX2- glut", "NEFM-/CUX2+ glut", "GABA",
#                            "NPC")]
# save(obj_chr1, file = "obj_chr1_cleaned.RData")
# rm(obj_complete)
# gc()


load("obj_chr1_cleaned.RData")
# obj_chr1_backup <- obj_chr1


# test regression using 0vs1 or 0vs6 condition, by cell type
# NEFM+ 0vs1/0vs6 hr
# obj_chr1 <-
#   obj_chr1_backup[, (obj_chr1_backup$fine.cell.type %in% "NEFM+/CUX2- glut") &
#                     ((obj_chr1_backup$time.ident %in% "0hr") |
#                        (obj_chr1_backup$time.ident %in% "6hr"))]

obj_chr1 <-
  obj_chr1_backup[, (obj_chr1_backup$fine.cell.type %in% "NEFM-/CUX2+ glut") &
                    ((obj_chr1_backup$time.ident %in% "0hr") |
                       (obj_chr1_backup$time.ident %in% "1hr"))]


# GABA 0vs1/0vs6 hr
# obj_chr1 <-
#   obj_chr1_backup[, (obj_chr1_backup$fine.cell.type %in% "GABA") &
#                     ((obj_chr1_backup$time.ident %in% "6hr") |
#                        (obj_chr1_backup$time.ident %in% "1hr"))]

# obj_chr1 <-
#   obj_chr1_backup[, (obj_chr1_backup$fine.cell.type %in% "NPC") &
#                     ((obj_chr1_backup$time.ident %in% "0hr") |
#                        (obj_chr1_backup$time.ident %in% "6hr"))]


# obj_chr16 <- obj_chr16@assays$ATAC
# obj_chr16@assays$ATAC@
ATAC_cds <- as.cell_data_set(obj_chr1)
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
# ccans <- generate_ccans(conns, 
#                         tolerance_digits = 3)

ATAC_links <- 
  ConnectionsToLinks(conns = conns,
                     ccans = ccans)
Links(obj_chr1) <- ATAC_links

Idents(obj_chr1) <- "time.ident"
DefaultAssay(obj_chr1)

CoveragePlot(obj_chr1, 
             region = "chr1-150059035-150269039") 


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

## extract barcodes for cell subsetting
library(stringr)

load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/ATAC_new_obj_with_MACS2_new_peaks_and_RNA_labels.RData")
obj_complete <- ATAC_new
# obj_complete@assays$ATAC <- NULL
# obj_complete@assays$peaks <- NULL
# DefaultAssay(obj_complete) <- "RNA"

unique(obj_complete$group.ident)
# [1] "51-0" "35-6" "5-0"  "33-1" "5-1"  "35-0" "63-1" "51-1" "33-6" "63-6" "51-6" "33-0" "5-6" 
# [14] "35-1" "63-0"
Idents(obj_complete) <- "group.ident"
obj_complete <- 
  obj_complete[, obj_complete$RNA.cell.type %in% 
                 c("NEFM_pos_glut", "NEFM_neg_glut", "GABA")]

obj_complete@active.assay
obj_complete@active.ident
sum(is.na(obj_complete@active.ident))

cell_types <- unique(obj_complete$RNA.cell.type)
cell_types <- as.character(cell_types)

library_name <- unique(obj_complete$group.ident)
# library_name <-
#   str_replace_all(string = library_name,
#                   pattern = '-',
#                   replacement = '_')

# test <- obj_complete@assays$RNA@data@Dimnames[[2]]


i <- ""
for (i in cell_types) {
  cell_barcodes_writeout <-
    obj_complete@assays$RNA@data@Dimnames[[2]][obj_complete@active.ident == i]
  cell_barcodes_writeout <-
    str_split(string = cell_barcodes_writeout,
              pattern = '-',
              simplify = T)[, 1]
  cell_barcodes_writeout <-
    paste(cell_barcodes_writeout,
          '-1', sep = "")
  cell_barcodes_writeout <-
    unique(sort(cell_barcodes_writeout))
  write.table(cell_barcodes_writeout,
              file = paste0('barcodes_4_demux/', i, '_barcodes.txt'),
              quote = F, sep = "\n",
              row.names = F, col.names = F)
}


### barcodes by library
### pre-create "GABA"          "NEFM_pos_glut" "NEFM_neg_glut"
# folder

i <- ""
j <- ""

for (i in cell_types) {
  for (j in library_name) {
    cell_barcodes_writeout <-
      obj_complete@assays$ATAC@data@Dimnames[[2]][(obj_complete$RNA.cell.type == i) &
                                                    (obj_complete$group.ident == j)]
    cell_barcodes_writeout <-
      str_split(string = cell_barcodes_writeout,
                pattern = '-',
                simplify = T)[, 1]
    cell_barcodes_writeout <-
      paste(cell_barcodes_writeout,
            '-1', sep = "")
    cell_barcodes_writeout <-
      unique(sort(cell_barcodes_writeout))
    write.table(cell_barcodes_writeout,
                file = paste0('barcodes_by_library/', 
                              i, '/', i, '_', j,
                              '_barcodes.txt'),
                quote = F, sep = "\n",
                row.names = F, col.names = F)
  }
}
