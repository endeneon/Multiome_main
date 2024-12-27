# 13 Sept 2019
# function for plotting rs78710909 site (BIN1)
# revised from plot_anywhere and plot_ASoC_composite
# Use OverlayTrack to combine data tracks

# init
library(Gviz)
# data("cpgIslands")
library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ensembldb)
library(org.Hs.eg.db)
library(grDevices)
library(gridExtra)
library(GenomicRanges)
###
options(ucscChromosomeNames = F)

plot_BIN1 <- function(chr, start, end, gene_name = "", mcols = 100, strand = "+", x_offset_1 = 0, x_offset_2 = 0, ylimit = 800) {
  cell_type <-
    GRanges(seqnames = Rle(chr),
            seqinfo = Seqinfo(seqnames = chr,
                              genome = "hg38"),
            ranges = IRanges(start = start,
                             end = end,
                             names = chr),
            mcols = as.data.frame(mcols),
            strand = Rle(strand(strand)))

  print(as.character(unique(seqnames(cell_type))))

  iTrack <- IdeogramTrack(genome = genome(cell_type),
                          chromosome = as.character(unique(seqnames(cell_type))),
                          fontcolor = "black",
                          fontsize = 18)
  gTrack <- GenomeAxisTrack(col = "black",
                            fontcolor = "black",
                            fontsize = 16
                            # scale = 0.1 # not using scale, draw full axis
                            )
  alTrack_CN <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/Galaxy_Drive/BAMs/merged/Glut_all_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                transformation=function(x) {x * 1.0},
                                background.title = "darkblue",
                                fill.coverage = "darkblue",
                                col.coverage = "darkblue",
                                ylim = c(0, ylimit),
                                cex.axis = 1.2,
                                name = "iN-Glut")
  alTrack_NPC <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/Galaxy_Drive/BAMs/merged/NSC_all_sorted.bam",
                                 isPaired = T, coverageOnly = T,
                                 chromosome = as.character(unique(seqnames(cell_type))),
                                 genome = "hg38", type = "coverage",
                                 transformation=function(x) {x * 1.012},
                                 background.title = "goldenrod3",
                                 fill.coverage = "goldenrod3",
                                 col.coverage = "goldenrod3",
                                 ylim = c(0, ylimit),
                                 cex.axis = 1.2,
                                 name = "NPC")
  alTrack_DN <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/Galaxy_Drive/BAMs/merged/DN_all_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                transformation=function(x) {x * 0.685},
                                ylim = c(0, ylimit),
                                background.title = "darkorchid2",
                                fill.coverage = "darkorchid2",
                                col.coverage = "darkorchid2",
                                cex.axis = 1.2,
                                name = "iN-DN")
  alTrack_GA <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/Galaxy_Drive/BAMs/merged/GA_all_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                transformation=function(x) {x * 0.708},
                                ylim = c(0, ylimit),
                                background.title = "magenta",
                                fill.coverage = "magenta",
                                col.coverage = "magenta",
                                cex.axis = 1.2,
                                name = "iN-GA")
  alTrack_iPS <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/Galaxy_Drive/BAMs/merged/iPS_all_sorted.bam",
                                 isPaired = T, coverageOnly = T,
                                 chromosome = as.character(unique(seqnames(cell_type))),
                                 genome = "hg38", type = "coverage",
                                 transformation=function(x) {x * 0.613},
                                 ylim = c(0, ylimit),
                                 background.title = "brown",
                                 fill.coverage = "brown",
                                 col.coverage = "brown",
                                 cex.axis = 1.2,
                                 name = "iPS")

  snpTrack <- AnnotationTrack(start = 127107346, end = 127107346, chromosome = "chr2", # rs78710909
                              id = "rs78710909", shape = "box",
                              name = "SNP", strand = "*",
                              group = c("rs78710909"),
                              fontcolor.group = "black", fontcolor.item = "black",
                              fontsize = 16,
                              col = "black", col.title = "black",
                              just.group = "below",
                              showID = TRUE,
                              cex.group = 1,
                              groupAnnotation = "id")

  ########### plotting
  ucscGenes <- UcscTrack(genome=genome(cell_type), table="ncbiRefSeq",
                         track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                         chromosome = as.character(unique(seqnames(cell_type))),
                         rstarts = "exonStarts", rends = "exonEnds",
                         gene = "name", symbol = 'name', transcript = "name",
                         strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)
  z <- ranges(ucscGenes)
  mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "",
                                                             mcols(z)$symbol), "SYMBOL","REFSEQ"))
  grTrack <- ucscGenes
  ranges(grTrack) <- z
  grTrack@dp@pars$col.line <- "black"
  grTrack@dp@pars$fontcolor<- "black"
  grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
  grTrack@dp@pars$fontcolor.title <- "black"
  grTrack@dp@pars$fontcolor.item <- "black"
  grTrack@dp@pars$fontcolor.group <- "black"
  grTrack@dp@pars$fontsize.group <- 18

  ######

  # htTrack <- HighlightTrack(trackList = list(alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS),
  #                           start = c(26212000, 26241000),
  #                           width = c(15000, 2000)
  #                           chromosome = as.character(unique(seqnames(cell_type))))



  plotTracks(list(iTrack, gTrack, alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS, snpTrack, grTrack),
             sizes = c(0.5,0.5,1,1,1,1,1,0.5,1),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,

}


###
options(ucscChromosomeNames = F)

plot_CLU <- function(chr, start, end, gene_name = "", mcols = 100, strand = "+", x_offset_1 = 0, x_offset_2 = 0, ylimit = 800) {
  cell_type <-
    GRanges(seqnames = Rle(chr),
            seqinfo = Seqinfo(seqnames = chr,
                              genome = "hg38"),
            ranges = IRanges(start = start,
                             end = end,
                             names = chr),
            mcols = as.data.frame(mcols),
            strand = Rle(strand(strand)))

  print(as.character(unique(seqnames(cell_type))))

  iTrack <- IdeogramTrack(genome = genome(cell_type),
                          chromosome = as.character(unique(seqnames(cell_type))),
                          fontcolor = "black",
                          fontsize = 18)
  gTrack <- GenomeAxisTrack(col = "black",
                            fontcolor = "black",
                            fontsize = 16
                            # scale = 0.1 # not using scale, draw full axis
  )
  alTrack_CN <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/Galaxy_Drive/BAMs/merged/Glut_all_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                transformation=function(x) {x * 1.0},
                                background.title = "darkblue",
                                fill.coverage = "darkblue",
                                col.coverage = "darkblue",
                                ylim = c(0, ylimit),
                                cex.axis = 1.2,
                                name = "iN-Glut")
  alTrack_NPC <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/Galaxy_Drive/BAMs/merged/NSC_all_sorted.bam",
                                 isPaired = T, coverageOnly = T,
                                 chromosome = as.character(unique(seqnames(cell_type))),
                                 genome = "hg38", type = "coverage",
                                 transformation=function(x) {x * 1.012},
                                 background.title = "goldenrod3",
                                 fill.coverage = "goldenrod3",
                                 col.coverage = "goldenrod3",
                                 ylim = c(0, ylimit),
                                 cex.axis = 1.2,
                                 name = "NPC")
  alTrack_DN <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/Galaxy_Drive/BAMs/merged/DN_all_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                transformation=function(x) {x * 0.685},
                                ylim = c(0, ylimit),
                                background.title = "darkorchid2",
                                fill.coverage = "darkorchid2",
                                col.coverage = "darkorchid2",
                                cex.axis = 1.2,
                                name = "iN-DN")
  alTrack_GA <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/Galaxy_Drive/BAMs/merged/GA_all_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                transformation=function(x) {x * 0.708},
                                ylim = c(0, ylimit),
                                background.title = "magenta",
                                fill.coverage = "magenta",
                                col.coverage = "magenta",
                                cex.axis = 1.2,
                                name = "iN-GA")
  alTrack_iPS <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/Galaxy_Drive/BAMs/merged/iPS_all_sorted.bam",
                                 isPaired = T, coverageOnly = T,
                                 chromosome = as.character(unique(seqnames(cell_type))),
                                 genome = "hg38", type = "coverage",
                                 transformation=function(x) {x * 0.613},
                                 ylim = c(0, ylimit),
                                 background.title = "brown",
                                 fill.coverage = "brown",
                                 col.coverage = "brown",
                                 cex.axis = 1.2,
                                 name = "iPS")

  snpTrack <- AnnotationTrack(start = 27608798, end = 27608798, chromosome = "chr8", # rs78710909
                              id = "rs1532278", shape = "box",
                              name = "SNP", strand = "*",
                              group = c("rs1532278"),
                              fontcolor.group = "black", fontcolor.item = "black",
                              fontsize = 16,
                              col = "black", col.title = "black",
                              just.group = "below",
                              showID = TRUE,
                              cex.group = 1,
                              groupAnnotation = "id")

  ########### plotting
  ucscGenes <- UcscTrack(genome=genome(cell_type), table="ncbiRefSeq",
                         track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                         chromosome = as.character(unique(seqnames(cell_type))),
                         rstarts = "exonStarts", rends = "exonEnds",
                         gene = "name", symbol = 'name', transcript = "name",
                         strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)
  z <- ranges(ucscGenes)
  mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "",
                                                             mcols(z)$symbol), "SYMBOL","REFSEQ"))
  grTrack <- ucscGenes
  ranges(grTrack) <- z
  grTrack@dp@pars$col.line <- "black"
  grTrack@dp@pars$fontcolor<- "black"
  grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
  grTrack@dp@pars$fontcolor.title <- "black"
  grTrack@dp@pars$fontcolor.item <- "black"
  grTrack@dp@pars$fontcolor.group <- "black"
  grTrack@dp@pars$fontsize.group <- 18

  ######

  # htTrack <- HighlightTrack(trackList = list(alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS),
  #                           start = c(26212000, 26241000),
  #                           width = c(15000, 2000)
  #                           chromosome = as.character(unique(seqnames(cell_type))))



  plotTracks(list(iTrack, gTrack, alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS, snpTrack, grTrack),
             sizes = c(0.5,2,1,1,1,1,1,0.5,1),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,

}
