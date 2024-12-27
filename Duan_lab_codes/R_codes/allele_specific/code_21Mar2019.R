# Siwei 21 Mar 2019

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
# data("twoGroups")

# FigS10e$..1 <- NULL

colnames(FigS10e) <- c("Geneid", "Gene_Symbol", "logFC", "logCPM", "LR", "Pvalue", "FDR", "CHR", "Strand", "txSTART", "txSTOP")


VPS45_plot <- makeGRangesFromDataFrame(FigS10e,
                                       keep.extra.columns = T,
                                       seqnames.field = FigS10e$CHR,
                                       start.field = as.character(FigS10e$txSTART),
                                       end.field = as.character(FigS10e$txSTOP),
                                       strand.field = as.character(FigS10e$Strand))

VPS45_plot <- makeGRangesFromDataFrame(FigS10e,
                                       keep.extra.columns = T,
                                       seqnames.field = "CHR",
                                       start.field = "txStart",
                                       end.field = "txStop",
                                       strand.field = "Strand")
VPS45_plot@seqinfo@genome <- "hg38"

iTrack <- IdeogramTrack(genome = genome(VPS45_plot),
                        chromosome = as.character(unique(seqnames(VPS45_plot))),
                        fontcolor = "black",
                        fontsize = 18)

gTrack <- GenomeAxisTrack(col = "black",
                          fontcolor = "black",
                          fontsize = 14,
                          scale = 0.1)

expression_fold <- DataTrack(VPS45_plot,
                             data = VPS45_plot$logFC,
                             strand = "+",
                             type = "h"
                             # col = rep(c("red", "green"), length.out = length(VPS45_plot$logFC/2))
                             )

ucscGenes <- UcscTrack(genome=genome(VPS45_plot), table="ncbiRefSeq",
                       track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome = as.character(unique(seqnames(VPS45_plot))),
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



plotTracks(iTrack, gTrack, expression_fold, grTrack)

plotTracks(grTrack)

##########


plot_anywhere(chr = "chr1", start = 149887589, end = 150476866, ylimit = 8000)


plot_AsoC_peaks_rev(chr = "chr1", start = 149887589, end = 150476866, ylim = 1000,
                cell_type_to_plot = "NPC", title_name = " ")
