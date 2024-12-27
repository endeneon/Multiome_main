# function for plotting BCL11B site
# SNPname = "rs12895055", chr = 14, SNPposition = 99246457,
# revised from plot_anywhere
# Use OverlayTrack to combine data tracks

plot_AsoC_BCL11B_2sites <- function(chr, start, end, gene_name = "",
                                       mcols = 100, strand = "+",
                                       x_offset_1 = 0, x_offset_2 = 0, ylimit = 400,
                                       cell_type_to_plot = "CN", title_name = "") {
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
                            fontsize = 14)

  CN_alTrack_Region <- AlignmentsTrack("BCL11B/BCL11B_het_merged_sorted.bam",
                                       isPaired = T, coverageOnly = T,
                                       chromosome = as.character(unique(seqnames(cell_type))),
                                       genome = "hg38", type = "coverage",
                                       background.title = "black",
                                       fill.coverage = "forestgreen",
                                       col.coverage = "forestgreen",
                                       transformation=function(x) {x * 4.078446},
                                       ylim = c(0, ylimit),
                                       cex.axis = 1,
                                       col.title="black",
                                       col.axis="black",
                                       # alpha = 0.9,
                                       name = "Glut",
                                       fontsize = 14,
                                       show.title = FALSE)

  rs12895055_C <- AlignmentsTrack("BCL11B/rs12895055.ref.C.bam",
                                  isPaired = T, coverageOnly = T,
                                  chromosome = as.character(unique(seqnames(cell_type))),
                                  genome = "hg38", type = "coverage",
                                  background.title = "black",
                                  fill.coverage = "darkred",
                                  col.coverage = "darkred",
                                  transformation=function(x) {x * 4.078446},
                                  ylim = c(0, ylimit),
                                  cex.axis = 1,
                                  alpha = 0.9,
                                  col.title="black",
                                  col.axis="black",
                                  name = "Allele_C")

  rs12895055_T <- AlignmentsTrack("BCL11B/rs12895055.alt.T.bam",
                                  isPaired = T, coverageOnly = T,
                                  chromosome = as.character(unique(seqnames(cell_type))),
                                  genome = "hg38", type = "coverage",
                                  background.title = "black",
                                  fill.coverage = "darkblue",
                                  col.coverage = "darkblue",
                                  transformation=function(x) {x * 4.078446},
                                  ylim = c(0, ylimit),
                                  alpha = 0.9,
                                  cex.axis = 1,
                                  col.title="black",
                                  name = "Allele_T")
  rs11624408_A <- AlignmentsTrack("BCL11B/rs11624408.ref.A.bam",
                                  isPaired = T, coverageOnly = T,
                                  chromosome = as.character(unique(seqnames(cell_type))),
                                  genome = "hg38", type = "coverage",
                                  background.title = "black",
                                  fill.coverage = "cyan",
                                  col.coverage = "cyan",
                                  transformation=function(x) {x * 4.078446},
                                  ylim = c(0, ylimit),
                                  cex.axis = 1,
                                  alpha = 0.9,
                                  col.title="black",
                                  col.axis="black",
                                  name = "Allele_A")

  rs11624408_G <- AlignmentsTrack("BCL11B/rs11624408.alt.G.bam",
                                  isPaired = T, coverageOnly = T,
                                  chromosome = as.character(unique(seqnames(cell_type))),
                                  genome = "hg38", type = "coverage",
                                  background.title = "black",
                                  fill.coverage = "darkorange",
                                  col.coverage = "darkorange",
                                  transformation=function(x) {x * 4.078446},
                                  ylim = c(0, ylimit),
                                  alpha = 0.9,
                                  cex.axis = 1,
                                  col.title="black",
                                  name = "Allele_G")

  CN_OvlTrack <- OverlayTrack(trackList = list(CN_alTrack_Region,
                                               rs12895055_C, rs12895055_T,
                                               rs11624408_A, rs11624408_G),
                              name = "REF_ALT_REGION",
                              fontsize = 14)

  snpTrack <- AnnotationTrack(start = c(99246457, 99246608),
                              end = c(99246457, 99246608),
                              chromosome = c("chr14", "chr14"),
                              id = c("rs12895055", "rs11624408"),
                              shape = "box",
                              name = "SNP", strand = "*",
                              group = c("rs12895055", "rs11624408"),
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



  plotTracks(list(iTrack, gTrack, CN_OvlTrack,
                  # NPC_OvlTrack,
                  snpTrack),
             sizes = c(0.5,0.5,4,0.5),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,

}
