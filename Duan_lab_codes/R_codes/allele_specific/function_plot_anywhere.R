# function for plotting
plot_anywhere <- function(chr, start, end, gene_name = "",
                          SNPname = "", SNPposition = 1L,
                          mcols = 100, strand = "+",
                          x_offset_1 = 0, x_offset_2 = 0, ylimit = 800) {
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
                            fontsize = 16,
                            scale = 0.1)
  alTrack_CN <- AlignmentsTrack("~/3TB/merged_BAM/CN_all_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                transformation=function(x) {x * 4.078446},
                                background.title = "darkblue",
                                fill.coverage = "darkblue",
                                col.coverage = "darkblue",
                                ylim = c(0, ylimit),
                                cex.axis = 1.2,
                                name = "iN-Glut")
  alTrack_NPC <- AlignmentsTrack("~/3TB/merged_BAM/NSC_all_sorted.bam",
                                 isPaired = T, coverageOnly = T,
                                 chromosome = as.character(unique(seqnames(cell_type))),
                                 genome = "hg38", type = "coverage",
                                 transformation=function(x) {x * 4.21563},
                                 background.title = "goldenrod3",
                                 fill.coverage = "goldenrod3",
                                 col.coverage = "goldenrod3",
                                 ylim = c(0, ylimit),
                                 cex.axis = 1.2,
                                 name = "NPC")
  alTrack_DN <- AlignmentsTrack("~/3TB/merged_BAM/DN_all_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                transformation=function(x) {x * 9.182318},
                                ylim = c(0, ylimit),
                                background.title = "darkorchid2",
                                fill.coverage = "darkorchid2",
                                col.coverage = "darkorchid2",
                                cex.axis = 1.2,
                                name = "iN-DN")
  alTrack_GA <- AlignmentsTrack("~/3TB/merged_BAM/GA_all_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                transformation=function(x) {x * 5.278966},
                                ylim = c(0, ylimit),
                                background.title = "magenta",
                                fill.coverage = "magenta",
                                col.coverage = "magenta",
                                cex.axis = 1.2,
                                name = "iN-GA")
  alTrack_iPS <- AlignmentsTrack("~/3TB/merged_BAM/ips_all_sorted.bam",
                                 isPaired = T, coverageOnly = T,
                                 chromosome = as.character(unique(seqnames(cell_type))),
                                 genome = "hg38", type = "coverage",
                                 transformation=function(x) {x * 4.560793},
                                 ylim = c(0, ylimit),
                                 background.title = "brown",
                                 fill.coverage = "brown",
                                 col.coverage = "brown",
                                 cex.axis = 1.2,
                                 name = "iPS")

  snpTrack <- AnnotationTrack(start = SNPposition, end = SNPposition, chromosome = chr,
                              id = SNPname, shape = "box",
                              name = SNPname, strand = "*",
                              group = SNPname,
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



  plotTracks(list(iTrack, gTrack, alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS,
                  snpTrack, grTrack),
             sizes = c(0.5, 0.5, 1,1,1,1,1, 0.5, 0.5),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,

}
