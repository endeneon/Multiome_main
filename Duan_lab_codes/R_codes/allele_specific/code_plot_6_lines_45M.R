# 20 Dec 2019 Siwei
# plot all 6 AA/AG/GG isogenic lines of VPS45

plot_VPS45_6_lines_q20 <- function(chr, start, end, gene_name = "VPS45",
                                   mcols = 100, strand = "+",
                                   x_offset_1 = 0, x_offset_2 = 0, ylimit = 400,
                                   transparency = 1,
                                   title_name = "") {
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
                            fontsize = 14,
                            scale = 0.1)

  AA_1_track_region <- AlignmentsTrack("VPS45_isogenic/Q20/down_45M/AA_A11_WASPed.45M.q20.bam",
                                       isPaired = T, coverageOnly = T,
                                       chromosome = as.character(unique(seqnames(cell_type))),
                                       genome = "hg38", type = "coverage",
                                       background.title = "transparent",
                                       fill.coverage = "orangered",
                                       col.coverage = "orangered",
                                       transformation=function(x) {x * 2},
                                       ylim = c(0, ylimit),
                                       cex.axis = 1,
                                       col.title="black",
                                       col.axis="black",
                                       alpha = transparency,
                                       name = "AA #1",
                                       fontsize = 10,
                                       show.title = FALSE)
  AA_2_track_region <- AlignmentsTrack("VPS45_isogenic/Q20/down_45M/AA_H12_WASPed.45M.q20.bam",
                                       isPaired = T, coverageOnly = T,
                                       chromosome = as.character(unique(seqnames(cell_type))),
                                       genome = "hg38", type = "coverage",
                                       background.title = "transparent",
                                       fill.coverage = "orangered",
                                       col.coverage = "orangered",
                                       transformation=function(x) {x * 2},
                                       ylim = c(0, ylimit),
                                       cex.axis = 1,
                                       col.title="black",
                                       col.axis="black",
                                       alpha = transparency,
                                       name = "AA #2",
                                       fontsize = 10,
                                       show.title = FALSE)
  GG_1_track_region <- AlignmentsTrack("VPS45_isogenic/Q20/down_45M/GG_B11_WASPed.45M.q20.bam",
                                       isPaired = T, coverageOnly = T,
                                       chromosome = as.character(unique(seqnames(cell_type))),
                                       genome = "hg38", type = "coverage",
                                       background.title = "transparent",
                                       fill.coverage = "blue",
                                       col.coverage = "blue",
                                       transformation=function(x) {x * 2},
                                       ylim = c(0, ylimit),
                                       cex.axis = 1,
                                       col.title="black",
                                       col.axis="black",
                                       alpha = transparency,
                                       name = "GG #1",
                                       fontsize = 10,
                                       show.title = FALSE)
  GG_2_track_region <- AlignmentsTrack("VPS45_isogenic/Q20/down_45M/GG_G7_WASPed.45M.q20.bam",
                                       isPaired = T, coverageOnly = T,
                                       chromosome = as.character(unique(seqnames(cell_type))),
                                       genome = "hg38", type = "coverage",
                                       background.title = "transparent",
                                       fill.coverage = "blue",
                                       col.coverage = "blue",
                                       transformation=function(x) {x * 2},
                                       ylim = c(0, ylimit),
                                       cex.axis = 1,
                                       col.title="black",
                                       col.axis="black",
                                       alpha = transparency,
                                       name = "GG #2",
                                       fontsize = 10,
                                       show.title = FALSE)
  AG_1_track_region <- AlignmentsTrack("VPS45_isogenic/Q20/down_45M/CD11_1_WASPed.45M.q20.bam",
                                       isPaired = T, coverageOnly = T,
                                       chromosome = as.character(unique(seqnames(cell_type))),
                                       genome = "hg38", type = "coverage",
                                       background.title = "transparent",
                                       fill.coverage = "forestgreen",
                                       col.coverage = "forestgreen",
                                       transformation=function(x) {x * 2},
                                       ylim = c(0, ylimit),
                                       cex.axis = 1,
                                       col.title="black",
                                       col.axis="black",
                                       alpha = transparency,
                                       name = "Het #1",
                                       fontsize = 10,
                                       show.title = FALSE)
  AG_2_track_region <- AlignmentsTrack("VPS45_isogenic/Q20/down_45M/CD11_2_WASPed.45M.q20.bam",
                                       isPaired = T, coverageOnly = T,
                                       chromosome = as.character(unique(seqnames(cell_type))),
                                       genome = "hg38", type = "coverage",
                                       background.title = "transparent",
                                       fill.coverage = "forestgreen",
                                       col.coverage = "forestgreen",
                                       transformation=function(x) {x * 2},
                                       ylim = c(0, ylimit),
                                       cex.axis = 1,
                                       col.title="black",
                                       col.axis="black",
                                       alpha = transparency,
                                       name = "Het #2",
                                       fontsize = 10,
                                       show.title = FALSE)

  ###
  # test_track <- AG_1_track_region + AG_2_track_region

  snpTrack <- AnnotationTrack(start = 150067621, end = 150067621, chromosome = "chr1",
                              id = "rs2027349", shape = "box",
                              name = "SNP", strand = "*",
                              group = c("rs2027349"),
                              fontcolor.group = "black", fontcolor.item = "black",
                              fontsize = 16,
                              col = "black", col.title = "black",
                              just.group = "below",
                              showID = TRUE,
                              cex.group = 1,
                              groupAnnotation = "id")
  ###
  htTrack <- HighlightTrack(trackList = list(AA_1_track_region, AA_2_track_region,
                                             AG_1_track_region, AG_2_track_region,
                                             GG_1_track_region, GG_2_track_region),
                            start = c(150067621),
                            width = c(1),
                            col = "grey30",
                            inBackground = F,
                            chromosome = as.character(unique(seqnames(cell_type))))
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
  ###

  plotTracks(list(iTrack, gTrack,
                  htTrack,
                  # test_track,
                  # AA_1_track_region, AA_2_track_region,
                  # AG_1_track_region, AG_2_track_region,
                  # GG_1_track_region, GG_2_track_region,
                  snpTrack),
             sizes = c(0.5,0.5,
                       # 12,
                       2,2,2,2,2,2,
                       0.5),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,


}
