# function for plotting rs2027349 site
# revised from plot_anywhere
# Use OverlayTrack to combine data tracks

plot_AsoC_composite <- function(chr, start, end, gene_name = "",
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
                            fontsize = 14,
                            scale = 0.1)

  CN_alTrack_Region <- AlignmentsTrack(paste("rs2027349_het_05Mar2019/", "CN", "_het_rs2027349.bam",
                                          collapse = "", sep = ""),
                                    isPaired = T, coverageOnly = T,
                                    chromosome = as.character(unique(seqnames(cell_type))),
                                    genome = "hg38", type = "coverage",
                                    background.title = "darkorange",
                                    fill.coverage = "darkorange",
                                    col.coverage = "darkorange",
                                    ylim = c(0, ylimit),
                                    cex.axis = 1,
                                    col.title="black",
                                    col.axis="black",
                                    alpha = 0.9,
                                    name = "Glut",
                                    fontsize = 14,
                                    show.title = FALSE)

  CN_alTrack_A <- AlignmentsTrack(paste("rs2027349_het_05Mar2019/SNP_split_output/", "CN", "_het_rs2027349.ref.A.bam",
                                     collapse = "", sep = ""),
                               isPaired = T, coverageOnly = T,
                               chromosome = as.character(unique(seqnames(cell_type))),
                               genome = "hg38", type = "coverage",
                               background.title = "darkred",
                               fill.coverage = "darkred",
                               col.coverage = "darkred",
                               ylim = c(0, ylimit),
                               cex.axis = 1,
                               alpha = 0.9,
                               # col.title="black",
                               name = "Allele_A")

  CN_alTrack_G <- AlignmentsTrack(paste("rs2027349_het_05Mar2019/SNP_split_output/", "CN", "_het_rs2027349.alt.G.bam",
                                     collapse = "", sep = ""),
                               isPaired = T, coverageOnly = T,
                               chromosome = as.character(unique(seqnames(cell_type))),
                               genome = "hg38", type = "coverage",
                               background.title = "darkblue",
                               fill.coverage = "darkblue",
                               col.coverage = "darkblue",
                               ylim = c(0, ylimit),
                               alpha = 0.9,
                               cex.axis = 1,
                               # col.title="black",
                               name = "Allele_G")

  CN_OvlTrack <- OverlayTrack(trackList = list(CN_alTrack_Region, CN_alTrack_G, CN_alTrack_A),
                           name = "REF_ALT_REGION",
                           fontsize = 14)

  NPC_alTrack_Region <- AlignmentsTrack(paste("rs2027349_het_05Mar2019/", "NPC", "_het_rs2027349.bam",
                                             collapse = "", sep = ""),
                                       isPaired = T, coverageOnly = T,
                                       chromosome = as.character(unique(seqnames(cell_type))),
                                       genome = "hg38", type = "coverage",
                                       background.title = "darkolivegreen4",
                                       fill.coverage = "darkolivegreen4",
                                       col.coverage = "darkolivegreen4",
                                       ylim = c(0, ylimit),
                                       cex.axis = 1,
                                       col.title="black",
                                       col.axis="black",
                                       alpha = 0.9,
                                       name = "NPC",
                                       fontsize = 14,
                                       show.title = FALSE)

  NPC_alTrack_A <- AlignmentsTrack(paste("rs2027349_het_05Mar2019/SNP_split_output/", "NPC", "_het_rs2027349.ref.A.bam",
                                        collapse = "", sep = ""),
                                  isPaired = T, coverageOnly = T,
                                  chromosome = as.character(unique(seqnames(cell_type))),
                                  genome = "hg38", type = "coverage",
                                  background.title = "darkred",
                                  fill.coverage = "darkred",
                                  col.coverage = "darkred",
                                  ylim = c(0, ylimit),
                                  alpha = 0.9,
                                  cex.axis = 1,
                                  # col.title="black",
                                  name = "Allele_A")

  NPC_alTrack_G <- AlignmentsTrack(paste("rs2027349_het_05Mar2019/SNP_split_output/", "NPC", "_het_rs2027349.alt.G.bam",
                                        collapse = "", sep = ""),
                                  isPaired = T, coverageOnly = T,
                                  chromosome = as.character(unique(seqnames(cell_type))),
                                  genome = "hg38", type = "coverage",
                                  background.title = "darkblue",
                                  fill.coverage = "darkblue",
                                  col.coverage = "darkblue",
                                  alpha = 0.9,
                                  ylim = c(0, ylimit),
                                  cex.axis = 1,
                                  # col.title="black",
                                  name = "Allele_G")

  NPC_OvlTrack <- OverlayTrack(trackList = list(NPC_alTrack_Region, NPC_alTrack_G, NPC_alTrack_A),
                              name = "REF_ALT_REGION",
                              fontsize = 14)

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

  htTrack <- HighlightTrack(trackList = list(CN_OvlTrack, NPC_OvlTrack),
                            start = 150067560,
                            width = 120,
                            chromosome = as.character(unique(seqnames(cell_type))))



  plotTracks(list(iTrack, gTrack, CN_OvlTrack, NPC_OvlTrack, snpTrack, grTrack),
             sizes = c(1,1,4,4,1,1),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,

}
