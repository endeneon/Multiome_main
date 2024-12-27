# 17 Dec 2019 Siwei
# plot AA/AG/GG isogenic lines of VPS45
# use the merged file downsampled to 45M

plot_VPS45_isogenic_45M <- function(chr, start, end, gene_name = "VPS45",
                                mcols = 100, strand = "+",
                                x_offset_1 = 0, x_offset_2 = 0, ylimit = 400,
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
                            fontsize = 14
                            # scale = 0.1
                            )

  AA_track_region <- AlignmentsTrack("VPS45_isogenic/Q20/down_45M/temp/AA_45M_q20_merged_sorted.bam",
                                     isPaired = T, coverageOnly = T,
                                     chromosome = as.character(unique(seqnames(cell_type))),
                                     genome = "hg38", type = "coverage",
                                     background.title = "black",
                                     fill.coverage = "orangered",
                                     col.coverage = "orangered",
                                     transformation=function(x) {x * 2},
                                     ylim = c(0, ylimit),
                                     cex.axis = 1,
                                     col.title="black",
                                     col.axis="black",
                                     alpha = 0.7,
                                     name = "in NPC",
                                     fontsize = 10,
                                     show.title = FALSE)
  GG_track_region <- AlignmentsTrack("VPS45_isogenic/Q20/down_45M/temp/GG_45M_q20_merged_sorted.bam",
                                     isPaired = T, coverageOnly = T,
                                     chromosome = as.character(unique(seqnames(cell_type))),
                                     genome = "hg38", type = "coverage",
                                     background.title = "black",
                                     fill.coverage = "blue",
                                     col.coverage = "blue",
                                     transformation=function(x) {x * 2},
                                     ylim = c(0, ylimit),
                                     cex.axis = 1,
                                     col.title="black",
                                     col.axis="black",
                                     alpha = 1,
                                     name = "in NPC",
                                     fontsize = 10,
                                     show.title = FALSE)
  AG_track_region <- AlignmentsTrack("VPS45_isogenic/Q20/down_45M/temp/CD_Het_45M_q20_merged_sorted.bam",
                                     isPaired = T, coverageOnly = T,
                                     chromosome = as.character(unique(seqnames(cell_type))),
                                     genome = "hg38", type = "coverage",
                                     background.title = "cyan",
                                     fill.coverage = "cyan",
                                     col.coverage = "cyan",
                                     transformation=function(x) {x * 2},
                                     ylim = c(0, ylimit),
                                     cex.axis = 1,
                                     col.title="black",
                                     col.axis="black",
                                     alpha = 0.5,
                                     name = "in NPC",
                                     fontsize = 10,
                                     show.title = FALSE)
  Het_AG_track_region <- AlignmentsTrack("VPS45_isogenic/Q20/down_45M/temp/CD_Het_45M_q20_merged_sorted.bam",
                                         isPaired = T, coverageOnly = T,
                                         chromosome = as.character(unique(seqnames(cell_type))),
                                         genome = "hg38", type = "coverage",
                                         background.title = "black",
                                         fill.coverage = "forestgreen",
                                         col.coverage = "forestgreen",
                                         transformation=function(x) {x * 2},
                                         ylim = c(0, ylimit),
                                         cex.axis = 1,
                                         col.title="black",
                                         col.axis="black",
                                         alpha = 0.9,
                                         name = "in Het",
                                         fontsize = 10,
                                         show.title = FALSE)
  Het_A_track_region <- AlignmentsTrack("VPS45_isogenic/Q20/down_45M/temp/CD_Het.ref.A.bam",
                                        isPaired = T, coverageOnly = T,
                                        chromosome = as.character(unique(seqnames(cell_type))),
                                        genome = "hg38", type = "coverage",
                                        background.title = "orangered",
                                        fill.coverage = "orangered",
                                        col.coverage = "orangered",
                                        transformation=function(x) {x * 2},
                                        ylim = c(0, ylimit),
                                        cex.axis = 1,
                                        col.title="black",
                                        col.axis="black",
                                        alpha = 0.9,
                                        name = "in Het",
                                        fontsize = 10,
                                        show.title = FALSE)
  Het_G_track_region <- AlignmentsTrack("VPS45_isogenic/Q20/down_45M/temp/CD_Het.alt.G.bam",
                                        isPaired = T, coverageOnly = T,
                                        chromosome = as.character(unique(seqnames(cell_type))),
                                        genome = "hg38", type = "coverage",
                                        background.title = "black",
                                        fill.coverage = "blue",
                                        col.coverage = "blue",
                                        transformation=function(x) {x * 2},
                                        ylim = c(0, ylimit),
                                        cex.axis = 1,
                                        col.title="black",
                                        col.axis="black",
                                        alpha = 0.9,
                                        name = "in Het",
                                        fontsize = 10,
                                        show.title = FALSE)
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

  Het_OvlTrack <- OverlayTrack(trackList = list(Het_AG_track_region,
                                                Het_G_track_region,
                                                Het_A_track_region),
                               name = "Pileup-Het",
                               fontsize = 10)

  Gross_OvlTrack <- OverlayTrack(trackList = list(GG_track_region,
                                                  # AG_track_region,
                                                  AA_track_region),
                                 name = "Pileup-Genotype",
                                 fontsize = 10)
  ############
  htTrack <- HighlightTrack(trackList = list(Het_OvlTrack, Gross_OvlTrack),
                            start = c(150067621),
                            width = c(1),
                            col = "black",
                            alpha = 0.7,
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
                  # Het_OvlTrack, Gross_OvlTrack,
                  snpTrack),
             sizes = c(0.5,0.5,4,4,0.5),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,


}
