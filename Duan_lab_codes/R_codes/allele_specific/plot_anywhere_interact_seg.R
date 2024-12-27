# function for plotting with interactionTrack
# Siwei 25 Dec 2019

plot_anywhere_interact_seg <- function(chr, start, end, gene_name = "",
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
                            # scale = 0.1
  )
  dtTrack_CN <- DataTrack(range = "~/3TB/merged_BAM/CN_RERE.bigWig",
                          genome = "hg38",
                          chromosome = "chr1",
                          name = "ATAC-Seq",
                          type = "l",
                          # fontcolor.group = "black", fontcolor.item = "black",
                          fontsize = 16,
                          col.title = "black",
                          col.axis = "black",
                          col = "darkblue",
                          transformation = function(x) {
                            x <- x / 10
                            x
                          },
                          ylim = c(0, 500))
  displayPars(dtTrack_CN) <- list(fontcolor.title = "black",
                                  fontsize.group = 18,
                                  baseline = 0#,
                                  # type = "mountain"
  )
  # alTrack_CN <- AlignmentsTrack("~/3TB/merged_BAM/CN_RERE_sorted.bam",
  #                               isPaired = T, coverageOnly = T,
  #                               chromosome = as.character(unique(seqnames(cell_type))),
  #                               genome = "hg38", type = "coverage",
  #                               transformation=function(x) {x * 4.078446},
  #                               background.title = "darkblue",
  #                               fill.coverage = "darkblue",
  #                               col.coverage = "darkblue",
  #                               ylim = c(0, ylimit),
  #                               cex.axis = 1.2,
  #                               name = "iN-Glut")
  # alTrack_NPC <- AlignmentsTrack("~/3TB/merged_BAM/NSC_all_sorted.bam",
  #                                isPaired = T, coverageOnly = T,
  #                                chromosome = as.character(unique(seqnames(cell_type))),
  #                                genome = "hg38", type = "coverage",
  #                                transformation=function(x) {x * 4.21563},
  #                                background.title = "goldenrod3",
  #                                fill.coverage = "goldenrod3",
  #                                col.coverage = "goldenrod3",
  #                                ylim = c(0, ylimit),
  #                                cex.axis = 1.2,
  #                                name = "NPC")
  # alTrack_DN <- AlignmentsTrack("~/3TB/merged_BAM/DN_all_sorted.bam",
  #                               isPaired = T, coverageOnly = T,
  #                               chromosome = as.character(unique(seqnames(cell_type))),
  #                               genome = "hg38", type = "coverage",
  #                               transformation=function(x) {x * 9.182318},
  #                               ylim = c(0, ylimit),
  #                               background.title = "darkorchid2",
  #                               fill.coverage = "darkorchid2",
  #                               col.coverage = "darkorchid2",
  #                               cex.axis = 1.2,
  #                               name = "iN-DN")
  # alTrack_GA <- AlignmentsTrack("~/3TB/merged_BAM/GA_all_sorted.bam",
  #                               isPaired = T, coverageOnly = T,
  #                               chromosome = as.character(unique(seqnames(cell_type))),
  #                               genome = "hg38", type = "coverage",
  #                               transformation=function(x) {x * 5.278966},
  #                               ylim = c(0, ylimit),
  #                               background.title = "magenta",
  #                               fill.coverage = "magenta",
  #                               col.coverage = "magenta",
  #                               cex.axis = 1.2,
  #                               name = "iN-GA")
  # alTrack_iPS <- AlignmentsTrack("~/3TB/merged_BAM/ips_all_sorted.bam",
  #                                isPaired = T, coverageOnly = T,
  #                                chromosome = as.character(unique(seqnames(cell_type))),
  #                                genome = "hg38", type = "coverage",
  #                                transformation=function(x) {x * 4.560793},
  #                                ylim = c(0, ylimit),
  #                                background.title = "brown",
  #                                fill.coverage = "brown",
  #                                col.coverage = "brown",
  #                                cex.axis = 1.2,
  #                                name = "iPS")
  dtTrack_RERE_eQTL <- DataTrack(gRanges_RERE_eQTL,
                                 name = "eQTL",
                                 fontcolor.group = "black", fontcolor.item = "black",
                                 # fontsize = 16,
                                 col.title = "black",
                                 col.axis = "black")

  interaction_track <- InteractionTrack(hic_data,
                                        name = "HiC",
                                        chromosome = chr)
  displayPars(interaction_track) <- list(fontcolor.title = "black",
                                         fontsize.group = 18,
                                         # fontcolor.item = "black",
                                         plot.trans = F)

  snpTrack <- AnnotationTrack(start = SNPposition, end = SNPposition, chromosome = chr,
                              id = SNPname, shape = "box",
                              name = SNPname, strand = "*",
                              group = "SNP",
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

  htTrack <- HighlightTrack(trackList = list(dtTrack_CN,
                                             interaction_track),
                            start = c(8408312, 8698883),
                            width = c(1, 1673),
                            col = "limegreen",
                            chromosome = as.character(unique(seqnames(cell_type))))



  plotTracks(list(
                  # gTrack,
                  # dtTrack_CN
                  # # alTrack_CN,
                  # # alTrack_NPC, #alTrack_DN, alTrack_GA, alTrack_iPS,
                  htTrack,
                  # interaction_track
                  # dtTrack_RERE_eQTL,
                  snpTrack
                  # grTrack
                  ),
             sizes = c(1,
               0.5,
               0.2
               # 0.5
               # 0.5, 0.5
               # 1,
               # 0.5, #1,1,1,1,
               # 0.5,
               # 0.2, 0.5
               ),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,

}
