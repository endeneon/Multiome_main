# function for plotting BCL11B site
# SNPname = "rs12895055", chr = 14, SNPposition = 99246457,
# revised from plot_anywhere
# Use OverlayTrack to combine data tracks

plot_AsoC_composite_BCL11B <- function(chr, start, end, gene_name = "",
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

  CN_alTrack_Region <- AlignmentsTrack("~/3TB/merged_BAM/CN_all_sorted.bam",
                                       isPaired = T, coverageOnly = T,
                                       chromosome = as.character(unique(seqnames(cell_type))),
                                       genome = "hg38", type = "coverage",
                                       background.title = "darkolivegreen4",
                                       fill.coverage = "darkolivegreen4",
                                       col.coverage = "darkolivegreen4",
                                       transformation=function(x) {x * 4.078446},
                                       ylim = c(0, ylimit),
                                       cex.axis = 1,
                                       col.title="black",
                                       col.axis="black",
                                       alpha = 0.9,
                                       name = "Glut",
                                       fontsize = 14,
                                       show.title = FALSE)

  CN_alTrack_A <- AlignmentsTrack(paste("BCL11B_4_het_bams/", "CN", "_all_sorted.ref.C.bam",
                                        collapse = "", sep = ""),
                                  isPaired = T, coverageOnly = T,
                                  chromosome = as.character(unique(seqnames(cell_type))),
                                  genome = "hg38", type = "coverage",
                                  background.title = "darkred",
                                  fill.coverage = "darkred",
                                  col.coverage = "darkred",
                                  transformation=function(x) {x * 4.078446},
                                  ylim = c(0, ylimit),
                                  cex.axis = 1,
                                  alpha = 0.9,
                                  col.title="black",
                                  col.axis="black",
                                  name = "Allele_C")

  CN_alTrack_G <- AlignmentsTrack(paste("BCL11B_4_het_bams/", "CN", "_all_sorted.alt.T.bam",
                                        collapse = "", sep = ""),
                                  isPaired = T, coverageOnly = T,
                                  chromosome = as.character(unique(seqnames(cell_type))),
                                  genome = "hg38", type = "coverage",
                                  background.title = "darkblue",
                                  fill.coverage = "darkblue",
                                  col.coverage = "darkblue",
                                  transformation=function(x) {x * 4.078446},
                                  ylim = c(0, ylimit),
                                  alpha = 0.9,
                                  cex.axis = 1,
                                  col.title="black",
                                  name = "Allele_T")

  CN_OvlTrack <- OverlayTrack(trackList = list(CN_alTrack_Region, CN_alTrack_A, CN_alTrack_G),
                              name = "REF_ALT_REGION",
                              fontsize = 14)

  NPC_alTrack_Region <- AlignmentsTrack("~/3TB/merged_BAM/NSC_all_sorted.bam",
                                        isPaired = T, coverageOnly = T,
                                        chromosome = as.character(unique(seqnames(cell_type))),
                                        genome = "hg38", type = "coverage",
                                        background.title = "darkolivegreen4",
                                        fill.coverage = "darkolivegreen4",
                                        col.coverage = "darkolivegreen4",
                                        transformation=function(x) {x * 4.21563},
                                        ylim = c(0, ylimit),
                                        cex.axis = 1,
                                        col.title="black",
                                        col.axis="black",
                                        alpha = 0.9,
                                        name = "NPC",
                                        fontsize = 14,
                                        show.title = FALSE)

  NPC_alTrack_A <- AlignmentsTrack(paste("BCL11B_4_het_bams/", "NSC", "_all_sorted.ref.C.bam",
                                         collapse = "", sep = ""),
                                   isPaired = T, coverageOnly = T,
                                   chromosome = as.character(unique(seqnames(cell_type))),
                                   genome = "hg38", type = "coverage",
                                   background.title = "darkred",
                                   fill.coverage = "darkred",
                                   col.coverage = "darkred",
                                   transformation=function(x) {x * 4.21563},
                                   ylim = c(0, ylimit),
                                   alpha = 0.9,
                                   cex.axis = 1,
                                   col.title = "black",
                                   col.axis = "black",
                                   name = "Allele_C")

  NPC_alTrack_G <- AlignmentsTrack(paste("BCL11B_4_het_bams/", "NSC", "_all_sorted.alt.T.bam",
                                         collapse = "", sep = ""),
                                   isPaired = T, coverageOnly = T,
                                   chromosome = as.character(unique(seqnames(cell_type))),
                                   genome = "hg38", type = "coverage",
                                   background.title = "darkblue",
                                   fill.coverage = "darkblue",
                                   col.coverage = "darkblue",
                                   transformation=function(x) {x * 4.21563},
                                   alpha = 0.9,
                                   ylim = c(0, ylimit),
                                   cex.axis = 1,
                                   col.title="black",
                                   col.axis = "black",
                                   name = "Allele_T")

  NPC_OvlTrack <- OverlayTrack(trackList = list(NPC_alTrack_G, NPC_alTrack_A, NPC_alTrack_Region),
                               name = "REF_ALT_REGION",
                               fontsize = 14)

  snpTrack <- AnnotationTrack(start = 99246457, end = 99246457, chromosome = "chr14",
                              id = "rs12895055", shape = "box",
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
