# function for plotting rs10792832 site
# revised from plot_anywhere
# Use OverlayTrack to combine data tracks

plot_AsoC_peaks <- function(chr, start, end, gene_name = "",
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

  gTrack <- GenomeAxisTrack(col = "black",
                            fontcolor = "black",
                            fontsize = 16,
                            scale = 0.1)

  alTrack_Region <- AlignmentsTrack("~/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/downsample_100M/MG/rs10792832_1MB/rs10792832_1MB_plot_split_AandG_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                background.title = "palegreen",
                                fill.coverage = "palegreen",
                                col.coverage = "palegreen",
                                ylim = c(0, ylimit),
                                cex.axis = 1,
                                col.title = "black",
                                col.axis = "black",
                                transformation = function(x) {x / 15},
                                name = title_name,
                                fontsize = 14,
                                show.title = FALSE)

  alTrack_A <- AlignmentsTrack("~/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/downsample_100M/MG/rs10792832_1MB/MG_A_merged_sorted.bam",
                               isPaired = T, coverageOnly = T,
                               chromosome = as.character(unique(seqnames(cell_type))),
                               genome = "hg38", type = "coverage",
                               background.title = "darkred",
                               fill.coverage = "darkred",
                               col.coverage = "darkred",
                               ylim = c(0, ylimit),
                               cex.axis = 1,
                               transformation = function(x) {x / 15},
                               # col.title="black",
                               name = "Allele_A")

  alTrack_G <- AlignmentsTrack("~/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/downsample_100M/MG/rs10792832_1MB/MG_G_merged_sorted.bam",
                              isPaired = T, coverageOnly = T,
                              chromosome = as.character(unique(seqnames(cell_type))),
                              genome = "hg38", type = "coverage",
                              background.title = "darkblue",
                              fill.coverage = "darkblue",
                              col.coverage = "darkblue",
                              ylim = c(0, ylimit),
                              cex.axis = 1,
                              transformation = function(x) {x / 15},
                              # col.title="black",
                              name = "Allele_G")

  OvlTrack <- OverlayTrack(trackList = list(alTrack_Region,
                                            alTrack_A,
                                            alTrack_G),
                           name = "REF_ALT_REGION",
                           fontsize = 14)

  snpTrack <- AnnotationTrack(start = 86156833,
                              end = 86156833, chromosome = "chr11",
                              id = "rs10792832", shape = "box",
                              name = "SNP", strand = "*",
                              group = c("rs10792832"),
                              fontcolor.group = "black", fontcolor.item = "black",
                              fontsize = 14,
                              col = "black", col.title = "black",
                              just.group = "below",
                              showID = TRUE,
                              cex.group = 1,
                              groupAnnotation = "id")
  ########### plotting


  plotTracks(list(gTrack, OvlTrack, snpTrack),
             sizes = c(0.5, 4, 1),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,

}
