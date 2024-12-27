# function for plotting rs16966337 site
# revised from plot_anywhere
# Use OverlayTrack to combine data tracks

# chr16:9939960

# init #####
library(Gviz)

library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ensembldb)
library(org.Hs.eg.db)

library(grDevices)
library(gridExtra)

library(RColorBrewer)

library(readr)
library(stringr)
#
options(ucscChromosomeNames = F)
options(Gviz.ucscUrl = "https://genome.ucsc.edu/cgi-bin/hgGateway")

plot_AsoC_peaks <- function(chr, start, end, gene_name = "",
                            alTrackRegion,
                            alTrackRefT,
                            alTrackAltC,
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

  # iTrack <- IdeogramTrack(genome = genome(cell_type),
  #                         chromosome = as.character(unique(seqnames(cell_type))),
  #                         fontcolor = "black",
  #                         fontsize = 18)

  gTrack <- GenomeAxisTrack(col = "black",
                            fontcolor = "black",
                            fontsize = 16,
                            scale = 0.1)

  alTrack_Region <-
    AlignmentsTrack(alTrackRegion,
                    # "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337/subsampled_0hr_GABA.bam",
                    isPaired = F, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    background.title = "palegreen",
                    fill.coverage = "palegreen",
                    col.coverage = "palegreen",
                    ylim = c(0, ylimit),
                    cex.axis = 1,
                    col.title = "black",
                    col.axis = "black",
                    # transformation = function(x) {x / 15},
                    name = title_name,
                    fontsize = 14,
                    show.title = FALSE)

  alTrack_T <-
    AlignmentsTrack(alTrackRefT,
                    # "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337/SNP_split_output/subsampled_0hr_GABA.ref.T.bam",
                    isPaired = F, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    background.title = "darkblue",
                    fill.coverage = "darkblue",
                    col.coverage = "darkblue",
                    ylim = c(0, ylimit),
                    cex.axis = 1,
                    # transformation = function(x) {x / 15},
                    # col.title="black",
                    name = "Allele_T")

  alTrack_C <- AlignmentsTrack(alTrackAltC,
                               # "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337/SNP_split_output/subsampled_0hr_GABA.alt.C.bam",
                               isPaired = F, coverageOnly = T,
                               chromosome = as.character(unique(seqnames(cell_type))),
                               genome = "hg38", type = "coverage",
                               background.title = "darkred",
                               fill.coverage = "darkred",
                               col.coverage = "darkred",
                               ylim = c(0, ylimit),
                               cex.axis = 1,
                               # transformation = function(x) {x / 15},
                               # col.title="black",
                               name = "Allele_C")

  OvlTrack <- OverlayTrack(trackList = list(alTrack_Region,
                                            alTrack_C,
                                            alTrack_T),
                           name = "REF_ALT_REGION",
                           fontsize = 14)

  snpTrack <- AnnotationTrack(start = 9788096,
                              end = 9788096, chromosome = "chr16",
                              id = "rs16966337", shape = "box",
                              name = "SNP", strand = "*",
                              group = c("rs16966337"),
                              fontcolor.group = "black", fontcolor.item = "black",
                              fontsize = 14,
                              col = "black", col.title = "black",
                              just.group = "below",
                              showID = TRUE,
                              cex.group = 1,
                              groupAnnotation = "id")

  # ucscGenes <-
  #   UcscTrack(genome = genome(cell_type),
  #             track = "NCBI RefSeq",
  #             table = "ncbiRefSeqCurated",
  #             trackType = "GeneRegionTrack",
  #             chromosome = as.character(unique(seqnames(cell_type))),
  #             rstarts = "exonStarts", rends = "exonEnds",
  #             gene = "name", symbol = 'name', transcript = "name",
  #             strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)
  # z <- ranges(ucscGenes)
  # mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "",
  #                                                            mcols(z)$symbol), "SYMBOL","REFSEQ"))
  # grTrack <- ucscGenes
  # ranges(grTrack) <- z
  # grTrack@dp@pars$col.line <- "black"
  # grTrack@dp@pars$fontcolor <- "black"
  # grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
  # grTrack@dp@pars$fontcolor.title <- "black"
  # grTrack@dp@pars$fontcolor.item <- "black"
  # grTrack@dp@pars$fontcolor.group <- "black"
  # grTrack@dp@pars$fontsize.group <- 18
  ########### plotting


  plotTracks(list(gTrack,
                  OvlTrack,
                  snpTrack),
             sizes = c(0.5,
                       1,
                       0.5),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,

}

TrackRegionBam <-
  list.files(path = "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337",
             pattern = "*.bam$",
             full.names = T)
TrackRefTBam <-
  list.files(path = "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337/SNP_split_output",
             pattern = "*ref.T.bam$",
             full.names = T)
TrackAltCBam <-
  list.files(path = "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337/SNP_split_output",
             pattern = "*alt.C.bam$",
             full.names = T)
BAM_names <-
  list.files(path = "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337",
             pattern = "*.bam$",
             full.names = F)
BAM_names <-
  str_split(string = BAM_names,
            pattern = "led_",
            simplify = T)[, 2]
BAM_names <-
  str_split(string = BAM_names,
            pattern = ".bam",
            simplify = T)[, 1]


# chr16:9788096
for (i in 1:9) {
  print(i)
  pdf(file = paste0("ASHG_2023_",
                    BAM_names[i]),
      width = 5,
      height = 2,
      title = BAM_names[i])
  plot_AsoC_peaks(chr = "chr16",
                  start = 9787096,
                  end = 9789096,
                  alTrackRegion = TrackRegionBam[i],
                  alTrackRefT = TrackRefTBam[i],
                  alTrackAltC = TrackAltCBam[i],
                  ylimit = 100,
                  title_name = BAM_names[i])
  dev.off()
}


plot_AsoC_peaks <- function(chr, start, end, gene_name = "",
                            alTrackRegion,
                            alTrackRefT,
                            alTrackAltC,
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

  # iTrack <- IdeogramTrack(genome = genome(cell_type),
  #                         chromosome = as.character(unique(seqnames(cell_type))),
  #                         fontcolor = "black",
  #                         fontsize = 18)

  gTrack <- GenomeAxisTrack(col = "black",
                            fontcolor = "black",
                            fontsize = 16,
                            scale = 0.1)

  alTrack_Region <-
    AlignmentsTrack(alTrackRegion,
                    # "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337/subsampled_0hr_GABA.bam",
                    isPaired = F, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    background.title = "palegreen",
                    fill.coverage = "palegreen",
                    col.coverage = "palegreen",
                    ylim = c(0, ylimit),
                    cex.axis = 1,
                    col.title = "black",
                    col.axis = "black",
                    # transformation = function(x) {x / 15},
                    name = title_name,
                    fontsize = 14,
                    show.title = FALSE)

  alTrack_T <-
    AlignmentsTrack(alTrackRefT,
                    # "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337/SNP_split_output/subsampled_0hr_GABA.ref.T.bam",
                    isPaired = F, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    background.title = "darkblue",
                    fill.coverage = "darkblue",
                    col.coverage = "darkblue",
                    ylim = c(0, ylimit),
                    cex.axis = 1,
                    # transformation = function(x) {x / 15},
                    # col.title="black",
                    name = "Allele_T")

  alTrack_C <- AlignmentsTrack(alTrackAltC,
                               # "/home/zhangs3/NVME/Bulk_ATACseq_ASoC/Gviz_code/plot_rs16966337/SNP_split_output/subsampled_0hr_GABA.alt.C.bam",
                               isPaired = F, coverageOnly = T,
                               chromosome = as.character(unique(seqnames(cell_type))),
                               genome = "hg38", type = "coverage",
                               background.title = "darkred",
                               fill.coverage = "darkred",
                               col.coverage = "darkred",
                               ylim = c(0, ylimit),
                               cex.axis = 1,
                               # transformation = function(x) {x / 15},
                               # col.title="black",
                               name = "Allele_C")

  OvlTrack <- OverlayTrack(trackList = list(alTrack_Region,
                                            alTrack_T,
                                            alTrack_C),
                           name = "REF_ALT_REGION",
                           fontsize = 14)

  snpTrack <- AnnotationTrack(start = 9788096,
                              end = 9788096, chromosome = "chr16",
                              id = "rs16966337", shape = "box",
                              name = "SNP", strand = "*",
                              group = c("rs16966337"),
                              fontcolor.group = "black", fontcolor.item = "black",
                              fontsize = 14,
                              col = "black", col.title = "black",
                              just.group = "below",
                              showID = TRUE,
                              cex.group = 1,
                              groupAnnotation = "id")

  # ucscGenes <-
  #   UcscTrack(genome = genome(cell_type),
  #             track = "NCBI RefSeq",
  #             table = "ncbiRefSeqCurated",
  #             trackType = "GeneRegionTrack",
  #             chromosome = as.character(unique(seqnames(cell_type))),
  #             rstarts = "exonStarts", rends = "exonEnds",
  #             gene = "name", symbol = 'name', transcript = "name",
  #             strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)
  # z <- ranges(ucscGenes)
  # mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "",
  #                                                            mcols(z)$symbol), "SYMBOL","REFSEQ"))
  # grTrack <- ucscGenes
  # ranges(grTrack) <- z
  # grTrack@dp@pars$col.line <- "black"
  # grTrack@dp@pars$fontcolor <- "black"
  # grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
  # grTrack@dp@pars$fontcolor.title <- "black"
  # grTrack@dp@pars$fontcolor.item <- "black"
  # grTrack@dp@pars$fontcolor.group <- "black"
  # grTrack@dp@pars$fontsize.group <- 18
  ########### plotting


  plotTracks(list(gTrack,
                  OvlTrack,
                  snpTrack),
             sizes = c(0.5,
                       1,
                       0.5),
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,

}

i <- 7

print(paste(i,
            BAM_names[i]))
pdf(file = paste0("ASHG_2023_",
                  BAM_names[i],
                  ".pdf"),
    width = 5,
    height = 2,
    title = BAM_names[i])
plot_AsoC_peaks(chr = "chr16",
                start = 9787096,
                end = 9789096,
                alTrackRegion = TrackRegionBam[i],
                alTrackRefT = TrackRefTBam[i],
                alTrackAltC = TrackAltCBam[i],
                ylimit = 100,
                title_name = BAM_names[i])
dev.off()


####
UcscTrack()


# chr16:9788096
ucscGenes <-
  UcscTrack(genome = "hg38",
            track = "NCBI RefSeq",
            table = "ncbiRefSeqCurated",
            trackType = "GeneRegionTrack",
            chromosome = "chr16",
            rstarts = "exonStarts",
            rends = "exonEnds",
            gene = "name",
            symbol = 'name',
            transcript = "name",
            strand = "strand",
            stacking = 'pack', showID = T, geneSymbol = T)
z <- ranges(ucscGenes)
mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "",
                                                           mcols(z)$symbol), "SYMBOL","REFSEQ"))
grTrack <- ucscGenes
ranges(grTrack) <- z
grTrack@dp@pars$col.line <- "black"
grTrack@dp@pars$fontcolor <- "black"
grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
grTrack@dp@pars$fontcolor.title <- "black"
grTrack@dp@pars$fontcolor.item <- "black"
grTrack@dp@pars$fontcolor.group <- "black"
grTrack@dp@pars$fontsize.group <- 18

test_UCSC_session <-
  browserSession("ucsc", url = "http://genome.ucsc.edu/cgi-bin")

