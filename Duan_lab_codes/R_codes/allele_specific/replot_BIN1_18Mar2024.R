# Siwei 09 Aug 2023
# Siwei 05 Jul 2023
# plot 1MB proximal region of rs1532278 (CLU)
# chr8:27608798

# init #####
{
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
  library(RCurl)
  library(stringr)

  library(EnsDb.Hsapiens.v86)
}

{
  options(ucscChromosomeNames = F)
  # !! use HTTPS by redesignate URL here!! ####
  options(Gviz.ucscUrl = "https://genome-euro.ucsc.edu/cgi-bin/")

}

ens <-
  EnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86@ensdb)

EnsTrack <-
  getGeneRegionTrackForGviz(ens)
# EnsTrack <-
#   getGeneRegionTrackForGviz(ens,
#                             filter = GRangesFilter(GRanges("8:27600000-27620000")))
EnsTrack <-
  GeneRegionTrack(EnsTrack)
EnsTrack@range@seqnames@values <-
  as.factor(paste0("chr",
                   EnsTrack@range@seqnames@values))
EnsTrack@range@seqinfo@seqnames <-
  paste0("chr",
         EnsTrack@range@seqinfo@seqnames)

EnsTrack@stacking <- "dense"


df_GWASTrack <-
  read_delim("~/backuped_space/Siwei_python_workspaces/AD_plot/data/Alz_GWAS_track_hg38_subsetted.scatter.cons.txt",
             delim = " ", escape_double = FALSE,
             col_names = FALSE, trim_ws = TRUE)
colnames(df_GWASTrack) <-
  c("chr", "start", "end", "score")
df_GWASTrack <-
  df_GWASTrack[df_GWASTrack$chr == "chr11", ]

# df_GWASTrack$chr
# df_GWASTrack$score <-
#   0 - log10(df_GWASTrack$score)

# df_GWASTrack <-
#   df_GWASTrack[df_GWASTrack$chr == "chr2", ]

gRangesGWASTrack <-
  makeGRangesFromDataFrame(df = df_GWASTrack,
                           keep.extra.columns = T,
                           ignore.strand = T)


df_ASoCTrack <-
  read_delim("~/backuped_space/Siwei_python_workspaces/AD_plot/data/MG_all_SNP_4_circos.txt",
             delim = " ", escape_double = FALSE,
             col_names = FALSE, trim_ws = TRUE)
colnames(df_ASoCTrack) <-
  c("chr", "start", "end", "score")
df_ASoCTrack <-
  df_ASoCTrack[df_ASoCTrack$chr == "chr11", ]

gRangesASoCTrack <-
  makeGRangesFromDataFrame(df = df_ASoCTrack,
                           keep.extra.columns = T,
                           ignore.strand = T)
# Sys.getenv("http_proxy")
# curl_handle <-
#   getCurlHandle()
# curlSetOpt(.opts = list(proxy = 'http://proxy.enh.org:8080'),
#            curl = curl_handle)
# ans <-
#   getURL("http://www.google.com",
#          curl = curl_handle)
# #

# getURL("https://genome-euro.ucsc.edu/cgi-bin/")
# load data


# function for plotting
plot_5_types <-
  function(chr, start, end, gene_name = "",
           SNPname = "", SNPposition = 1L,
           mcols = 100, strand = "+",
           GWASTrack = "",
           lineWidth = 1,
           minHeight = 10,
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

  iTrack <-
    IdeogramTrack(genome = genome(cell_type),
                          chromosome = as.character(unique(seqnames(cell_type))),
                          fontcolor = "black",
                          fontsize = 18)
  gTrack <-
    GenomeAxisTrack(col = "black",
                            fontcolor = "black",
                            fontsize = 16,
                            scale = 0.1)
  # "GWAS_Alz_2022_Jansen_et_al_rs10792832_hg38.bedGraph"
  gwasTrack <-
    DataTrack(range = gRangesGWASTrack,
              genome = "hg38",
              type = "p",
              # chromosome = "chr11",
              name = "Alz\nGWAS",
              background.title = "darkblue",
              fill.coverage = "darkblue",
              col.coverage = "darkblue",
              cex.axis = 0.5,
              cex.title = 0.8)

  # "Alz_GWAS_Bellenguez_etal_Nat_Genet_2022")
  # alTrack_iMG <-
  #   AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/MG_downsampled_100M_het_rs10792832.bam",
  #                   isPaired = T, coverageOnly = T,
  #                   chromosome = as.character(unique(seqnames(cell_type))),
  #                   genome = "hg38", type = "coverage",
  #                   transformation = function(x) {x * 2},
  #                   ylim = c(0, ylimit),
  #                   background.title = "tomato4",
  #                   fill.coverage = "tomato4",
  #                   col.coverage = "tomato4",
  #                   cex.axis = 1,
  #                   cex.title = 0.8,
  #                   lwd = lineWidth,
  #
  #                   name = "iMG")
  # alTrack_iAst <-
  #   AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/Ast_downsampled_100M_het_rs10792832.bam",
  #                   isPaired = T, coverageOnly = T,
  #                   chromosome = as.character(unique(seqnames(cell_type))),
  #                   genome = "hg38", type = "coverage",
  #                   transformation=function(x) {x * 2},
  #                   ylim = c(0, ylimit),
  #                   background.title = "#1B9E77",
  #                   fill.coverage = "#1B9E77",
  #                   col.coverage = "#1B9E77",
  #                   cex.axis = 1,
  #                   cex.title = 0.8,
  #                   lwd = lineWidth,
  #                   lwd.border = lineWidth,
  #                   min.width = lineWidth,
  #                   min.height = minHeight,
  #                   minCoverageHeight = minHeight,
  #                   name = "iAst")

  alTrack_iMG <-
    AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/MG_downsampled_100M_het_rs10792832.bam",
                    isPaired = T, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    transformation = function(x) {x * 1},
                    ylim = c(0, ylimit),
                    background.title = "brown",
                    fill.coverage = "brown",
                    col.coverage = "brown",
                    cex.axis = 0.5,
                    cex.title = 0.8,
                    lwd = lineWidth,

                    name = "iMG")
  alTrack_iAst <-
    AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/Ast_downsampled_100M_het_rs10792832.bam",
                    isPaired = T, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    transformation=function(x) {x * 1},
                    ylim = c(0, ylimit),
                    background.title = "deepskyblue4",
                    fill.coverage = "deepskyblue4",
                    col.coverage = "deepskyblue4",
                    cex.axis = 0.5,
                    cex.title = 0.8,
                    lwd = lineWidth,
                    lwd.border = lineWidth,
                    min.width = lineWidth,
                    min.height = minHeight,
                    minCoverageHeight = minHeight,
                    name = "iAst")
  alTrack_NGN2 <-
    AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/NGN2-glut_downsampled_100M_het_rs10792832.bam",
                    isPaired = T, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    transformation=function(x) {x * 1},
                    background.title = "#D95F02",
                    fill.coverage = "#D95F02",
                    col.coverage = "#D95F02",
                    ylim = c(0, ylimit),
                    cex.axis = 0.5,
                    cex.title = 0.8,
                    lwd = lineWidth,
                    lwd.border = lineWidth,
                    min.width = lineWidth,
                    min.height = minHeight,
                    minCoverageHeight = minHeight,
                    name = "NGN2\nGlut")

  alTrack_DN <-
    AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/DN_downsampled_100M_het_rs10792832.bam",
                    isPaired = T, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    transformation = function(x) {x * 1},
                    ylim = c(0, ylimit),
                    background.title = "#7570B3",
                    fill.coverage = "#7570B3",
                    col.coverage = "#7570B3",
                    cex.axis = 0.5,
                    cex.title = 0.8,
                    lwd = lineWidth,
                    min.height = minHeight,
                    minCoverageHeight = minHeight,
                    name = "DN")
  alTrack_GA <-
    AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/GA_downsampled_100M_het_rs10792832.bam",
                    isPaired = T, coverageOnly = T,
                    chromosome = as.character(unique(seqnames(cell_type))),
                    genome = "hg38", type = "coverage",
                    transformation = function(x) {x * 1},
                    ylim = c(0, ylimit),
                    background.title = "#E7298A",
                    fill.coverage = "#E7298A",
                    col.coverage = "#E7298A",
                    cex.axis = 0.5,
                    cex.title = 0.8,
                    lwd = lineWidth,
                    lwd.border = lineWidth,
                    min.width = lineWidth,
                    min.height = minHeight,
                    minCoverageHeight = minHeight,
                    name = "GA")

  snpTrack <-
    AnnotationTrack(start = SNPposition, end = SNPposition, chromosome = chr,
                    id = SNPname, shape = "box",
                    name = "SNP", strand = "*",
                    group = SNPname,
                    fontcolor.group = "black", fontcolor.item = "black",
                    fontsize = 18,
                    col = "black", col.title = "black",
                    just.group = "below",
                    showID = TRUE,
                    cex.group = 0.8,
                    cex.title = 0.6,
                    rotate.title = 0,
                    groupAnnotation = "id")

  ########## plotting
  ucscGenes <-
    UcscTrack(genome = genome(cell_type),
              track = "NCBI RefSeq",
              table = "ncbiRefSeqCurated",
              trackType = "GeneRegionTrack",
              chromosome = as.character(unique(seqnames(cell_type))),
              rstarts = "exonStarts", rends = "exonEnds",
              gene = "name", symbol = 'name', transcript = "name",
              strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)
  z <- ranges(ucscGenes)
  mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "",
                                                             mcols(z)$symbol), "SYMBOL","REFSEQ"))
  # EnsTrack <-
  #   getGeneRegionTrackForGviz(ens,
  #                             filter = GRangesFilter(GRanges(paste0(str_remove_all(cell_type@ranges@NAMES,
  #                                                                                  pattern = "chr"),
  #                                                                   ":",
  #                                                                   as.character(cell_type@ranges@start - x_offset_1),
  #                                                                   "-",
  #                                                                   as.character(cell_type@ranges@start + cell_type@ranges@width + x_offset_2)))))
  # EnsTrack <-
  #   getGeneRegionTrackForGviz(ens)

  # grTrack <- GeneRegionTrack(EnsTrack)
  # grTrack <- EnsTrack
    grTrack <- ucscGenes
  ranges(grTrack) <- z
  grTrack@dp@pars$col.line <- "black"
  grTrack@dp@pars$fontcolor <- "black"
  grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
  grTrack@name <- "Gene"
  grTrack@dp@pars$fontcolor.title <- "black"
  grTrack@dp@pars$fontcolor.item <- "black"
  grTrack@dp@pars$fontcolor.group <- "black"
  grTrack@dp@pars$fontsize.group <- 16
  grTrack@dp@pars$cex.title = 0.8
  grTrack@dp@pars$rotate.title = 0

  ######

  # htTrack <- HighlightTrack(trackList = list(alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS),
  #                           start = c(26212000, 26241000),
  #                           width = c(15000, 2000)
  #                           chromosome = as.character(unique(seqnames(cell_type))))



  plotTracks(list(
    iTrack, gTrack,
    # gwasTrack,
    # alTrack_hMG, alTrack_hAst,
    alTrack_iMG, alTrack_iAst,
    alTrack_NGN2,
    alTrack_DN, alTrack_GA,
    # snpTrack,
    grTrack
  ),
  sizes = c(
    0.5, 0.5,
    # 1,
    1, 1,
    1,
    1, 1,
    # 0.5,
    0.5
  ),
  chromosome = cell_type@ranges@NAMES,
  from = (cell_type@ranges@start - x_offset_1),
  to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
  transcriptAnnotation = "gene_id",
  collapseTranscripts = "gene_id")#,

}


# make the plot #####
# BIN1 w/ SNP
plot_5_types(chr = "chr2",
               start = 127084000,
               end = 127142500,
               SNPposition = 127135234,
               SNPname = "rs6733839",
               # GWASTrack = gRangesGWASTrack,
               lineWidth = 0,
               minHeight = 0,
               ylimit = 400)

## TMEM119
plot_5_types(chr = "chr12",
             start = 108580000,
             end = 108610000,
             # SNPposition = 127135234,
             # SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 100)

## CD45 (PTPRC)
plot_5_types(chr = "chr1",
             start = 198600000,
             end = 198800000,
             # SNPposition = 127135234,
             # SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 100)

## CD11B (ITGAM)
plot_5_types(chr = "chr16",
             start = 31240000,
             end = 31380000,
             # SNPposition = 127135234,
             # SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 100)

## GFAP
plot_5_types(chr = "chr17",
             start = 44890000,
             end = 44918000,
             # SNPposition = 127135234,
             # SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 200)

brewer.pal(8, name = "Dark2")


## test BiomartGeneRegionTrack #####
BiomartGenes <-
  BiomartGeneRegionTrack(genome = genome(cell_type),
            track = "NCBI RefSeq",
            table = "ncbiRefSeqCurated",
            trackType = "GeneRegionTrack",
            chromosome = as.character(unique(seqnames(cell_type))),
            rstarts = "exonStarts", rends = "exonEnds",
            gene = "name", symbol = 'name', transcript = "name",
            strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)


library(EnsDb.Hsapiens.v86)
ens <-
  EnsDb(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86@ensdb)
EnsTrack <-
  getGeneRegionTrackForGviz(ens)
# EnsTrack <-
#   getGeneRegionTrackForGviz(ens,
#                             filter = GRangesFilter(GRanges("8:27600000-27620000")))
EnsTrack <-
  GeneRegionTrack(EnsTrack)
EnsTrack@range@seqnames@values <-
  as.factor(paste0("chr",
                   EnsTrack@range@seqnames@values))
EnsTrack@range@seqinfo@seqnames <-
  paste0("chr",
                   EnsTrack@range@seqinfo@seqnames)
# EnsTrack@chromosome <-
#   paste0("chr",
#          EnsTrack@chromosome)
# EnsTrack@genome <-
#   paste0("chr",
#          EnsTrack@genome)
EnsTrack@stacking <- "dense"

plotTracks(list(EnsTrack),
           chromosome = "chr8",
           from = 27600000,
           to = 27620000,
           collapseTranscripts = "symbol",
           transcriptAnnotation = "symbol")

bmTrack <-
  BiomartGeneRegionTrack(
  start = 100, end = 1660,
  chromosome = "chr1", genome = "hg38")

## Use UCSC tracks #####
# function for plotting
plot_5_types <-
  function(chr, start, end, gene_name = "",
           SNPname = "", SNPposition = 1L,
           mcols = 100, strand = "+",
           GWASTrack = "",
           lineWidth = 1,
           minHeight = 10,
           x_offset_1 = 0,
           x_offset_2 = 0,
           ylimit = 800,
           gene_track_height = 1) {
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

    iTrack <-
      IdeogramTrack(genome = genome(cell_type),
                    chromosome = as.character(unique(seqnames(cell_type))),
                    fontcolor = "black",
                    fontsize = 18)
    gTrack <-
      GenomeAxisTrack(col = "black",
                      fontcolor = "black",
                      fontsize = 16,
                      scale = 0.1)
    # "GWAS_Alz_2022_Jansen_et_al_rs10792832_hg38.bedGraph"
    gwasTrack <-
      DataTrack(range = gRangesGWASTrack,
                genome = "hg38",
                type = "p",
                # chromosome = "chr11",
                name = "Alz\nGWAS",
                background.title = "darkblue",
                fill.coverage = "darkblue",
                col.coverage = "darkblue",
                cex.axis = 0.5,
                cex.title = 0.8)

    # "Alz_GWAS_Bellenguez_etal_Nat_Genet_2022")
    # alTrack_iMG <-
    #   AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/MG_downsampled_100M_het_rs10792832.bam",
    #                   isPaired = T, coverageOnly = T,
    #                   chromosome = as.character(unique(seqnames(cell_type))),
    #                   genome = "hg38", type = "coverage",
    #                   transformation = function(x) {x * 2},
    #                   ylim = c(0, ylimit),
    #                   background.title = "tomato4",
    #                   fill.coverage = "tomato4",
    #                   col.coverage = "tomato4",
    #                   cex.axis = 1,
    #                   cex.title = 0.8,
    #                   lwd = lineWidth,
    #
    #                   name = "iMG")
    # alTrack_iAst <-
    #   AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/Ast_downsampled_100M_het_rs10792832.bam",
    #                   isPaired = T, coverageOnly = T,
    #                   chromosome = as.character(unique(seqnames(cell_type))),
    #                   genome = "hg38", type = "coverage",
    #                   transformation=function(x) {x * 2},
    #                   ylim = c(0, ylimit),
    #                   background.title = "#1B9E77",
    #                   fill.coverage = "#1B9E77",
    #                   col.coverage = "#1B9E77",
    #                   cex.axis = 1,
    #                   cex.title = 0.8,
    #                   lwd = lineWidth,
    #                   lwd.border = lineWidth,
    #                   min.width = lineWidth,
    #                   min.height = minHeight,
    #                   minCoverageHeight = minHeight,
    #                   name = "iAst")

    alTrack_iMG <-
      AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/MG_downsampled_100M_het_rs10792832.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation = function(x) {x * 1},
                      ylim = c(0, ylimit),
                      background.title = "brown",
                      fill.coverage = "brown",
                      col.coverage = "brown",
                      cex.axis = 0.5,
                      cex.title = 0.8,
                      lwd = lineWidth,

                      name = "iMG")
    alTrack_iAst <-
      AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/Ast_downsampled_100M_het_rs10792832.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation=function(x) {x * 1},
                      ylim = c(0, ylimit),
                      background.title = "deepskyblue4",
                      fill.coverage = "deepskyblue4",
                      col.coverage = "deepskyblue4",
                      cex.axis = 0.5,
                      cex.title = 0.8,
                      lwd = lineWidth,
                      lwd.border = lineWidth,
                      min.width = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "iAst")
    alTrack_NGN2 <-
      AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/NGN2-glut_downsampled_100M_het_rs10792832.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation=function(x) {x * 1},
                      background.title = "#D95F02",
                      fill.coverage = "#D95F02",
                      col.coverage = "#D95F02",
                      ylim = c(0, ylimit),
                      cex.axis = 0.5,
                      cex.title = 0.8,
                      lwd = lineWidth,
                      lwd.border = lineWidth,
                      min.width = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "NGN2\nGlut")

    alTrack_DN <-
      AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/DN_downsampled_100M_het_rs10792832.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation = function(x) {x * 1},
                      ylim = c(0, ylimit),
                      background.title = "#7570B3",
                      fill.coverage = "#7570B3",
                      col.coverage = "#7570B3",
                      cex.axis = 0.5,
                      cex.title = 0.8,
                      lwd = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "DN")
    alTrack_GA <-
      AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/GA_downsampled_100M_het_rs10792832.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation = function(x) {x * 1},
                      ylim = c(0, ylimit),
                      background.title = "#E7298A",
                      fill.coverage = "#E7298A",
                      col.coverage = "#E7298A",
                      cex.axis = 0.5,
                      cex.title = 0.8,
                      lwd = lineWidth,
                      lwd.border = lineWidth,
                      min.width = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "GA")

    snpTrack <-
      AnnotationTrack(start = SNPposition, end = SNPposition, chromosome = chr,
                      id = SNPname, shape = "box",
                      name = "SNP", strand = "*",
                      group = SNPname,
                      fontcolor.group = "black", fontcolor.item = "black",
                      fontsize = 18,
                      col = "black", col.title = "black",
                      just.group = "below",
                      showID = TRUE,
                      cex.group = 0.8,
                      cex.title = 0.6,
                      rotate.title = 0,
                      groupAnnotation = "id")

    ########### plotting
    ucscGenes <-
      UcscTrack(genome = genome(cell_type),
                track = "NCBI RefSeq",
                table = "ncbiRefSeqCurated",
                trackType = "GeneRegionTrack",
                chromosome = as.character(unique(seqnames(cell_type))),
                rstarts = "exonStarts", rends = "exonEnds",
                gene = "name", symbol = 'name', transcript = "name",
                strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)
    z <- ranges(ucscGenes)
    mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "",
                                                               mcols(z)$symbol), "SYMBOL","REFSEQ"))
    # EnsTrack <-
    #   getGeneRegionTrackForGviz(ens,
    #                             filter = GRangesFilter(GRanges(paste0(str_remove_all(cell_type@ranges@NAMES,
    #                                                                                  pattern = "chr"),
    #                                                                   ":",
    #                                                                   as.character(cell_type@ranges@start - x_offset_1),
    #                                                                   "-",
    #                                                                   as.character(cell_type@ranges@start + cell_type@ranges@width + x_offset_2)))))
    # EnsTrack <-
    #   getGeneRegionTrackForGviz(ens)

    # grTrack <- GeneRegionTrack(EnsTrack)
    # grTrack <- EnsTrack
    grTrack <- ucscGenes
    ranges(grTrack) <- z
    grTrack@dp@pars$col.line <- "black"
    grTrack@dp@pars$fontcolor <- "black"
    grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
    grTrack@name <- "Gene"
    grTrack@dp@pars$fontcolor.title <- "black"
    grTrack@dp@pars$fontcolor.item <- "black"
    grTrack@dp@pars$fontcolor.group <- "black"
    grTrack@dp@pars$fontsize.group <- 16
    grTrack@dp@pars$cex.title = 0.8
    grTrack@dp@pars$rotate.title = 0

    ######

    # htTrack <- HighlightTrack(trackList = list(alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS),
    #                           start = c(26212000, 26241000),
    #                           width = c(15000, 2000)
    #                           chromosome = as.character(unique(seqnames(cell_type))))



    plotTracks(list(
      iTrack, gTrack,
      # gwasTrack,
      # alTrack_hMG, alTrack_hAst,
      alTrack_iMG, alTrack_iAst,
      alTrack_NGN2,
      alTrack_DN, alTrack_GA,
      # snpTrack,
      grTrack
    ),
    sizes = c(
      0.5, 0.5,
      # 1,
      1, 1,
      1,
      1, 1,
      # 0.5,
      gene_track_height
    ),
    chromosome = cell_type@ranges@NAMES,
    from = (cell_type@ranges@start - x_offset_1),
    to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
    transcriptAnnotation = "transcript",
    collapseTranscripts = "transcript")#,

  }



## GFAP
plot_5_types(chr = "chr17",
             start = 44895000,
             end = 44928000,
             # SNPposition = 127135234,
             # SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 200)

## S100B
plot_5_types(chr = "chr21",
             start = 46507000,
             end = 46684000,
             # SNPposition = 127135234,
             # SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 100,
             gene_track_height = 0.5)

## SLC1A3 EAAT1
plot_5_types(chr = "chr5",
             start = 36600000,
             end = 36700000,
             # SNPposition = 127135234,
             # SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 200)

## SLC1A2 EAAT2
plot_5_types(chr = "chr11",
             start = 35200000,
             end = 35450000,
             # SNPposition = 127135234,
             # SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 200)

## IBA1 AIF1
plot_5_types(chr = "chr6",
             start = 31600000,
             end = 31650000,
             # SNPposition = 127135234,
             # SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 400)

## TREM2
plot_5_types(chr = "chr6",
             start = 41090000,
             end = 41200000,
             # SNPposition = 127135234,
             # SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 200)

## SPI1
plot_5_types(chr = "chr11",
             start = 47320000,
             end = 47420000,
             # SNPposition = 127135234,
             # SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 200)

## VIM
plot_5_types(chr = "chr10",
             start = 17110000,
             end = 17350000,
             # SNPposition = 127135234,
             # SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 400)

# BIN1 w/ SNP
plot_5_types(chr = "chr2",
             start = 127000000,
             end = 127202500,
             SNPposition = 127135234,
             SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 400)





## plot w GWAS and SNP tracks, UCSC genes #####
plot_5_types <-
  function(chr, start, end, gene_name = "",
           SNPname = "", SNPposition = 1L,
           mcols = 100, strand = "+",
           GWASTrack = "",
           ASoCTrack = "",
           lineWidth = 1,
           minHeight = 10,
           x_offset_1 = 0,
           x_offset_2 = 0,
           ylimit = 800,
           gene_track_height = 1) {
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

    iTrack <-
      IdeogramTrack(genome = genome(cell_type),
                    chromosome = as.character(unique(seqnames(cell_type))),
                    fontcolor = "black",
                    fontsize = 18)
    gTrack <-
      GenomeAxisTrack(col = "black",
                      fontcolor = "black",
                      fontsize = 16)
    # "GWAS_Alz_2022_Jansen_et_al_rs10792832_hg38.bedGraph"
    gwasTrack <-
      DataTrack(range = gRangesGWASTrack,
                genome = "hg38",
                type = "p",
                # chromosome = "chr11",
                name = "Alz\nGWAS",
                background.title = "darkblue",
                fill.coverage = "darkblue",
                col.coverage = "darkblue",
                cex = 0.5,
                cex.axis = 0.5,
                cex.title = 0.8)

    ASoC_track <-
      DataTrack(range = ASoCTrack,
                genome = "hg38",
                type = "p",
                # chromosome = "chr11",
                name = "ASoC",
                background.title = "seagreen4",
                fill.coverage = "seagreen4",
                col.coverage = "seagreen4",
                col = "forestgreen",
                ylim = c(0, 10),
                cex = 0.5,
                cex.axis = 0.5,
                cex.title = 0.8)

    # "Alz_GWAS_Bellenguez_etal_Nat_Genet_2022")
    # alTrack_iMG <-
    #   AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/MG_downsampled_100M_het_rs10792832.bam",
    #                   isPaired = T, coverageOnly = T,
    #                   chromosome = as.character(unique(seqnames(cell_type))),
    #                   genome = "hg38", type = "coverage",
    #                   transformation = function(x) {x * 2},
    #                   ylim = c(0, ylimit),
    #                   background.title = "tomato4",
    #                   fill.coverage = "tomato4",
    #                   col.coverage = "tomato4",
    #                   cex.axis = 1,
    #                   cex.title = 0.8,
    #                   lwd = lineWidth,
    #
    #                   name = "iMG")
    # alTrack_iAst <-
    #   AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/Ast_downsampled_100M_het_rs10792832.bam",
    #                   isPaired = T, coverageOnly = T,
    #                   chromosome = as.character(unique(seqnames(cell_type))),
    #                   genome = "hg38", type = "coverage",
    #                   transformation=function(x) {x * 2},
    #                   ylim = c(0, ylimit),
    #                   background.title = "#1B9E77",
    #                   fill.coverage = "#1B9E77",
    #                   col.coverage = "#1B9E77",
    #                   cex.axis = 1,
    #                   cex.title = 0.8,
    #                   lwd = lineWidth,
    #                   lwd.border = lineWidth,
    #                   min.width = lineWidth,
    #                   min.height = minHeight,
    #                   minCoverageHeight = minHeight,
    #                   name = "iAst")

    alTrack_iMG <-
      AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/MG_downsampled_100M_het_rs10792832.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation = function(x) {x * 1},
                      ylim = c(0, ylimit),
                      background.title = "brown",
                      fill.coverage = "brown",
                      col.coverage = "brown",
                      cex.axis = 0.5,
                      cex.title = 0.8,
                      lwd = lineWidth,

                      name = "iMG")
    alTrack_iAst <-
      AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/Ast_downsampled_100M_het_rs10792832.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation=function(x) {x * 1},
                      ylim = c(0, ylimit),
                      background.title = "deepskyblue4",
                      fill.coverage = "deepskyblue4",
                      col.coverage = "deepskyblue4",
                      cex.axis = 0.5,
                      cex.title = 0.8,
                      lwd = lineWidth,
                      lwd.border = lineWidth,
                      min.width = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "iAst")
    alTrack_NGN2 <-
      AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/NGN2-glut_downsampled_100M_het_rs10792832.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation=function(x) {x * 1},
                      background.title = "#D95F02",
                      fill.coverage = "#D95F02",
                      col.coverage = "#D95F02",
                      ylim = c(0, ylimit),
                      cex.axis = 0.5,
                      cex.title = 0.8,
                      lwd = lineWidth,
                      lwd.border = lineWidth,
                      min.width = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "iGlut")

    alTrack_DN <-
      AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/DN_downsampled_100M_het_rs10792832.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation = function(x) {x * 1},
                      ylim = c(0, ylimit),
                      background.title = "#7570B3",
                      fill.coverage = "#7570B3",
                      col.coverage = "#7570B3",
                      cex.axis = 0.5,
                      cex.title = 0.8,
                      lwd = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "iDN")
    alTrack_GA <-
      AlignmentsTrack("/home/zhangs3/Data/FASTQ/Unified_peak_count_Roussos_et_al_Nat_Genet_2022/merged_100M_BAMs/GA_downsampled_100M_het_rs10792832.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation = function(x) {x * 1},
                      ylim = c(0, ylimit),
                      background.title = "#E7298A",
                      fill.coverage = "#E7298A",
                      col.coverage = "#E7298A",
                      cex.axis = 0.5,
                      cex.title = 0.8,
                      lwd = lineWidth,
                      lwd.border = lineWidth,
                      min.width = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "iGA")

    snpTrack <-
      AnnotationTrack(start = SNPposition,
                      end = SNPposition,
                      chromosome = chr,
                      id = SNPname, shape = "box",
                      name = "SNP", strand = "*",
                      group = SNPname,
                      fontcolor.group = "black", fontcolor.item = "black",
                      fontsize = 18,
                      col = "black", col.title = "black",
                      just.group = "below",
                      showID = TRUE,
                      cex.group = 0.8,
                      cex.title = 0.6,
                      rotate.title = 0,
                      groupAnnotation = "id")

    ########### plotting
    ucscGenes <-
      UcscTrack(genome = genome(cell_type),
                track = "NCBI RefSeq",
                table = "ncbiRefSeqCurated",
                trackType = "GeneRegionTrack",
                chromosome = as.character(unique(seqnames(cell_type))),
                rstarts = "exonStarts", rends = "exonEnds",
                gene = "name", symbol = 'name', transcript = "name",
                strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)
    z <- ranges(ucscGenes)
    mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "",
                                                               mcols(z)$symbol), "SYMBOL","REFSEQ"))
    # EnsTrack <-
    #   getGeneRegionTrackForGviz(ens,
    #                             filter = GRangesFilter(GRanges(paste0(str_remove_all(cell_type@ranges@NAMES,
    #                                                                                  pattern = "chr"),
    #                                                                   ":",
    #                                                                   as.character(cell_type@ranges@start - x_offset_1),
    #                                                                   "-",
    #                                                                   as.character(cell_type@ranges@start + cell_type@ranges@width + x_offset_2)))))
    # EnsTrack <-
    #   getGeneRegionTrackForGviz(ens)

    # grTrack <- GeneRegionTrack(EnsTrack)
    # grTrack <- EnsTrack
    grTrack <- ucscGenes
    ranges(grTrack) <- z
    grTrack@dp@pars$col.line <- "black"
    grTrack@dp@pars$fontcolor <- "black"
    grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
    grTrack@name <- "Gene"
    grTrack@dp@pars$fontcolor.title <- "black"
    grTrack@dp@pars$fontcolor.item <- "black"
    grTrack@dp@pars$fontcolor.group <- "black"
    grTrack@dp@pars$fontsize.group <- 16
    grTrack@dp@pars$cex.title = 0.8
    grTrack@dp@pars$rotate.title = 0

    ######

    # htTrack <- HighlightTrack(trackList = list(alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS),
    #                           start = c(26212000, 26241000),
    #                           width = c(15000, 2000)
    #                           chromosome = as.character(unique(seqnames(cell_type))))



    plotTracks(list(
      iTrack, gTrack,
      gwasTrack,
      # ASoC_track,
      # alTrack_hMG, alTrack_hAst,
      alTrack_iMG, alTrack_iAst,
      alTrack_NGN2,
      alTrack_DN, alTrack_GA,
      snpTrack,
      grTrack
    ),
    sizes = c(
      0.5, 0.5,
      1,
      # 1,
      1, 1,
      1,
      1, 1,
      0.5,
      0.5
    ),
    chromosome = cell_type@ranges@NAMES,
    from = (cell_type@ranges@start - x_offset_1),
    to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
    transcriptAnnotation = "transcript",
    collapseTranscripts = "transcript")#,

  }

# BIN1 w/ SNP
plot_5_types(chr = "chr2",
             start = 127040000,
             end = 127145000,
             SNPposition = 127135484,
             SNPname = "rs6733839",
             GWASTrack = gRangesGWASTrack,
             ASoCTrack = gRangesASoCTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 200,
             gene_track_height = 0.5)

# PICALM w/ SNP
plot_5_types(chr = "chr11",
             start = 85894483,
             end = 86190000,
             SNPposition = 86156833,
             SNPname = "rs10792832",
             GWASTrack = gRangesGWASTrack,
             ASoCTrack = gRangesASoCTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 200,
             gene_track_height = 0.5)





# BIN1 w/ SNP
plot_5_types(chr = "chr2",
             start = 127050000,
             end = 127225000,
             SNPposition = 127135234,
             SNPname = "rs6733839",
             GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 300,
             gene_track_height = 0.5)

# BIN1 w/ SNP
plot_5_types(chr = "chr2",
             start = 127000000,
             end = 127202500,
             SNPposition = 127135234,
             SNPname = "rs6733839",
             # GWASTrack = gRangesGWASTrack,
             lineWidth = 0,
             minHeight = 0,
             ylimit = 400)
