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
}
#
options(ucscChromosomeNames = F)
options(Gviz.ucscUrl = "https://genome-euro.ucsc.edu/cgi-bin/")

# load data
df_GWASTrack <-
  read_delim("COX7_GWAS_Alz_2022_Jansen_et_al_rs62375391_hg38.bedGraph",
             delim = "\t", escape_double = FALSE,
             col_names = FALSE, trim_ws = TRUE)
colnames(df_GWASTrack) <-
  c("chr", "start", "end", "score")
df_GWASTrack$score <-
  0 - log10(df_GWASTrack$score)

gRangesGWASTrack <-
  makeGRangesFromDataFrame(df = df_GWASTrack,
                           keep.extra.columns = T,
                           ignore.strand = T)

# test DisplayPars #####
alTrack_test <-
  AlignmentsTrack("rs1532278_2MB_15Aug2023/MG_downsampled_50M_het_rs1532278_2MB.bam",
                  chromosome = "chr8",
                  start = 27600000,
                  end = 27620000,
                  isPaired = T, coverageOnly = T,
                  # chromosome = as.character(unique(seqnames(cell_type))),
                  genome = "hg38", type = "coverage",
                  transformation = function(x) {x * 1},
                  ylim = c(0, 100),
                  background.title = "brown",
                  fill.coverage = "brown",
                  col.coverage = "brown",
                  cex.axis = 1.2,
                  name = "MG")

alTrack_dp <-
  displayPars(alTrack_test)
alTrack_dp$minCoverageHeight


# make the plot #####
plot_rs1532278(chr = "chr8",
               start = 27600000,
               end = 27635000,
               SNPposition = 27608798,
               SNPname = "rs1532278",
               GWASTrack = gRangesGWASTrack,
               lineWidth = 0,
               minHeight = 0,
               ylimit = 150)

plot_rs62375391(chr = "chr5",
               start = 86600000,
               end = 87400000,
               SNPposition = 86997803,
               SNPname = "rs62375391",
               GWASTrack = gRangesGWASTrack,
               lineWidth = 0,
               minHeight = 0,
               ylimit = 150)


plot_rs10792832(chr = "chr11",
                start = 86015483,
                end = 86175000,
                SNPposition = 86156833,
                SNPname = "rs10792832",
                GWASTrack = gRangesGWASTrack,
                ylimit = 300)

plot_rs10792832(chr = "chr11",
                start = 86060124,
                end = 86175000,
                SNPposition = 86156833,
                SNPname = "rs10792832",
                GWASTrack = gRangesGWASTrack,
                ylimit = 300)

plot_rs10792832(chr = "chr11",
                start = 85985483,
                end = 86200000,
                SNPposition = 86156833,
                SNPname = "rs10792832",
                GWASTrack = gRangesGWASTrack,
                ylimit = 300)

plot_rs10792832(chr = "chr11",
                start = 85885483,
                end = 86200000,
                SNPposition = 86156833,
                SNPname = "rs10792832",
                GWASTrack = gRangesGWASTrack,
                ylimit = 300)

plot_rs10792832(chr = "chr11",
                start = 86156500,
                end = 86157000,
                SNPposition = 86156833,
                SNPname = "rs10792832",
                GWASTrack = gRangesGWASTrack,
                ylimit = 300)

plot_rs10792832(chr = "chr11",
                start = 86150500,
                end = 86169000,
                SNPposition = 86156833,
                SNPname = "rs10792832",
                GWASTrack = gRangesGWASTrack,
                ylimit = 300)

plot_AsoC_peaks(chr = "chr11",
                start = 86156600,
                end = 86157000,
                ylimit = 80,
                title_name = " ")

brewer.pal(8, name = "Dark2")
