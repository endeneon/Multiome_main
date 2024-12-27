# Siwei 05 Jul 2023
# plot 1MB proximal region of rs10792832
# chr11:86156833

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
#
options(ucscChromosomeNames = F)

# load data
df_GWASTrack <-
  read_delim("GWAS_Alz_2022_Jansen_et_al_rs10792832_hg38.bedGraph",
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

# make the plot #####
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
