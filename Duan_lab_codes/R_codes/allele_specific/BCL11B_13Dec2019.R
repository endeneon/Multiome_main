# plot BCL11B-related peaks
# Siwei 13 Dec 2019


# init
library(Gviz)
# data("cpgIslands")
library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ensembldb)
library(org.Hs.eg.db)
library(grDevices)
library(gridExtra)
##########
options(ucscChromosomeNames = F)
##########

plot_anywhere(chr = "chr14", start = 99159000, end = 99425000,
              SNPname = "rs12895055", SNPposition = 99246457,
              ylimit = 6000)

plot_AsoC_composite_BCL11B(chr = "chr14", start = 99246200, end = 99246600,
                           ylimit = 6500)
plot_VPS45_isogenic(chr = "chr1",
              start = 150066791,
              end = 150068025, ylimit = 600)

plot_VPS45_isogenic_45M(chr = "chr1",
                    start = 150066791,
                    end = 150068025, ylimit = 600)

plot_VPS45_isogenic_45M(chr = "chr1",
                        start = 150066891,
                        end = 150067600, ylimit = 600)

plot_VPS45_isogenic_45M(chr = "chr1", # designated region
                        start = 150066901,
                        end = 150067900, ylimit = 600)

plot_VPS45_6_lines_q20(chr = "chr1",
                    start = 150066701,
                    end = 150068025, ylimit = 300)

plot_VPS45_6_lines_q20_45M(chr = "chr1",
                       start = 150066901,
                       end = 150067900, ylimit = 300)


plot_VPS45_isogenic_45M(chr = "chr1",
                        start = 150067611,
                        end = 150067631, ylimit = 600)
##############

plot_AsoC_BCL11B_2sites(chr = "chr14",
                        start = 99246300,
                        end = 99246750,
                        ylimit = 1500)
