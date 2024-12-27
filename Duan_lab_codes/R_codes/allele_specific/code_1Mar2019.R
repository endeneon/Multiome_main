# plot ASoC peaks
# need to overlay 3 data tracks
# Siwei 1 March 2019

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

# rs2027349 coordination: chr1:150067621

plot_AsoC_fine_dtails(chr = "chr1", start = 150067560, end = 150067680, ylim = 450,
                cell_type_to_plot = "CN", title_name = " ") # details


  plot_AsoC_fine_dtails(chr = "chr1", start = 150067560, end = 150067680, ylim = 450,
                        cell_type_to_plot = "CN", title_name = " ") # details
# NPC_fine_detail <-
  plot_AsoC_fine_dtails(chr = "chr1", start = 150067560, end = 150067680, ylim = 450,
                        cell_type_to_plot = "NPC", title_name = " ") # details


plot_AsoC_composite(chr = "chr1", start = 149940560, end = 150200000, ylim = 500,
                        cell_type_to_plot = "NPC", title_name = " ") # details


plot_AsoC_peaks(chr = "chr1", start = 150067200, end = 150067950, ylim = 1000,
                cell_type_to_plot = "NPC", title_name = " ")


plot_AsoC_peaks(chr = "chr1", start = 150007200, end = 150100000, ylimit = 300,
                cell_type_to_plot = "CN", title_name = " ")










rs2027349_coordinate <- as.data.frame(marker_gene_list[1, ])
rs2027349_coordinate[1, ] <- c("chr1", as.numeric(145000001), as.numeric(155000000-1), "rs2027349", 100, "*")
rs2027349_coordinate$X2 <- as.numeric(rs2027349_coordinate$X2)
rs2027349_coordinate$X3 <- as.numeric(rs2027349_coordinate$X3)
rs2027349_coordinate$X5 <- as.numeric(rs2027349_coordinate$X5)


rs2027349_GRange <-
  GRanges(seqnames = Rle(rs2027349_coordinate$X1[1]),
          seqinfo = Seqinfo(seqnames = unique(rs2027349_coordinate$X1[1]),
                            # isCircular = F,
                            genome = "hg38"),
          # chromosomes <- Rle(marker_gene_list$X1),
          ranges = IRanges(start = rs2027349_coordinate$X2[1],
                           end = rs2027349_coordinate$X3[1],
                           names = rs2027349_coordinate$X1[1]),
          mcols = as.data.frame(rs2027349_coordinate$X5[1]),
          strand = Rle(strand(rs2027349_coordinate$X6[1])))


