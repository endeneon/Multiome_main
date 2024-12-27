# Siwei 11 Mar 2019

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


plot_AsoC_fine_dtails(chr = "chr1", start = 149887889, end = 150476566, ylim = 5000,
                      cell_type_to_plot = "CN", title_name = " ", transformation_factor = 4.078446) # details

plot_AsoC_fine_dtails(chr = "chr1", start = 150067560, end = 150067680, ylim = 400,
                      cell_type_to_plot = "NSC", title_name = " ", transformation_factor = 4.21563) # details

plot_AsoC_fine_dtails(chr = "chr1", start = 150067560, end = 150067680, ylim = 1200,
                      cell_type_to_plot = "DN", title_name = " ", transformation_factor = 9.182318) # details

plot_AsoC_fine_dtails(chr = "chr1", start = 150067560, end = 150067680, ylim = 1000,
                      cell_type_to_plot = "GA", title_name = " ", transformation_factor = 5.278966) # details

plot_AsoC_fine_dtails(chr = "chr1", start = 150067560, end = 150067680, ylim = 300,
                      cell_type_to_plot = "ips", title_name = " ", transformation_factor = 4.560793) # details


plot_AsoC_composite(chr = "chr1", start = 149940560, end = 150200000, ylim = 500,
                    cell_type_to_plot = "NPC", title_name = " ") # details


plot_AsoC_peaks(chr = "chr1", start = 150067200, end = 150067950, ylim = 1000,
                cell_type_to_plot = "NPC", title_name = " ")

plot_anywhere(chr = "chr1", start = 149887589, end = 150476866, ylimit = 8000)


plot_anywhere(chr = "chr12", start = 7786000, end = 7800000, ylimit = 3000)

