# plot practise
# Siwei 6 Feb 2019

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
##########

options(ucscChromosomeNames = F)

make_assembled_plot(cell_type = DN_GRange) # pass, need fine adjustment on location
make_assembled_plot(cell_type = CN_GRange, x_offset_1 = 0, x_offset_2 = 100000) # need to choose another
make_assembled_plot(cell_type = GA_GRange) # need to reposition
make_assembled_plot(cell_type = NPC_GRange) # need to reposition
make_assembled_plot(cell_type = iPS_GRange, x_offset_1 = -6000, x_offset_2 = 2000, ylimit = 3000) # pass

plot_anywhere(chr = "chr10", start = 26200000, end = 26280000, ylimit = 3000) # GAD2



plot_anywhere(chr = "chr12", start = 7785000, end = 7800000, ylimit = 3000) # NANOG
plot_anywhere(chr = "chr15", start = 76336044, end = 76339000, ylimit = 3000) # ISL2
plot_anywhere(chr = "chr10", start = 87088997, end = 87097220, ylimit = 6000) # GLUD1

##########






