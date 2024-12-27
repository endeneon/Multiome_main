# plot HiC track together with sites proximal to RERE gene
# input as BEDPE format
# need to lift hg19 to hg38 genome coordination
# process hippocampal data

# init
library(Gviz)
library(GenomicInteractions)
# library(GenomicRanges)
library(InteractionSet)
library(rtracklayer)
# library(ChIPseeker)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ensembldb)
library(org.Hs.eg.db)
# library(grDevices)
# library(gridExtra)
##########
options(ucscChromosomeNames = F)
##########

####
# hic_file <- file("RERE_hippocampal_chr1_ibed_hg38_final.bedpe")

# RERE_hippocampal_chr1_ibed_hg38_final$X8 <- RERE_hippocampal_chr1_ibed_hg38_final$X8 * 100
# write.table(RERE_hippocampal_chr1_ibed_hg38_final, file = "RERE_hippocampal_chr1_ibed_hg38_integer.bedpe",
#             quote = F, row.names = F, col.names = F, sep = "\t")
RERE_hippocampal_chr1_ibed_hg38_integer <- RERE_hippo_HiC_input_hg38_input
RERE_hippocampal_chr1_ibed_hg38_integer$X8 <- RERE_hippocampal_chr1_ibed_hg38_integer$X8 * 100

RERE_hippocampal_chr1_ibed_hg38_integer_section <- RERE_hippocampal_chr1_ibed_hg38_integer[(RERE_hippocampal_chr1_ibed_hg38_integer$X2 > 800000) |
                                                                                       (RERE_hippocampal_chr1_ibed_hg38_integer$X3 > 8000000) |
                                                                                       (RERE_hippocampal_chr1_ibed_hg38_integer$X5 > 8000000) |
                                                                                       (RERE_hippocampal_chr1_ibed_hg38_integer$X6 > 8000000), ]
RERE_hippocampal_chr1_ibed_hg38_integer_section <- RERE_hippocampal_chr1_ibed_hg38_integer_section[(RERE_hippocampal_chr1_ibed_hg38_integer_section$X2 < 9000000) |
                                                                                               (RERE_hippocampal_chr1_ibed_hg38_integer_section$X3 < 9000000) |
                                                                                               (RERE_hippocampal_chr1_ibed_hg38_integer_section$X5 < 9000000) |
                                                                                               (RERE_hippocampal_chr1_ibed_hg38_integer_section$X6 < 9000000), ]
RERE_hippocampal_chr1_ibed_hg38_integer_section <- RERE_hippocampal_chr1_ibed_hg38_integer[(RERE_hippocampal_chr1_ibed_hg38_integer$X2 > 8240000) &
                                                                                       (RERE_hippocampal_chr1_ibed_hg38_integer$X3 > 8240000) &
                                                                                       (RERE_hippocampal_chr1_ibed_hg38_integer$X5 > 8240000) &
                                                                                       (RERE_hippocampal_chr1_ibed_hg38_integer$X6 > 8240000), ]
RERE_hippocampal_chr1_ibed_hg38_integer_section <- RERE_hippocampal_chr1_ibed_hg38_integer_section[(RERE_hippocampal_chr1_ibed_hg38_integer_section$X2 < 9000000) &
                                                                                               (RERE_hippocampal_chr1_ibed_hg38_integer_section$X3 < 9000000) &
                                                                                               (RERE_hippocampal_chr1_ibed_hg38_integer_section$X5 < 9000000) &
                                                                                               (RERE_hippocampal_chr1_ibed_hg38_integer_section$X6 < 9000000), ]

RERE_hippocampal_chr1_ibed_hg38_integer_section_score <- RERE_hippocampal_chr1_ibed_hg38_integer_section[
  RERE_hippocampal_chr1_ibed_hg38_integer_section$X8 > 0,
  ]
write.table(RERE_hippocampal_chr1_ibed_hg38_integer_section_score, file = "RERE_hippocampal_chr1_ibed_hg38_integer_section_score.bedpe",
            quote = F, row.names = F, col.names = F, sep = "\t")

###
# RERE_hippocampal_chr1_ibed_hg38_final_section <- RERE_hippocampal_chr1_ibed_hg38_final[RERE_hippocampal_chr1_ibed_hg38_final$X8 > 4, ]
RERE_hippocampal_chr1_ibed_hg38_final_section_dedup <- RERE_hippocampal_chr1_ibed_hg38_final_section[!duplicated(RERE_hippocampal_chr1_ibed_hg38_final_section$X7), ]

RERE_hippocampal_chr1_ibed_hg38_final_section <- RERE_hippocampal_chr1_ibed_hg19
RERE_hippocampal_chr1_ibed_hg38_final_section <- RERE_hippocampal_chr1_ibed_hg38_final_section[
  RERE_hippocampal_chr1_ibed_hg38_final_section$X7 == "baitmap262", ]
# RERE_hippocampal_chr1_ibed_hg38_final_section <-
#   RERE_hippocampal_chr1_ibed_hg38_final_section[(RERE_hippocampal_chr1_ibed_hg38_final_section$X2 > 8000000) &
#                                                (RERE_hippocampal_chr1_ibed_hg38_final_section$X3 > 8000000) &
#                                                (RERE_hippocampal_chr1_ibed_hg38_final_section$X5 > 8000000) &
#                                                (RERE_hippocampal_chr1_ibed_hg38_final_section$X6 > 8000000), ]
# RERE_hippocampal_chr1_ibed_hg38_final_section <-
#   RERE_hippocampal_chr1_ibed_hg38_final_section[(RERE_hippocampal_chr1_ibed_hg38_final_section$X2 < 9000000) &
#                                                (RERE_hippocampal_chr1_ibed_hg38_final_section$X3 < 9000000) &
#                                                (RERE_hippocampal_chr1_ibed_hg38_final_section$X5 < 9000000) &
#                                                (RERE_hippocampal_chr1_ibed_hg38_final_section$X6 < 9000000), ]
RERE_hippocampal_chr1_ibed_hg38_final_section$X8 <- RERE_hippocampal_chr1_ibed_hg38_final_section$X8 * 100
RERE_hippocampal_chr1_ibed_hg38_final_section$X9 <- "*"
RERE_hippocampal_chr1_ibed_hg38_final_section$X10 <- "*"
RERE_hippocampal_chr1_ibed_hg38_final_section$X2 <- RERE_hippocampal_chr1_ibed_hg38_final_section$X2 - 60060
RERE_hippocampal_chr1_ibed_hg38_final_section$X3 <- RERE_hippocampal_chr1_ibed_hg38_final_section$X3 - 60060
RERE_hippocampal_chr1_ibed_hg38_final_section$X5 <- RERE_hippocampal_chr1_ibed_hg38_final_section$X5 - 60060
RERE_hippocampal_chr1_ibed_hg38_final_section$X6 <- RERE_hippocampal_chr1_ibed_hg38_final_section$X6 - 60060
RERE_hippocampal_chr1_ibed_hg38_final_section <- RERE_hippocampal_chr1_ibed_hg38_final_section[
  RERE_hippocampal_chr1_ibed_hg38_final_section$X8 > 500, ] # 400
# RERE_hippocampal_chr1_ibed_hg38_final_section <- RERE_hippocampal_chr1_ibed_hg38_final_section[-42, ]
write.table(RERE_hippocampal_chr1_ibed_hg38_final_section, file = "RERE_HiC_input_hg38_chr1_input_section_score.bedpe",
            quote = F, row.names = F, col.names = F, sep = "\t")
######

# RERE_HiC_input_hg38_chr1_input_section <-
#   RERE_HiC_input_hg38_chr1_input[(RERE_HiC_input_hg38_chr1_input$X2 > 8000000) &
#                                    (RERE_HiC_input_hg38_chr1_input$X3 > 8000000) &
#                                    (RERE_HiC_input_hg38_chr1_input$X5 > 8000000) &
#                                    (RERE_HiC_input_hg38_chr1_input$X6 > 8000000), ]
# RERE_HiC_input_hg38_chr1_input_section <-
#   RERE_HiC_input_hg38_chr1_input_section[(RERE_HiC_input_hg38_chr1_input_section$X2 < 9000000) &
#                                            (RERE_HiC_input_hg38_chr1_input_section$X3 < 9000000) &
#                                            (RERE_HiC_input_hg38_chr1_input_section$X5 < 9000000) &
#                                            (RERE_HiC_input_hg38_chr1_input_section$X6 < 9000000), ]
# RERE_HiC_input_hg38_chr1_input_section_score <- RERE_HiC_input_hg38_chr1_input_section[
#   RERE_HiC_input_hg38_chr1_input_section$X8 > 2,
#   ]
# write.table(RERE_HiC_input_hg38_chr1_input_section_score, file = "RERE_HiC_input_hg38_chr1_input_section_score.bedpe",
#             quote = F, row.names = F, col.names = F, sep = "\t")
#####
hic_data <- makeGenomicInteractionsFromFile("RERE_HiC_input_hg38_chr1_input_section_score.bedpe",
                                            type="bedpe",
                                            experiment_name = "RERE_HiC",
                                            description = "RERE_HiC_hg38")
colnames(RERE_BED6_eQTL) <- c("CHR", "START", "END", "ID", "SCORE", "STRAND")
RERE_BED6_eQTL$SCORE <- 0 - log10(RERE_BED6_eQTL$SCORE)

# gRanges_RERE_eQTL <- GRanges(seqnames = Rle("chr1", 1),
#                              seqinfo = Seqinfo(seqnames = RERE_BED6_eQTL$ID,
#                                                genome = "hg38"),
#                              ranges = IRanges(start = RERE_BED6_eQTL$START,
#                                               end = RERE_BED6_eQTL$END,
#                                               names = RERE_BED6_eQTL$CHR),
#                              mcols = as.data.frame(RERE_BED6_eQTL$SCORE),
#                              strand = Rle(strand(RERE_BED6_eQTL$STRAND)))
gRanges_RERE_eQTL <- GRanges(RERE_BED6_eQTL)
# write.table(RERE_BED6_eQTL,
#             file = "RERE_BED6_eQTL_log.txt",
#             quote = F,
#             row.names = F, col.names = F,
#             sep = "\t")
# gRanges_RERE_eQTL <-
#   lapply(split(RERE_BED6_eQTL, RERE_BED6_eQTL$ID), function(i){
#     GRanges(seqnames = i$ID,
#             ranges = IRanges(start = i$START,
#                              end = i$END,
#                              names = i$ID))
#   })
# gRanges_RERE_eQTL <- GRangesList(lapply(split(RERE_BED6_eQTL, RERE_BED6_eQTL$ID), function(i){
#   GRanges(seqnames = i$ID,
#           ranges = IRanges(start = i$START,
#                            end = i$END,
#                            names = i$ID))
# }))

# gRanges_RERE_eQTL <- readPeakFile("RERE_BED6_eQTL_log.txt")
###
# interaction_track <- InteractionTrack(hic_data, name = "HiC", chromosome = "chr1")
# plotTracks(interaction_track, chromosome="chr1",
#            from = 8340000, to = 8950000)
# 60060
######
hic_data <- makeGenomicInteractionsFromFile("RERE_HiC_input_hg38_chr1_input_section_score.bedpe",
                                            type="bedpe",
                                            experiment_name = "RERE_hippo_HiC",
                                            description = "RERE_HiC_hippo_hg38")
# plot_anywhere_interact(chr = "chr1", start = 8250000, end = 8900000,
#                        SNPname = "rs301791", SNPposition = 8408312,
#                        ylimit = 1200)


plot_anywhere_interact_seg(chr = "chr1", start = 8250000, end = 8900000,
                           SNPname = "rs301791", SNPposition = 8408312,
                           ylimit = 1200)

