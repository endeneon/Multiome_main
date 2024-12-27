# Feb 5 2019 Siwei
# plot the pileup figures of
# cell-type-specific genes
# ATAC-Seq

# init
library(Gviz)
# data("cpgIslands")
library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ensembldb)
library(org.Hs.eg.db)
##########
marker_gene_list <- data.frame()
i <- 1
options(ucscChromosomeNames = F)


marker_gene_list$X1[5] <- "chr3"

marker_gene_GRange <-
    GRanges(seqnames = Rle(marker_gene_list$X1),
            seqinfo = Seqinfo(seqnames = unique(marker_gene_list$X1),
                              # isCircular = F,
                              genome = "hg38"),
            # chromosomes <- Rle(marker_gene_list$X1),
            ranges = IRanges(start = marker_gene_list$X2,
                             end = marker_gene_list$X3,
                             names = marker_gene_list$X1),
                             mcols = as.data.frame(marker_gene_list$X5),
                             strand = Rle(strand(marker_gene_list$X6)))

# marker_gene_Grange_list <-
#   lapply(split(marker_gene_list, marker_gene_list$X4), function(i) {
#     # colnames(i)
#     GRanges(seqnames = Rle(i$X1),
#             ranges = IRanges(start = i$X2,
#                              end = i$X3,
#                              names = i$X4),
#             seqinfo = Seqinfo(genome = "hg38"))
#   })
# gencode_v28_TxDb <- makeTxDbFromGFF("~/1TB/Databases/hg38/gencode.v28.annotation.gtf",
#                                     organism = "Homo sapiens")
TxDb.Hsapiens.UCSC.refGene.hg38 <- makeTxDbFromUCSC(genome = "hg38",
                                              tablename = "refGene")




marker_gene_genome <- genome(marker_gene_GRange)
geneSymbol_track <- AnnotationTrack(marker_gene_GRange, name = "Marker Genes")

annoTrack <- BiomartGeneRegionTrack(genome="GRCh38.p12", name="ENSEMBL", symbol="NANOG")

chr <- as.character(unique(seqnames(marker_gene_GRange)))
gTrack <- GenomeAxisTrack()
iTrack <- IdeogramTrack(genome = "hg38", chromosome = chr)
plotTracks(list(geneSymbol_track))
plotTracks(list(geneSymbol_track, gTrack))
plotTracks(list(geneSymbol_track, gTrack, iTrack))

############

GA_GRange <-
  GRanges(seqnames = Rle(marker_gene_list$X1[1]),
          seqinfo = Seqinfo(seqnames = unique(marker_gene_list$X1[1]),
                            # isCircular = F,
                            genome = "hg38"),
          # chromosomes <- Rle(marker_gene_list$X1),
          ranges = IRanges(start = marker_gene_list$X2[1],
                           end = marker_gene_list$X3[1],
                           names = marker_gene_list$X1[1]),
          mcols = as.data.frame(marker_gene_list$X5[1]),
          strand = Rle(strand(marker_gene_list$X6[1])))

DN_GRange <-
  GRanges(seqnames = Rle(marker_gene_list$X1[2]),
          seqinfo = Seqinfo(seqnames = unique(marker_gene_list$X1[2]),
                            # isCircular = F,
                            genome = "hg38"),
          # chromosomes <- Rle(marker_gene_list$X1),
          ranges = IRanges(start = marker_gene_list$X2[2],
                           end = marker_gene_list$X3[2],
                           names = marker_gene_list$X1[2]),
          mcols = as.data.frame(marker_gene_list$X5[2]),
          strand = Rle(strand(marker_gene_list$X6[2])))

CN_GRange <-
  GRanges(seqnames = Rle(marker_gene_list$X1[3]),
          seqinfo = Seqinfo(seqnames = unique(marker_gene_list$X1[3]),
                            # isCircular = F,
                            genome = "hg38"),
          # chromosomes <- Rle(marker_gene_list$X1),
          ranges = IRanges(start = marker_gene_list$X2[3],
                           end = marker_gene_list$X3[3],
                           names = marker_gene_list$X1[3]),
          mcols = as.data.frame(marker_gene_list$X5[3]),
          strand = Rle(strand(marker_gene_list$X6[3])))

iPS_GRange <-
  GRanges(seqnames = Rle(marker_gene_list$X1[4]),
          seqinfo = Seqinfo(seqnames = unique(marker_gene_list$X1[4]),
                            # isCircular = F,
                            genome = "hg38"),
          # chromosomes <- Rle(marker_gene_list$X1),
          ranges = IRanges(start = marker_gene_list$X2[4],
                           end = marker_gene_list$X3[4],
                           names = marker_gene_list$X1[4]),
          mcols = as.data.frame(marker_gene_list$X5[4]),
          strand = Rle(strand(marker_gene_list$X6[4])))

NPC_GRange <-
  GRanges(seqnames = Rle(marker_gene_list$X1[5]),
          seqinfo = Seqinfo(seqnames = unique(marker_gene_list$X1[5]),
                            # isCircular = F,
                            genome = "hg38"),
          # chromosomes <- Rle(marker_gene_list$X1),
          ranges = IRanges(start = marker_gene_list$X2[5],
                           end = marker_gene_list$X3[5],
                           names = marker_gene_list$X1[5]),
          mcols = as.data.frame(marker_gene_list$X5[5]),
          strand = Rle(strand(marker_gene_list$X6[5])))
#######
aTrack <- AnnotationTrack(range = DN_GRange, genome = genome(DN_GRange),
                                    name = "Marker Genes",
                                    chromosome = as.character(unique(seqnames(DN_GRange)))
                                    )


iTrack <- IdeogramTrack(genome = genome(DN_GRange),
                        chromosome = as.character(unique(seqnames(DN_GRange))))
# head(displayPars(grTrack))

###


# plotTracks(list(ucscGenes2), chromosome = "chr2",
#            from = 1e+06,
#            to = 2e+06,
#            transcriptAnnotation = "symbol",
#            collapseTranscripts = "transcript")

#####
# bamFile <- system.file("~/3TB/05Feb2019_Gviz_plot/peak_plotting_4_Gviz/DN_all_sorted_gviz_plot.bam",
#                        package = "Gviz")
#
# dTrack <- DataTrack(range = ,
#                     genome = "hg38", type = "h",
#                     name = "Coverate",
#                     window = -1,
#                     chromosome = as.character(unique(seqnames(DN_GRange))))


alTrack_CN <- AlignmentsTrack("~/3TB/05Feb2019_Gviz_plot/peak_plotting_4_Gviz/CN_all_sorted_gviz_plot.bam",
                           isPaired = T, coverageOnly = T,
                           chromosome = as.character(unique(seqnames(DN_GRange))),
                           genome = "hg38", type = "coverage",
                           transformation=function(x) {x * 4.078446},
                           ylim = c(0, 700),
                           name = "Glut")
alTrack_NPC <- AlignmentsTrack("~/3TB/05Feb2019_Gviz_plot/peak_plotting_4_Gviz/NSC_all_sorted_gviz_plot.bam",
                               isPaired = T, coverageOnly = T,
                               chromosome = as.character(unique(seqnames(DN_GRange))),
                               genome = "hg38", type = "coverage",
                               transformation=function(x) {x * 4.21563},
                               ylim = c(0, 700),
                               name = "NPC")
alTrack_DN <- AlignmentsTrack("~/3TB/05Feb2019_Gviz_plot/peak_plotting_4_Gviz/DN_all_sorted_gviz_plot.bam",
                              isPaired = T, coverageOnly = T,
                              chromosome = as.character(unique(seqnames(DN_GRange))),
                              genome = "hg38", type = "coverage",
                              transformation=function(x) {x * 9.182318},
                              ylim = c(0, 700),
                              name = "Dopa")
alTrack_GA <- AlignmentsTrack("~/3TB/05Feb2019_Gviz_plot/peak_plotting_4_Gviz/GA_all_sorted_gviz_plot.bam",
                              isPaired = T, coverageOnly = T,
                              chromosome = as.character(unique(seqnames(DN_GRange))),
                              genome = "hg38", type = "coverage",
                              transformation=function(x) {x * 5.278966},
                              ylim = c(0, 700),
                              name = "GABA")
alTrack_iPS <- AlignmentsTrack("~/3TB/05Feb2019_Gviz_plot/peak_plotting_4_Gviz/ips_all_sorted_gviz_plot.bam",
                              isPaired = T, coverageOnly = T,
                              chromosome = as.character(unique(seqnames(DN_GRange))),
                              genome = "hg38", type = "coverage",
                              transformation=function(x) {x * 4.560793},
                              ylim = c(0, 700),
                              name = "iPS")


########### DN
ucscGenes <- UcscTrack(genome=genome(DN_GRange), table="ncbiRefSeq",
                       track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome = as.character(unique(seqnames(DN_GRange))),
                       rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)
z <- ranges(ucscGenes)
mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "",
                                                           mcols(z)$symbol), "SYMBOL","REFSEQ"))
grTrack <- ucscGenes
ranges(grTrack) <- z

plotTracks(list(iTrack, gTrack, alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS, grTrack),
           sizes = c(1,1,1,1,1,1,1,1),
           chromosome = DN_GRange@ranges@NAMES,
           from = (DN_GRange@ranges@start + 7000),
           to = (DN_GRange@ranges@start + DN_GRange@ranges@width + 8000),
           transcriptAnnotation = "transcript",
           collapseTranscripts = "transcript")#,

########### GA
ucscGenes <- UcscTrack(genome=genome(GA_GRange), table="ncbiRefSeq",
                       track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome = as.character(unique(seqnames(GA_GRange))),
                       rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "",
                                                           mcols(z)$symbol), "SYMBOL","REFSEQ"))
grTrack <- ucscGenes
ranges(grTrack) <- z

plotTracks(list(iTrack, gTrack, alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS, grTrack),
           sizes = c(1,1,1,1,1,1,1,1),
           chromosome = GA_GRange@ranges@NAMES,
           from = (GA_GRange@ranges@start - 10000),
           to = (GA_GRange@ranges@start + GA_GRange@ranges@width + 100000),
           transcriptAnnotation = "transcript",
           collapseTranscripts = "transcript")

########### iPS
ucscGenes <- UcscTrack(genome=genome(iPS_GRange), table="ncbiRefSeq",
                       track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome = as.character(unique(seqnames(iPS_GRange))),
                       rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$transcript <- as.vector(mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "",
                                                           mcols(z)$symbol), "SYMBOL","REFSEQ"))
grTrack <- ucscGenes
ranges(grTrack) <- z

plotTracks(list(iTrack, gTrack, alTrack_CN, alTrack_NPC, alTrack_DN, alTrack_GA, alTrack_iPS, grTrack),
           sizes = c(1,1,1,1,1,1,1,1),
           chromosome = iPS_GRange@ranges@NAMES,
           from = (iPS_GRange@ranges@start - 10000),
           to = (iPS_GRange@ranges@start + iPS_GRange@ranges@width + 10000),
           transcriptAnnotation = "transcript",
           collapseTranscripts = "transcript")
##########
alTrack_iPS <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/4TB_1/core_8_peaks_FDR001/merged_BAM/ips_all_sorted.bam",
                               isPaired = T, coverageOnly = T,
                               chromosome = as.character(unique(seqnames(DN_GRange))),
                               genome = "hg38", type = "coverage",
                               transformation=function(x) {x * 4.560793},
                               ylim = c(0, 700),
                               name = "iPS")

plotTracks(list(alTrack_iPS),
           # sizes = c(1,1,1,1,1,1,1,1),
           chromosome = iPS_GRange@ranges@NAMES,
           from = (iPS_GRange@ranges@start - 10000),
           to = (iPS_GRange@ranges@start + iPS_GRange@ranges@width + 10000),
           transcriptAnnotation = "transcript",
           collapseTranscripts = "transcript")
##########

library(scales)
show_col(hue_pal()(6))



#################
1/(840923614/3429661562) #CN 0.2451914, 4.078446
1/(373507158/3429661562) #DN 0.108905, 9.182318
1/(649684371/3429661562) #GA 0.189431, 5.278966
1/(751987961/3429661562) #iPS 0.2192601, 4.560793
1/(813558458/3429661562)  #NPC 0.2372125, 4.21563
