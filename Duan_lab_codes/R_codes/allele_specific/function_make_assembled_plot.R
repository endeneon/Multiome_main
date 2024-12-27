# function for plotting
make_assembled_plot <- function(cell_type, x_offset_1 = 10000, x_offset_2 = 10000, ylimit = 800) {
  # aTrack <- AnnotationTrack(range = cell_type, genome = genome(cell_type),
  #                           name = "Marker Genes",
  #                           chromosome = as.character(unique(seqnames(cell_type))))

  print(as.character(unique(seqnames(cell_type))))

  iTrack <- IdeogramTrack(genome = genome(cell_type),
                          chromosome = as.character(unique(seqnames(cell_type))))

  alTrack_CN <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/4TB_1/core_8_peaks_FDR001/merged_BAM/CN_all_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                transformation=function(x) {x * 4.078446},
                                ylim = c(0, ylimit),
                                name = "Glut")
  alTrack_NPC <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/4TB_1/core_8_peaks_FDR001/merged_BAM/NSC_all_sorted.bam",
                                 isPaired = T, coverageOnly = T,
                                 chromosome = as.character(unique(seqnames(cell_type))),
                                 genome = "hg38", type = "coverage",
                                 transformation=function(x) {x * 4.21563},
                                 ylim = c(0, ylimit),
                                 name = "NPC")
  alTrack_DN <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/4TB_1/core_8_peaks_FDR001/merged_BAM/DN_all_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                transformation=function(x) {x * 9.182318},
                                ylim = c(0, ylimit),
                                name = "Dopa")
  alTrack_GA <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/4TB_1/core_8_peaks_FDR001/merged_BAM/GA_all_sorted.bam",
                                isPaired = T, coverageOnly = T,
                                chromosome = as.character(unique(seqnames(cell_type))),
                                genome = "hg38", type = "coverage",
                                transformation=function(x) {x * 5.278966},
                                ylim = c(0, ylimit),
                                name = "GABA")
  alTrack_iPS <- AlignmentsTrack("/run/user/1000/gvfs/sftp:host=192.168.2.5/home/cpg/4TB_1/core_8_peaks_FDR001/merged_BAM/ips_all_sorted.bam",
                                 isPaired = T, coverageOnly = T,
                                 chromosome = as.character(unique(seqnames(cell_type))),
                                 genome = "hg38", type = "coverage",
                                 transformation=function(x) {x * 4.560793},
                                 ylim = c(0, ylimit),
                                 name = "iPS")

  ########### plotting
  ucscGenes <- UcscTrack(genome=genome(cell_type), table="ncbiRefSeq",
                         track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                         chromosome = as.character(unique(seqnames(cell_type))),
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
             chromosome = cell_type@ranges@NAMES,
             from = (cell_type@ranges@start - x_offset_1),
             to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
             transcriptAnnotation = "transcript",
             collapseTranscripts = "transcript")#,

}
