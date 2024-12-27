# function for plotting
plot_rs1532278 <-
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
      DataTrack(range = GWASTrack,
                genome = "hg38",
                type = "p",
                # chromosome = "chr11",
                name = "Alz\nGWAS",
                background.title = "darkblue",
                fill.coverage = "darkblue",
                col.coverage = "darkblue",
                cex.axis = 0.8,
                cex.axis = 0.75)

    # "Alz_GWAS_Bellenguez_etal_Nat_Genet_2022")
    alTrack_hMG <-
      AlignmentsTrack("rs1532278_2MB_15Aug2023/hMG_downsampled_25M_het_rs1532278_2MB.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation = function(x) {x * 2},
                      ylim = c(0, ylimit),
                      background.title = "tomato4",
                      fill.coverage = "tomato4",
                      col.coverage = "tomato4",
                      cex.axis = 0.75,
                      cex.title = 0.8,
                      lwd = lineWidth,

                      name = "hMG")
    alTrack_hAst <-
      AlignmentsTrack("rs1532278_2MB_15Aug2023/hAst_downsampled_25M_het_rs1532278_2MB.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation=function(x) {x * 2},
                      ylim = c(0, ylimit),
                      background.title = "#1B9E77",
                      fill.coverage = "#1B9E77",
                      col.coverage = "#1B9E77",
                      cex.axis = 0.75,
                      cex.title = 0.8,
                      lwd = lineWidth,
                      lwd.border = lineWidth,
                      min.width = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "hAst")

    alTrack_iMG <-
      AlignmentsTrack("rs1532278_2MB_15Aug2023/MG_downsampled_50M_het_rs1532278_2MB.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation = function(x) {x * 2},
                      ylim = c(0, ylimit),
                      background.title = "brown",
                      fill.coverage = "brown",
                      col.coverage = "brown",
                      cex.axis = 0.75,
                      cex.title = 0.8,
                      lwd = lineWidth,

                      name = "iMG")
    alTrack_iAst <-
      AlignmentsTrack("rs1532278_2MB_15Aug2023/Ast_downsampled_50M_het_rs1532278_2MB.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation=function(x) {x * 2},
                      ylim = c(0, ylimit),
                      background.title = "deepskyblue4",
                      fill.coverage = "deepskyblue4",
                      col.coverage = "deepskyblue4",
                      cex.axis = 0.75,
                      cex.title = 0.8,
                      lwd = lineWidth,
                      lwd.border = lineWidth,
                      min.width = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "iAst")
    alTrack_NGN2 <-
      AlignmentsTrack("rs1532278_2MB_15Aug2023/NGN2_downsampled_50M_het_rs1532278_2MB.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation = function(x) {x * 2},
                      background.title = "#D95F02",
                      fill.coverage = "#D95F02",
                      col.coverage = "#D95F02",
                      ylim = c(0, ylimit),
                      cex.axis = 0.75,
                      cex.title = 0.8,
                      lwd = lineWidth,
                      lwd.border = lineWidth,
                      min.width = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "iGlut")

    alTrack_DN <-
      AlignmentsTrack("rs1532278_2MB_15Aug2023/DN_downsampled_50M_het_rs1532278_2MB.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation = function(x) {x * 2},
                      ylim = c(0, ylimit),
                      background.title = "#7570B3",
                      fill.coverage = "#7570B3",
                      col.coverage = "#7570B3",
                      cex.axis = 0.75,
                      cex.title = 0.8,
                      lwd = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "iDN")
    alTrack_GA <-
      AlignmentsTrack("rs1532278_2MB_15Aug2023/GA_downsampled_50M_het_rs1532278_2MB.bam",
                      isPaired = T, coverageOnly = T,
                      chromosome = as.character(unique(seqnames(cell_type))),
                      genome = "hg38", type = "coverage",
                      transformation = function(x) {x * 2},
                      ylim = c(0, ylimit),
                      background.title = "#E7298A",
                      fill.coverage = "#E7298A",
                      col.coverage = "#E7298A",
                      cex.title = 0.8,
                      lwd = lineWidth,
                      lwd.border = lineWidth,
                      min.width = lineWidth,
                      min.height = minHeight,
                      minCoverageHeight = minHeight,
                      name = "iGA")

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
    grTrack <- ucscGenes
    ranges(grTrack) <- z
    grTrack@dp@pars$col.line <- "black"
    grTrack@dp@pars$fontcolor <- "black"
    # grTrack@name <- paste("RefSeq", "Gene", collapse = "\n")
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
      1, 1,
      # 1, 1,
      1,
      1,1,
      0.5,
      0.5
    ),
    chromosome = cell_type@ranges@NAMES,
    from = (cell_type@ranges@start - x_offset_1),
    to = (cell_type@ranges@start + cell_type@ranges@width + x_offset_2),
    transcriptAnnotation = "transcript",
    collapseTranscripts = "transcript")#,

  }
