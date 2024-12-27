### An example of making annotation with exon regions.
ex <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, "gene")
exon <- exons(TxDb.Hsapiens.UCSC.hg38.knownGene)
exon1 <- exon[!exon %over% unlist(ex),]
exon1 <- exon1[exon1@seqnames %in% paste0("chr",1:22),]
exon.query <- countOverlaps(snprange,exon1)
eq1 <- exon.query
eq1[eq1>0] <- 1
write.table(data.frame(snp=snprange$snp,exon_d=eq1),"eQTL_enrichment/annotation/exon.txt" ,row.names = F, quote=F,sep="\t")
system("gzip eQTL_enrichment/annotation/exon.txt")

### Execute command line with TORUS in R and extract the enrichment results.
torus.eqtl <- function(annot_file, matrixeqtl_res, smap, gmap,
                       torus_path="/path/to/torus") {
  # annot_file: location of the annotation file
  # matrixeqtl_res: location of the eQTL summary statistics
  # smap: snp position file
  # gmap: gene position file
  torus_args <- c("-d", matrixeqtl_res, "-annot", annot_file, 
                  #"-smap", smap, "-gmap", gmap,
                  "-est")
  res <- processx::run(command = torus_path, args = torus_args, 
                       echo_cmd = TRUE, echo = TRUE)
  enrich <- read.table(file = textConnection(res$stdout), 
                       skip = 1, header = FALSE, stringsAsFactors = FALSE)
  colnames(enrich) <- c("term", "estimate", "low", "high")
  return(enrich)
}


