library(ggplot2)
library(grid)
library(AnnotationHub)
library(ggrepel)
library(patchwork)
library(gtable)
library(data.table)

results_dir <- "/project/xinhe/lifan/neuron_stim/ctwas_SCZ/results/cellgroup_SCZ/"
analysis_id <- c("GABA","npglut","nmglut")
setwd("/project/xinhe/lifan/neuron_stim/ctwas_SCZ")

ctwas.res <- list()
ctwas.res[[analysis_id[1]]] <- read.table(paste0(results_dir,analysis_id[1],"/", analysis_id[1],".susieIrss.txt"), header=T)
ctwas.res[[analysis_id[2]]] <- read.table(paste0(results_dir,analysis_id[2],"/", analysis_id[2],".susieIrss.txt"), header=T)
ctwas.res[[analysis_id[3]]] <- read.table(paste0(results_dir,analysis_id[3],"/", analysis_id[3],".susieIrss.txt"), header=T)

snps <- fread("../PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv")
ctwas.res <- read.table("results/Joint_caQTL_ASoC_SCZ_09112024/joint_celltype_cover9.susieIrss.txt")

load.zgene <- function(celltype) {
  ld_exprfs <- paste0("results/cellgroup_caQTL_SCZ/",celltype,"/",celltype, "_chr", 1:22, ".expr.gz")
  z_gene <- list()
  for (i in 1:22){
    load(paste0("results/cellgroup_caQTL_SCZ/",celltype,"/",celltype, "_chr", i, ".exprqc.Rd"))
    z_gene[[i]] <- z_gene_chr
  }
  z_gene <- do.call(rbind, z_gene)
  #rm(qclist, wgtlist, z_gene_chr)
  #rownames(z_gene) <- NULL
  z_gene
}
z_gene <- list()
z_gene[["npglut"]] <- load.zgene("npglut")
z_gene[["nmglut"]] <- load.zgene("nmglut")
z_gene[["GABA"]] <- load.zgene("GABA")

construct.dat <- function(res.sub) {
  dat <- res.sub
  dat$name <- sapply(strsplit(dat$id,"\\|"), function(x) {x[1]})
  dat$time_point <- sapply(strsplit(dat$type,"_"), function(x) {x[1]})
  # dat$celltype <- factor(dat$celltype, levels=c("GABA","npglut","nmglut"))
  dat$pos <- dat$pos / 1e6
  genes <- unique(dat[dat$type!="SNP","name"])
  loc$TX <- loc$TX[loc$TX$gene_name %in% genes,]
  dat$label <- ""
  dat$label[dat$name %in% focus] <- dat$id[dat$name %in% focus]
  dat
}

focus <- "chr8_143933798_143934298"
focus <- "chr17_19317206_19317706"
focus <- "chr20_46051914_46052414"

focus <- "chr7_2010548_2011048"

locus_plot_caqtl <- function(snp, assoc, focus, celltype) {
  # snp: the SNP to be investigated
  # assoc: the associated peak ID in certain time point
  # focus: the peaks
  # celltype: the cell type for plotting
  
  chr.str <- strsplit(focus,"_")[[1]][1]
  chr <- as.numeric(gsub("chr","",chr.str))
  xrange <- as.numeric(strsplit(focus,"_")[[1]][2:3])
  xrange[1] <- xrange[1] - 2.5e5
  xrange[2] <- xrange[2] + 2.5e5
  res.sub <- ctwas.res[[celltype]][(ctwas.res[[celltype]]$region_tag1==chr) &
                                   (ctwas.res[[celltype]]$pos>=xrange[1]) &
                                   (ctwas.res[[celltype]]$pos<=xrange[2]),]
  
  dat.ct <- construct.dat(res.sub)
  dat.ct$type[dat.ct$id %in% assoc] <- "linked"
  
  snp.sub <- snps[snps$CHROM==chr,c("ID","PVAL")]
  z_gene2 <- merge(z_gene[[celltype]], dat.ct, by="id", all.x=F, all.y=F)
  snp.sub1 <- merge(snp.sub, dat.ct[,1:3], by.x="ID", by.y="id", all.x=F, all.y=F)
  snp.sub1$logP <- -log10(snp.sub1$PVAL)
  snp.sub1$type <- "SNP"
  #merge(z_snp[,c(1,4)], datNP, by="id", all.x=F, all.y=F)
  z_gene2$PVAL <- pnorm(abs(z_gene2$z), lower.tail = F)
  z_gene2$logP <- -log10(z_gene2$PVAL)
  colnames(z_gene2)[1] <- "ID"
  datP <- rbind(snp.sub1, z_gene2[,c("ID","PVAL","chrom","pos","logP","type")])
  datP$label <- ""
  datP$label[startsWith(datP$ID,focus)] <- datP$ID[startsWith(datP$ID,focus)]
  datP$type[datP$ID %in% assoc] <- "linked"
  
  
  alpha_vec <- rep(2.0,length(unique(datP$type)))
  names(alpha_vec) <- unique(datP$type)
  alpha_vec["SNP"] <- 1
  
  size_vec <- rep(2,length(unique(datP$type)))
  names(size_vec) <- unique(datP$type)
  size_vec["SNP"] <-0.5
  
  shape_vec <- rep(17,length(unique(datP$type)))
  names(shape_vec) <- unique(datP$type)
  shape_vec["SNP"] <-20
  shape_vec["linked"] <- 11
  
  snp.line <- dat.ct[dat.ct$id==snp,]
  
  p1 <- ggplot(datP, aes(x=pos, y=logP, shape=type, color=type, size=type, alpha=type, label=label)) + 
    geom_point() +
    # facet_grid(celltype~., scales = "free_y") + 
    scale_alpha_manual(values=alpha_vec, guide="none") +
    scale_size_manual(values=size_vec, guide="none") +
    scale_shape_manual(values=shape_vec, guide="none") + 
    ylab("Schizophrenia GWAS -log10 P value ") + xlim(xrange/1e6) + 
    xlab("") + theme_bw() + 
    geom_text_repel(data=datP[datP$type!="SNP",], size=3) + 
    geom_vline(data=snp.line, aes(xintercept = pos), linetype="dotted") 
  #theme(legend.position = c(0.2,0.6))
  
  p2 <- ggplot(dat.ct, aes(x=pos, y=susie_pip, shape=type, color=type, size=type, alpha=type, label=label)) + 
    geom_point() +
    #  facet_grid(celltype~., scales = "free_y") + 
    scale_alpha_manual(values = alpha_vec, guide="none") +
    scale_size_manual(values = size_vec, guide="none") +
    scale_shape_manual(values = shape_vec, guide="none") +
    ylab(paste0("cTWAS PIP (",celltype,")")) + xlim(xrange/1e6) + ylim(0,1) + 
    xlab("") + theme_bw() + 
    geom_text_repel(data=dat.ct[dat.ct$type!="SNP",], size=3) + 
    geom_vline(data=snp.line, aes(xintercept = pos), linetype="dotted") 
    
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  
  g <- rbind(g1, g2, size = "first")
  pdf(paste0("caQTL_cTWAS_",snp,"_", celltype,".pdf"), width=10, height=6)
  grid.newpage()
  grid.draw(g)
  dev.off()
}

focus <- "chr20_46051914_46052414"
assoc <- paste0(focus,"|",c("1hr_npglut", "6hr_npglut", "6hr_nmglut"))
locus_plot_caqtl("rs12624433", assoc ,focus, "npglut")
locus_plot_caqtl("rs12624433", assoc, focus, "nmglut")

focus <- "chr17_19317206_19317706"
assoc <- paste0(focus,"|", c("1hr_npglut", "6hr_npglut", "6hr_nmglut", "1hr_GABA", "6hr_GABA"))
locus_plot_caqtl("rs1472932", assoc ,focus, "npglut")
locus_plot_caqtl("rs1472932", assoc, focus, "nmglut")
locus_plot_caqtl("rs1472932", assoc, focus, "GABA")

focus <- c("chr8_143939404_143939904", "chr8_143933798_143934298", "chr8-143940066-143940566")
assoc <- c("chr8_143933798_143934298|0hr_npglut", "chr8_143939404_143939904|1hr_GABA", "chr8_143940066_143940566|1hr_GABA",
           "chr8_143933798_143934298|1hr_nmglut", "chr8_143933798_143934298|1hr_npglut", "chr8_143939404_143939904|6hr_GABA",
           "chr8_143940066_143940566|6hr_GABA", "chr8_143948394_143948894|6hr_GABA", "chr8_143933798_143934298|6hr_nmglut",
           "chr8_143933798_143934298|6hr_npglut")
locus_plot_caqtl("rs11782331", assoc ,focus, "npglut")
locus_plot_caqtl("rs11782331", assoc, focus, "nmglut")
locus_plot_caqtl("rs11782331", assoc, focus, "GABA")


### Plotting overlap with dynamic ASoC
focus <- "chr7_2010548_2011048"
assoc <- c(paste0("chr7_2010548_2011048|", c(0,1,6), "hr_GABA"))
locus_plot_caqtl("rs11771625", assoc, focus, "GABA")

focus <- "chr10_86244301_86244801"
assoc <- c(paste0(focus,"|", c(0,1,6), "hr_npglut"))
locus_plot_caqtl("rs74716128", assoc, focus, "npglut")

focus <- "chr18_29920961_29921461"
assoc <- c(paste0(focus,"|", c(0,1,6), "hr_npglut"))
locus_plot_caqtl("rs9964412", assoc, focus, "npglut")

focus <- "chr2_172098931_172099431"
assoc <- c(paste0(focus,"|", c(0,1,6), "hr_nmglut"))
locus_plot_caqtl("rs172099203", assoc, focus, "nmglut")




