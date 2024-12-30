### Plotting based on cTWAS 1.0.4 results

library(ggplot2)
library(grid)
library(AnnotationHub)
library(ggrepel)
library(patchwork)
library(gtable)
library(data.table)
library(locuszoomr)
library(latex2exp)

ah <- AnnotationHub(localHub=T)
ensdb <- ah[["AH104864"]] #EnsDb.Hsapiens.v106
#results_dir <- "/project/xinhe/lifan/neuron_stim/ctwas_SCZ/results/cellgroup_SCZ/"
results_dir <- "/folder/with/multigroup/cTWAS/results"
analysis_id <- c("GABA","npglut","nmglut")
setwd("/working/directory")

### Load cTWAS finemapping results
ctwas.res <- read.table(paste0(results_dir, "joint_celltype_cover9",".susieIrss.txt"), header=T)
ctwas.res$susie_pip[ctwas.res$cs_index==0] <- 0

### Information about SNPs
snps <- fread("../PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv")


### Load Imputed gene Z score
load.zgene <- function(celltype) {
  ld_exprfs <- paste0(results_dir, "_chr", 1:22, ".expr.gz")
  z_gene <- list()
  for (i in 1:22){
    load(paste0(results_dir,"joint_celltype", "_chr", i, ".exprqc.Rd"))
    z_gene[[i]] <- z_gene_chr
  }
  z_gene <- do.call(rbind, z_gene)
  #rm(qclist, wgtlist, z_gene_chr)
  #rownames(z_gene) <- NULL
  z_gene
}
zgene.all <- load.zgene("joint_celltype")

res <- ctwas.res[(ctwas.res$cs_index>0) & (ctwas.res$type!="SNP"),]
cond <- names(table(res$type))
res$genename <- sapply(strsplit(res$id, "\\|"), function(x) {x[1]})
res2 <- aggregate(res$susie_pip, list(res$genename), sum)

#genes <- c("cPeak_chr7_1986364_1986864")
#near.gene <- "MAD1L1"

genes <- c("cPeak_chr12_57096240_57096740")
near.gene <- "STAT6"

#genes <- c("ASOC_rs2157591")
#near.gene <- "TCF20"

#genes <- c("cPeak_chr5_90929346_90929846")
#near.gene <- "ADGRV1"

for(focus in genes) {
  temp <- ctwas.res[ctwas.res$id %in% paste0(focus,"|",cond),]
  
  chr = unique(temp$region_tag1)
  tag.count <- table(temp$region_tag2)
  tag = names(tag.count)[order(tag.count,decreasing=T)[1]]
  
  data(SLE_gwas_sub)  ## Inplace dataset to construct genetracks
  loc <- locus(xrange=c(temp$pos[1]-2.5e5,temp$pos[1]+2.5e5), seqname=chr, data=SLE_gwas_sub, ens_db=ensdb)
  loc$xrange <- c(loc$xrange[1], loc$xrange[length(loc$xrange)]) ## Sometime duplicate genes yields duplicate xranges
  ctwas.sub <- ctwas.res[(ctwas.res$region_tag1==chr) & (ctwas.res$pos>=loc$xrange[1]) & (ctwas.res$pos<=loc$xrange[2]),]
  
  # Format the dataframe for plotting
  
  construct.dat <- function(res.sub) {
    dat <- res.sub
    dat$name <- sapply(strsplit(dat$id,"\\|"), function(x) {x[1]})
    dat$time_point <- sapply(strsplit(dat$type,"_"), function(x) {x[1]})
    # dat$celltype <- factor(dat$celltype, levels=c("GABA","npglut","nmglut"))
    dat$pos <- dat$pos / 1e6
    genes <- unique(dat[dat$type!="SNP","name"])
    loc$TX <- loc$TX[loc$TX$gene_name %in% genes,]
    dat$label <- ""
    #dat$label[(dat$susie_pip>0.1) & (dat$type!="SNP")] <- paste0(dat$name[(dat$susie_pip>0.1) & (dat$type!="SNP")],"(",dat$type,")")
    dat
  }
  
  datAll <- construct.dat(ctwas.sub)
  datAll$cell_type <- "SNP"
  datAll$cell_type[datAll$type!="SNP"] <- sapply(strsplit(datAll$type[datAll$type!="SNP"],"_"), function(x){x[2]})
  datAll$label[datAll$type!="SNP"] <- paste0(datAll$name,"(",datAll$time_point,",",datAll$cell_type,")")[datAll$type!="SNP"]
  
  # Obtain position information from cTWAS results.
  z_gene2 <- merge(zgene.all, datAll, by="id", all.x=F, all.y=F)
  z_gene2$PVAL <- pnorm(abs(z_gene2$z), lower.tail = F)
  z_gene2$logP <- -log10(z_gene2$PVAL)
  colnames(z_gene2)[1] <- "ID"
  
  snp.sub <- snps[snps$CHROM==chr,c("ID","PVAL")]
  snp.sub1 <- merge(snp.sub, datAll[,1:3], by.x="ID", by.y="id", all.x=F, all.y=F)
  snp.sub1$logP <- -log10(snp.sub1$PVAL)
  snp.sub1$type <- "SNP"
  snp.sub1$time_point <- "SNP"
  snp.sub1$cell_type <- "SNP"
  snp.sub1$label <- ""
  datP <- rbind(snp.sub1, z_gene2[,c("ID","PVAL","chrom","pos","logP","type","time_point","cell_type","label")])
  #datP$label[startsWith(datP$ID,focus)] <- datP$ID[startsWith(datP$ID,focus)]
  
  colnames(datAll)[11] <- "color"
  colnames(datAll)[13] <- "shape"
  colnames(datP)[7] <- "color"
  colnames(datP)[8] <- "shape"
  datP$color[datP$color!="SNP"] <- paste0(datP$color[datP$color!="SNP"],"_","gene")
  datP$shape[datP$shape!="SNP"] <- paste0(datP$shape[datP$shape!="SNP"],"_","gene")
  datAll$color[datAll$color!="SNP"] <- paste0(datAll$color[datAll$color!="SNP"],"_","gene")
  datAll$shape[datAll$shape!="SNP"] <- paste0(datAll$shape[datAll$shape!="SNP"],"_","gene")
  
  ## Construct Z score locus plot
  alpha_vec <- rep(1,length(unique(datAll$type)))
  names(alpha_vec) <- unique(datAll$type)
  alpha_vec["SNP"] <- 0.5
  
  size_vec <- rep(4,length(unique(datAll$type)))
  names(size_vec) <- unique(datAll$type)
  size_vec["SNP"] <-1
  
  shape_vec <- c(17,15,18,16)
  names(shape_vec) <- c("GABA_gene","nmglut_gene","npglut_gene","SNP")
  
  color_vec <- c("orange","red","purple","blue")
  names(color_vec) <- c("0hr_gene","1hr_gene","6hr_gene","SNP")
  
  datP$label[datP$logP<3] <- ""
  datAll$label[!(datAll$label %in% datP$label)] <- ""
  
  p1 <- ggplot(datAll, aes(x=pos, y=susie_pip, shape=shape, color=color, size=type, alpha=type, label=label)) + 
    geom_point() + #geom_vline(xintercept = temp$pos[1]/1e6, linetype="dashed", color = "black", size=1, alpha=0.7) + 
    #  facet_grid(celltype~., scales = "free_y") + 
    #  annotate("text", temp$pos[1]/1e6, 1, hjust = -.25, 
    #           label = paste(near.gene,"eQTL")) +
    scale_alpha_manual(values=alpha_vec, guide="none") +
    scale_size_manual(values=size_vec, guide="none") +
    scale_shape_manual(values=shape_vec, guide="none") + 
    scale_color_manual(values=color_vec) + 
    ylab("cTWAS PIP") + xlim(loc$xrange/1e6) + ylim(0,1) + 
    xlab("") + theme_bw() +
    geom_label_repel(data=datAll[datAll$type!="SNP",], size=4) + 
    theme(legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 8))
  
  p2 <- ggplot(datP, aes(x=pos, y=logP, shape=shape, color=color, size=type, alpha=type, label=label)) + 
    geom_point() + guides(color="none") + 
    #  geom_vline(xintercept = temp$pos[1]/1e6, linetype="dashed", color = "black", size=1) + 
    #  annotate("text", temp$pos[1]/1e6, max(datP$logP), hjust = -.25, 
    #           label = paste(near.gene,"eQTL")) +
    # facet_grid(celltype~., scales = "free_y") + 
    scale_alpha_manual(values=alpha_vec, guide="none") +
    scale_size_manual(values=size_vec, guide="none") +
    scale_shape_manual(values=shape_vec) + 
    scale_color_manual(values=color_vec, guide="none") + 
    
    ylab(TeX(r'(GWAS $-log_{10}($p-value$)$)')) + xlim(loc$xrange/1e6) + 
    xlab("") + theme_bw() + 
    geom_label_repel(data=datP[datP$type!="SNP",], size=2) + 
    theme(legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 8))
  
  ## Adding gene tracks
  height = unit(3, "cm") # Heights of the track
  g1 <- ggplotGrob(p1)
  g2 <- ggplotGrob(p2)
  g <- rbind(g2,g1, size = "first")
  panels_extent <- g %>% find_panel()
  pg <- g %>%
    gtable_add_rows(heights = height) %>%
    gtable_add_grob(gg_genetracks(loc, filter_gene_biotype = "protein_coding"),
                    t = -1, b = -1,
                    l = panels_extent$l, r = panels_extent$l +1)
  pg$widths <- unit.pmax(g2$widths, g1$widths)
  pdf(paste0("ctwas_locus_plots/",focus,"_PIPsum_",sprintf("%.2f",res2[res2$Group.1==focus,"x"]),near.gene,"_eQTL","_joint_celltype.pdf"),width=8,height=5)
  grid.newpage()
  grid.draw(pg)
  dev.off()
}
