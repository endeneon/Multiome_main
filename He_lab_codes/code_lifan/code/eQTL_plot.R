### Plot examples of eqtls
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggbio)
library(Homo.sapiens)
library(cowplot)

setwd("/Users/lifanliang/Documents/neuron_stim/mateqtl_input_100lines/")

gwas <- fread("../SCZ_GWAS/PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv")
finemap <- read.table("../SCZ_GWAS/S11d_finemapping.txt",header=T)
snps <- fread("snps.txt")
snpinfo <- fread("../geno_info_100lines.txt")
colnames(snpinfo) <- c("chr","pos","ref","alt","rsid")
eqtl.store <- "../final_eqtl_08222023/"
eqtl.files <- list.files(eqtl.store,pattern=".rds")
eqtls <- list()
for(f in eqtl.files) {
  n <- substr(f, 1, nchar(f)-47)
  eqtls[[n]] <- readRDS(paste0(eqtl.store,f))
}
eqtls1 <- lapply(eqtls, function(x){x$cis$eqtls})

### Retain significant eQTLs
eqtls2 <- lapply(eqtls1, function(x){x[x$FDR<0.05,]})
for(n in names(eqtls2)) {
  write.table(eqtls2[[n]],paste0("../supp_table/TableS21_signficant_eQTL_",n,".txt"), quote=F,row.names = F)
}
Negenes <- sapply(eqtls2,function(x){length(unique(x$gene))})
#write.table(Negenes, "../supp_table/TableS21_egene_count.txt", quote=F, col.names=F)

### Retain 
target <- "CPT1C"
eqtls2 <- lapply(eqtls1, function(x){head(x[x$gene==target,],1)})
cpt1c <- do.call(rbind, eqtls2)

expr.files <- list.files(pattern="^expr_")
exprs <- list()
for(f in expr.files) {
  n <- substr(f, 6, nchar(f)-4)
  exprs[[n]] <- read.table(f)
}

  ### For box plot
plot.qtl <- function(celltype="GABA",gene,snp) {
  gwasP <- gwas[gwas$ID==snp,]$PVAL
  #PIP <- finemap[finemap$rsid==snp,"finemap_marginal_posterior_inclusion_probability"]
  center <- gwas[gwas$ID==snp,c("CHROM","POS")]
  
  expr0 <- exprs[[paste0("0hr_",celltype)]][gene,]
  expr1 <- exprs[[paste0("1hr_",celltype)]][gene,]
  expr6 <- exprs[[paste0("6hr_",celltype)]][gene,]
  #resid0 <- lm(`gene`~.,cbind(t(expr0), covars1[[paste0("0hr_",celltype)]][-1]))$residual
  #resid1 <- lm(`gene`~.,cbind(t(expr1), covars1[[paste0("1hr_",celltype)]][-1]))$residual
  #resid6 <- lm(`gene`~.,cbind(t(expr6), covars1[[paste0("6hr_",celltype)]][-1]))$residual
  eqtl0 <- eqtls[[paste0("0hr_",celltype)]]$cis$eqtl
  eqtl1 <- eqtls[[paste0("1hr_",celltype)]]$cis$eqtl
  eqtl6 <- eqtls[[paste0("6hr_",celltype)]]$cis$eqtl
  beta0 <- eqtl0[(eqtl0$snps==snp)&(eqtl0$gene==gene),"beta"]
  beta1 <- eqtl1[(eqtl1$snps==snp)&(eqtl1$gene==gene),"beta"]
  beta6 <- eqtl6[(eqtl6$snps==snp)&(eqtl6$gene==gene),"beta"]
  pval0 <- eqtl0[(eqtl0$snps==snp)&(eqtl0$gene==gene),"pvalue"]
  pval1 <- eqtl1[(eqtl1$snps==snp)&(eqtl1$gene==gene),"pvalue"]
  pval6 <- eqtl6[(eqtl6$snps==snp)&(eqtl6$gene==gene),"pvalue"]
  
  #if(length(PIP)==0) PIP <- NA
  #else PIP <- formatC(PIP, digits=2)
  info <- snpinfo[snpinfo$rsid==snp,]
  indiv0 <- colnames(expr0)
  indiv1 <- colnames(expr1)
  indiv6 <- colnames(expr6)
  dosage <- c(as.character(snps[snps$V1==snp,..indiv0]),as.character(snps[snps$V1==snp,..indiv1]),as.character(snps[snps$V1==snp,..indiv6]))
  dosage <- sub("0", paste0(info$ref,info$ref), dosage)
  dosage <- sub("1", paste0(info$ref,info$alt), dosage)
  dosage <- sub("2", paste0(info$alt,info$alt), dosage)
  #print(dosage)
  dat <- data.frame(snp = dosage,
                    expr = as.numeric(cbind(expr0,expr1,expr6)),
                    condition = c(rep(paste0("0hr_",celltype," | beta:",formatC(beta0,digits=2)," | P:",formatC(pval0,digits=2)),ncol(expr0)),
                                  rep(paste0("1hr_",celltype," | beta:",formatC(beta1,digits=2)," | P:",formatC(pval1,digits=2)),ncol(expr1)), 
                                  rep(paste0("6hr_",celltype," | beta:",formatC(beta6,digits=2)," | P:",formatC(pval6,digits=2)),ncol(expr6))))
  dat$snp <- factor(dat$snp, levels=c(paste0(info$ref,info$ref),paste0(info$ref,info$alt),paste0(info$alt,info$alt)))
  #print(dat$snp)
  p2 <- ggplot(dat, aes(y=expr, x=snp)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~condition) + 
    xlab(snp) + ylab(gene) + geom_jitter(color="blue",alpha=0.6, width=0.2) + theme_bw() +
    geom_smooth(method='lm')
    labs(title=paste0("Schizophrenia GWAS P value: ", formatC(gwasP,digits=2)))
  #grid.arrange(arrangeGrob(p1,p2, ncol=1, nrow=2))
  p2
}

plot.qtl("GABA","CROT","rs13233308")
plot.qtl("npglut","ADAM10","rs6494026")
plot.qtl("nmglut","CPT1C","rs12104272")

