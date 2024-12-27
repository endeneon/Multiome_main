### Pi1 analysis estimating the non-null proportion of neuron activity eQTL
### in GTEx brain tissues

library(qvalue)
library(data.table)
library(dplyr)
#test <- read.table("Documents/neuron_stim/GTEX_eqtl/GTEx_Analysis_v8_eQTL_independent/Brain_Amygdala.v8.independent_eqtls.txt.gz", header=T)

setwd("/project/xinhe/lifan/neuron_stim/")

genes <- read.table("symbol2ENSG.txt", header=T)
snpinfo <- fread("/project2/xinhe/lifanl/neuron_stim/geno_info_100lines.txt")

gtex.store <- "/scratch/midway3/lifanl/"
gtexfs <- list.files(gtex.store, pattern = "Brain")
gtex.regions <- substr(gtexfs, 1, nchar(gtexfs)-16)

eqtl.store <- "mateqtl_100lines_output/final_eqtl_08222023/"
eqtl.files <- list.files(eqtl.store,pattern="eqtl.rds")
eqtls <- list()
for(f in eqtl.files[-10]) {
  n <- substr(f, 9, nchar(f)-22)
  temp <- readRDS(paste0(eqtl.store,f))$cis$eqtls
  eqtls[[n]] <- temp[temp$FDR<0.05,]
}

# Transform matrix eQTL summary statistics into a format compatible with GTEx
mapping <- merge(data.table(aggr.eqtl[,c("SNP","gene","beta")]), genes, 
                 by.x="gene", by.y="hgnc_symbol", all.x=T, all.y=F)
mapping1 <- merge(mapping, snpinfo, by.x="SNP", by.y="V5", all.x=T, all.y=F)
mapping1$snpid <- paste(mapping1$V1, mapping1$V2, mapping1$V3, mapping1$V4, "b38", sep="_")
sig.map <- mapping1[,c("ensembl_gene_id","snpid","beta")]

res <- list()
for(r in gtex.regions) {
  gtex1 <- fread(paste0(gtex.store, r,".allpairs.txt.gz"))
  gtex1$gene <- nth(tstrsplit(gtex1$gene, split ="\\."),n=1)
  inters <- merge(sig.map, gtex1[,c("gene","variant_id","pval_nominal")],
                  by.x=c("ensembl_gene_id","snpid"), by.y=c("gene","variant_id"), all.y=F)
  res[[r]] <- 1 - qvalue(inters$pval_nominal)$pi0
}

gtex.store <- "/scratch/midway3/lifanl/"
gtexfs <- list.files(gtex.store, pattern = "Brain")
test <- fread(paste0(gtex.store,gtexfs[6]), header=T, stringsAsFactors = F)
test$gene <- nth(tstrsplit(test$gene, split ="\\."),n=1)
dym_GABA$id <- paste(dym_GABA$gene, dym_GABA$snps, sep="|")
genes <- read.table("symbol2ENSG.txt", header=T)
snpinfo <- fread("/project2/xinhe/lifanl/neuron_stim/geno_info_100lines.txt")
head(snpinfo)
pi1.res <- list()

sig_dym <- dym_nmglut[dym_nmglut$FDR<0.05,]
mapping <- merge(data.table(sig_dym[,c("snps","gene","FDR")]), genes, 
                 by.x="gene", by.y="hgnc_symbol", all.x=T, all.y=F)
mapping1 <- merge(mapping, snpinfo, by.x="snps", by.y="V5", all.x=T, all.y=F)
mapping1$snpid <- paste(mapping1$V1, mapping1$V2, mapping1$V3, mapping1$V4, "b38", sep="_")
sig.map2 <- mapping1[,c("ensembl_gene_id","snpid")]
inters <- merge(sig.map2, test[,c("gene","variant_id","pval_nominal")],
                by.x=c("ensembl_gene_id","snpid"), by.y=c("gene","variant_id"), all.y=F)
pi1.res[["nmglut"]] <- 1 - qvalue(inters$pval_nominal)$pi0

sig_dym <- dym_npglut[dym_npglut$FDR<0.05,]
mapping <- merge(data.table(sig_dym[,c("snps","gene","FDR")]), genes, 
                 by.x="gene", by.y="hgnc_symbol", all.x=T, all.y=F)
mapping1 <- merge(mapping, snpinfo, by.x="snps", by.y="V5", all.x=T, all.y=F)
mapping1$snpid <- paste(mapping1$V1, mapping1$V2, mapping1$V3, mapping1$V4, "b38", sep="_")
sig.map2 <- mapping1[,c("ensembl_gene_id","snpid")]
inters <- merge(sig.map2, test[,c("gene","variant_id","pval_nominal")],
                by.x=c("ensembl_gene_id","snpid"), by.y=c("gene","variant_id"), all.y=F)
pi1.res[["npglut"]] <- 1 - qvalue(inters$pval_nominal)$pi0

pi1.cond <- read.table("/project/xinhe/lifan/neuron_stim/pi1_with_gtex.txt")
pi1.cortex <- pi1.cond[c(1,4,7,2,5,8,3,6,9),'Brain_Cortex']
pi1.cortex <- c(pi1.cortex[1:3], pi1.res[["GABA"]], 
                pi1.cortex[4:6], pi1.res[["nmglut"]],
                pi1.cortex[7:9], pi1.res[["npglut"]])
names(pi1.cortex) <- c(paste0(c("0hr","1hr","6hr","dynamic"),"_GABA"),
                       paste0(c("0hr","1hr","6hr","dynamic"),"_nmglut"),
                       paste0(c("0hr","1hr","6hr","dynamic"),"_npglut"))
library(ggplot2)
dat <- data.frame(pi1.cortex)
dat$celltype <- 
  par(mar=c(4,8,2,2))
barplot(pi1.cortex, las=2, col=c(rep("red",4),rep("green",4),rep("blue",4)), horiz = T, xlab="Pi1 with GTEx_Cortex")



