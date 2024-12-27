# Chuxuan Li 02/24/2022 
# MAST differential expression tutorial

suppressPackageStartupMessages({
  library(ggplot2)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
})
#options(mc.cores = detectCores() - 1) #if you have multiple cores to spin
options(mc.cores = 1)
knitr::opts_chunk$set(message = FALSE,error = FALSE,warning = FALSE,cache = FALSE,fig.width=8,fig.height=6)

# load data ####

# define constants 
freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)

# load data: matrix of aligned reads from T cell experiment with stimulation
data(maits, package='MAST')
dim(maits$expressionmat) # 96 cells, 16302 genes
head(maits$cdat)
head(maits$fdat)

# make object
scaRaw <- FromMatrix(t(maits$expressionmat), # count matrix transposed
                     maits$cdat, # cell level covariates
                     maits$fdat # feature level covariates
                     )

# make heatmap
aheatmap(x = assay(scaRaw[1:1000,]), 
         labRow = '', # row label
         # colored bar on top labeling which ones are stim/unstim and passed our filter
         annCol = as.data.frame(colData(scaRaw)[, c('condition', 'ourfilter')]), 
         distfun='spearman' # default distance method
         )

# plot PCs
set.seed(123)
plotPCA <- function(sca_obj){
  projection <- rpca(t(assay(sca_obj)), retx=TRUE, k=4)$x
  colnames(projection)=c("PC1","PC2","PC3","PC4")
  pca <- data.table(projection,  as.data.frame(colData(sca_obj)))
  print(ggpairs(pca, columns=c('PC1', 'PC2', 'PC3', 'libSize', 'PercentToHuman', 'nGeneOn', 'exonRate'),
                mapping=aes(color=condition), upper=list(continuous='blank')))
  invisible(pca)
}

plotPCA(scaRaw)

# filter out non mRNAs and low diversity libraries
filterCrit <- with(colData(scaRaw), 
                   pastFastqc == "PASS" & exonRate > 0.3 & 
                     PercentToHuman > 0.6 & nGeneOn > 4000)
sca <- subset(scaRaw,filterCrit) # 73 cells left

# filter out reads that don't have EntrezGene id
eid <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,
              keys = mcols(sca)$entrez,
              keytype ="GENEID",
              columns = c("GENEID","TXNAME"))
ueid <- unique(na.omit(eid)$GENEID)
sca <- sca[mcols(sca)$entrez %in% ueid,]
# Remove invariant genes
sca <- sca[sample(which(freq(sca) > 0), 6000),] # 6000 genes left
# recalculate cellular detection rate in filtered data
cdr2 <- colSums(assay(sca) > 0) # cells with nonzero expression
qplot(x = cdr2, y = colData(sca)$nGeneOn) + 
  xlab('New CDR') + 
  ylab('Old CDR')

colData(sca)$cngeneson <- scale(cdr2)
plotPCA(sca)

colData(sca)$chiploc <- str_extract(string = colData(sca)$wellKey, 
                                    pattern = "C[0-9][0-9]")
head(colData(sca)$chiploc)


# adaptive thresholding ####

# single cell data is zero-inflated and bimodal
scaSample <- sca[sample(which(freq(sca) > .1), 20), ]
flat <- as(scaSample, 'data.table')
ggplot(flat, aes(x = value)) +
  geom_density() +
  facet_wrap(~symbolid, scale = 'free_y')

# determine threshold for each cell (?)
thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 20, min_per_bin = 30)
par(mfrow = c(5,4))
plot(thres)

# apply threshold to the genes
assays(sca, withDimnames = FALSE) <- list(thresh = thres$counts_threshold, 
                                          tpm = assay(sca))
expressed_genes <- freq(sca) > freq_expressed # include genes in at least 0.2 of the sample
sca <- sca[expressed_genes, ]

# differential expression ####
cond <- factor(colData(sca)$condition)
cond <- relevel(cond, "Unstim")
colData(sca)$condition <- cond

# For each gene in sca, fits the hurdle model in formula (linear for et>0), 
#logistic for et==0 vs et>0. Return an object of class ZlmFit containing slots 
#giving the coefficients, variance-covariance matrices, etc. After each gene, 
#optionally run the function on the fit named by 'hook'
zlmCond <- zlm(~condition + cngeneson, sca) # condition vs. cellular detection rate, LHS is variable, RHS is predictor

#only test the condition coefficient.
summaryCond <- summary(zlmCond, doLRT = 'conditionStim') 
#print the top 4 genes by contrast using the logFC
print(summaryCond, n = 4)

# make data.table with coeff, standar errors
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='conditionStim' & component=='H',
                            .(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='conditionStim' & component=='logFC', 
                            .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)


# visualization
entrez_to_plot <- fcHurdleSig[1:50,primerid]
symbols_to_plot <- fcHurdleSig[1:50,symbolid]
flat_dat <- as(sca[entrez_to_plot,], 'data.table')
ggbase <- ggplot(flat_dat, aes(x=condition, y=thresh,color=condition)) + 
  geom_jitter() +
  facet_wrap(~symbolid, scale='free_y') +
  ggtitle("DE Genes in Activated MAIT Cells")
ggbase +
  geom_violin() 

mat_to_plot <- assay(sca[entrez_to_plot, ])
rownames(mat_to_plot) <- symbols_to_plot
aheatmap(mat_to_plot,
         annCol = colData(sca)[ ,"condition"],
         main = "DE genes",
         col = rev(colorRampPalette(colors = brewer.pal(name = "PiYG", n = 10))(20)))
