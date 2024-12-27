# Chuxuan Li 02/25/2022
# Use MAST to do differential expression analysis directly on 5-line data

# init ####
library(Seurat)

library(GSEABase)
library(limma)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rsvd)
library(MAST)

library(ggplot2)
library(RColorBrewer)
library(GGally)
library(reshape2)
library(data.table)
library(knitr)
library(stringr)
library(NMF)


# load data and build data.table ####
load("/5line_codes_and_rdata/27Jan2022_5line_added_cellline_idents.RData")

unique(RNAseq_integrated_labeled$cell.line.ident)
sum(RNAseq_integrated_labeled$cell.line.ident == "unmatched") # 25173

pure_human <- subset(RNAseq_integrated_labeled, cell.line.ident %in% c("CD_54", "CD_27",
                                                                       "CD_26", "CD_08",
                                                                       "CD_25"))
pure_human$spec.cell.type <- "others"
pure_human$spec.cell.type[pure_human$broad.cell.type == "NEFM_pos_glut"] <- "NEFM_pos_glut"
pure_human$spec.cell.type[pure_human$broad.cell.type == "NEFM_neg_glut"] <- "NEFM_neg_glut"
pure_human$spec.cell.type[pure_human$broad.cell.type == "NPC"] <- "NPC"
pure_human$spec.cell.type[pure_human$broad.cell.type == "GABA"] <- "GABA"
unique(pure_human$spec.cell.type)

DefaultAssay(pure_human) <- "RNA"
pure_human <- ScaleData(pure_human)


# make object

feat_mat <- data.frame(geneid = pure_human@assays$RNA@counts@Dimnames[[1]], 
                       row.names = pure_human@assays$RNA@counts@Dimnames[[1]])
cell_mat <- pure_human@meta.data # has to be a data.frame
ncol(pure_human@assays$RNA@scale.data)
nrow(cell_mat)
nrow(pure_human@assays$RNA@scale.data)
nrow(feat_mat)
scaRaw <- FromMatrix(exprsArray = pure_human@assays$RNA@scale.data, # count matrix
                     cData = cell_mat, # cell level covariates
                     fData = feat_mat # feature level covariates
)

# filtering
expressed_genes <- freq(scaRaw) > 0.05
sca <- scaRaw[expressed_genes, ]
thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 10, min_per_bin = 20)
par(mfrow = c(5,4))
plot(thres)


# differential expression ####
cond <- factor(colData(sca)$time.ident)
cond <- relevel(cond, "0hr")
colData(sca)$condition <- cond

# For each gene in sca, fits the hurdle model in formula (linear for et>0), 
#logistic for et==0 vs et>0. Return an object of class ZlmFit containing slots 
#giving the coefficients, variance-covariance matrices, etc. After each gene, 
#optionally run the function on the fit named by 'hook'
zlmCond <- zlm(~condition + spec.cell.type, sca) # time vs. , LHS is variable, RHS is predictor

#only test the condition coefficient.
summaryCond <- summary(zlmCond, doLRT = 'condition1hr') 
#print the top 4 genes by contrast using the logFC
print(summaryCond, n = 4)

# make data.table with coeff, standar errors
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='condition1hr' & component=='H',
                            .(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='condition1hr' & component=='logFC', 
                            .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)


# visualization
entrez_to_plot <- fcHurdleSig[1:16, primerid]
genes_to_plot <- fcHurdleSig[1:16, geneid]
flat_dat <- as(sca[entrez_to_plot,], 'data.table')
ggbase <- ggplot(flat_dat, aes(x=condition, y=thresh,color=condition)) + 
  geom_jitter() +
  facet_wrap(~symbolid, scale='free_y') +
  ggtitle("DE Genes in Activated MAIT Cells")
ggbase +
  geom_violin() 

mat_to_plot <- assay(sca[entrez_to_plot, ])
rownames(mat_to_plot) <- genes_to_plot
aheatmap(mat_to_plot,
         annCol = colData(sca)[ ,"condition"],
         main = "DE genes",
         col = rev(colorRampPalette(colors = brewer.pal(name = "PiYG", n = 10))(20)))
