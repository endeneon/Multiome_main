# Chuxuan Li 11/04/2022
# Still using MAST directly, but using raw count -> log2(TPM + 1) -> MAST to 
#see if it is different from using SCTransform to normalize the counts

# init ####
library(ggplot2)
library(ggrepel)
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
library(Seurat)
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000_with_covars_from_30Sep2022_single_cell_DE.RData")

# make object ####
unique(obj$cell.line.ident)
unique(obj$cell.type)
npglut <- subset(obj, cell.type == "NEFM_neg_glut")

DefaultAssay(npglut) <- "RNA"

# calculate TPM ####
library(readr)
library(dplyr)
gene_length <- read_table("/data/Databases/Genomes/CellRanger_10x/refdata-cellranger-arc-GRCh38-2020-A/genes/gene_length.txt", 
                        col_names = c("start", "end", "gene_name"))
gene_length$length <- gene_length$end - gene_length$start + 1
gene_length <- as.data.frame(distinct_at(gene_length, c("gene_name"), .keep_all = T))
rownames(gene_length) <- unique(gene_length$gene_name)
gene_length["BDNF", 4]
count_dgcmat <- npglut@assays$RNA@counts
convertToMat <- function(row) {
  return(as.numeric(row))
}
library(parallel)
cl <- makeCluster(2)
count_mat <- parApply(cl = cl, X = count_dgcmat, MARGIN = 1, FUN = convertToMat)
stopCluster(cl)
count_mat <- rbind(count_mat, colnames(count_mat))
RPK <- function(col, gl) {
  rowname <- col[length(col)]
  length <- gl[rowname, 4]
  rpk <- as.numeric(col[1:(length(col) - 1)]) / length
  return(rpk)
}
cl <- makeCluster(2)
RPK_mat <- as.matrix(parApply(cl = cl, X = count_mat, MARGIN = 2, FUN = RPK, gl = gene_length))
stopCluster(cl)
libRPK <- sum(rowSums(RPK_mat))
TPM_mat <- RPK_mat / libRPK
save(TPM_mat, file = "18line_TPM_matrix.RData")

feat_mat <- data.frame(geneid = npglut@assays$RNA@counts@Dimnames[[1]], 
                       row.names = npglut@assays$RNA@counts@Dimnames[[1]])
cell_mat <- npglut@meta.data # has to be a data.frame, contains metadata for each cell
ncol(npglut@assays$RNA@counts) # 41536 cells
nrow(cell_mat) # 41536
nrow(npglut@assays$RNA@counts) # 29546 genes
nrow(feat_mat) # 29546
scaRaw <- FromMatrix(exprsArray = as.array(npglut@assays$RNA@counts), # count matrix
                     cData = cell_mat, # cell level covariates
                     fData = feat_mat, # feature level covariates
                     check_sanity = F
)


save(scaRaw, file = "./MAST_direct_DE/scaRaw_obj_npglut.RData")
aheatmap(x = assay(scaRaw[1:1000,]), 
         labRow = '', 
         annCol = as.data.frame(colData(scaRaw)[,c('time.ident', 'orig.ident', 'cell.type')]), 
         distfun = 'spearman')

# filtering
expressed_genes <- freq(scaRaw) > 0.1 # proportion of nonzero values across all cells
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
zlmCond <- zlm(~condition + age + sex + aff + orig.ident + cell.type.fraction, sca) # time vs. , LHS is variable, RHS is predictor
save(zlmCond, file = "./MAST_direct_DE/MAST_directly_zlmCond_npglut.RData")

# 1hr ####
#only test the condition coefficient.
summaryCond_1v0 <- summary(zlmCond, doLRT = 'condition1hr') 
save(summaryCond_1v0, file = "./MAST_direct_DE/summaryCond_1v0_npglut.RData")

#print the top 4 genes by contrast using the logFC
print(summaryCond_1v0, n = 4)

# make data.table with coeff, standar errors
summaryDt <- summaryCond_1v0$datatable
fcHurdle_1v0 <- merge(summaryDt[contrast=='condition1hr' & component=='H',
                                .(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='condition1hr' & component=='logFC', 
                                .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle_1v0[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
fcHurdleSig_1v0 <- merge(fcHurdle_1v0[fdr < 0.05], 
                         as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig_1v0, fdr)
save(fcHurdle_1v0, file = "./MAST_direct_DE/result_1v0_npglut.RData")

# 6hr ####
#only test the condition coefficient.
summaryCond_6v0 <- summary(zlmCond, doLRT = 'condition6hr') 
save(summaryCond_6v0, file = "./MAST_direct_DE/summaryCond_6v0_npglut.RData")

#print the top 4 genes by contrast using the logFC
print(summaryCond_6v0, n = 4)

# make data.table with coeff, standar errors
summaryDt <- summaryCond_6v0$datatable
fcHurdle_6v0 <- merge(summaryDt[contrast=='condition6hr' & component=='H',
                                .(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='condition6hr' & component=='logFC', 
                                .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle_6v0[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
fcHurdleSig_6v0 <- merge(fcHurdle_6v0[fdr < 0.05], 
                         as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig_6v0, fdr)
save(fcHurdle_6v0, file = "./MAST_direct_DE/result_6v0_npglut.RData")

