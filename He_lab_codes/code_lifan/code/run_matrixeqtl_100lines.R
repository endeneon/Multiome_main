library(MatrixEQTL)
library(gap)
library(optparse)

option_list <- list(
  make_option(
    "--input_folder", action = "store", default = NA, type = 'character',
    help = "Folder where the input dataset is [required]"
  ),
  make_option(
    "--cond", action = "store", default = "0hr_GABA", type = 'character',
    help = "The split index. Default is %default [optional]"
  ),
  make_option(
    "--cisdist", action = "store", default = "1e6", type = 'character',
    help = "Upstream and downstream distance to the gene body. Default is %default [optional]"
  ),
  make_option(
    "--genoPC", action = "store", default = 5, type = 'integer',
    help = "The number of genotype PCs as covariates. Default is %default [optional]"
  ),
  make_option(
    "--exprPC", action = "store", default = 10, type = 'integer',
    help = "The number of gene expression PCs as covariates. Default is %default [optional]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
ct <- opt$cond # from 1 to 12
folder <- opt$input_folder
cisDist <- as.numeric(opt$cisdist)
NgenoPC <- opt$genoPC
NexprPC <- opt$exprPC
DEST <- "/project/xinhe/lifan/neuron_stim/mateqtl_100lines_output/final_eqtl_08222023/"

useModel = modelLINEAR
SNP_file_name = paste0(folder,"snps.txt")
snps_location_file_name = paste0(folder, "snpsloc.txt")
gene_location_file_name = paste0(folder, "geneloc.txt")

errorCovariance = numeric()
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0; # Skip trans eQTL
# 
# celltypes = c("nmglut", "npglut", "GABA")
# timepoints = c("0hr", "1hr", "6hr")
# cellbytime <- expand.grid(timepoints, celltypes)
# cellbytime <- paste(cellbytime$Var1, cellbytime$Var2, sep="_")

geno = as.matrix(read.table(SNP_file_name, header = TRUE, stringsAsFactors = FALSE))

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFact=F);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFact=F);

# library(pheatmap)
# cvrt_1hr_GABA <- as.matrix(read.table(paste0(folder, "covar_1hr_GABA.txt"), header=T, row.names=1))
# cvrt_6hr_GABA <- as.matrix(read.table(paste0(folder, "covar_6hr_GABA.txt"), header=T, row.names=1))
# cvrt_6hr_npglut <- as.matrix(read.table(paste0(folder, "covar_6hr_npglut.txt"), header=T, row.names=1))
# pheatmap(cor(cvrt_1hr_GABA[,6:33], cvrt_1hr_GABA[,44:63]),cluster_rows = F, cluster_cols = F)
# pheatmap(cor(cvrt_6hr_GABA[,6:33], cvrt_6hr_GABA[,44:63]),cluster_rows = F, cluster_cols = F)
# pheatmap(cor(cvrt_6hr_npglut[,6:33], cvrt_6hr_npglut[,44:63]),cluster_rows = F, cluster_cols = F)

# ct.geno <- expand.grid(cellbytime, 0:1, 1:3)
# colnames(ct.geno) <- c("cellbytime","genoPC","exprPC")
# ct.geno$Negenes <- 0
# i <- 1
# ct <- ct.geno[i,1]
# cat("Running eQTL for", ct)

full.cvrt = as.matrix(read.table(paste0(folder, "covar_", ct, ".txt"), header=T, row.names=1))
#drop <- c(paste0("genoPC",6:10),paste0("exprPC",6:10))
#cvrt1 <- full.cvrt[, !colnames(full.cvrt) %in% drop] 

if (NgenoPC < 10) {
  drop <- paste0("genoPC", (NgenoPC+1):10)
  cvrt1 <- full.cvrt[, !colnames(full.cvrt) %in% drop] 
} else {
  cvrt1 <- full.cvrt
}

if (NexprPC < 20) {
  drop <- paste0("exprPC", (NexprPC+1):20)
  cvrt1 <- cvrt1[, !colnames(cvrt1) %in% drop] 
}

cvrt1 <- cvrt1[,apply(cvrt1,2,sd)>0] # Remove empty categories
#cvrt1 <- cvrt1[, !startsWith(colnames(cvrt1), "seq.batch")] # Removing sequencing batch because it correlates strongly with expression PC
#cvrt1 <- cvrt1[, -c(4,6)]  # Remove one coculture batch so they will not become colinear, also remove cell count and use proportion only
#cvrt1 <- cvrt1[, -4] # Remove cell count
cvrt1 <- cvrt1[, !grepl("batch",colnames(cvrt1))] # Remove all coculture batch for testing purpose.

expression_file_name = paste0(folder, "expr_", ct, ".txt")
expr = as.matrix(read.table(expression_file_name, header=T, row.names=1))
#expr1 = normalize.quantiles(log(expr+1))
#colnames(expr1) <- colnames(expr)
#rownames(expr1) <- rownames(expr)

snp1 <- geno[, colnames(expr)]

output_file_name_tra = paste0(DEST, ct, "_genoPC" , NgenoPC, "_exprPC", NexprPC, "_cis", opt$cisdist, "_eqtl.txt")
output_file_name_cis = paste0(DEST, ct, "_genoPC" , NgenoPC, "_exprPC", NexprPC, "_trans", opt$cisdist, "_eqtl.txt")

gene = SlicedData$new();
gene$CreateFromMatrix(expr);

cvrt = SlicedData$new();
cvrt$CreateFromMatrix(t(cvrt1))

snps = SlicedData$new()
snps$CreateFromMatrix(snp1)
snps$ResliceCombined(sliceSize = 1200L)

res = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

saveRDS(res, paste0(DEST, ct, "_rmCellCount", "_genoPC" , NgenoPC, "_exprPC", NexprPC, "_cis", opt$cisdist, "_eqtl.rds"))
