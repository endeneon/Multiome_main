# Chuxuan Li 07/18/2022
# DE analysis using Limma blocking out all covariates (not cell line, but include
#group) same as those in the data transfer to Lifan in June 2022

# init ####
library(Seurat)
library(limma)
library(edgeR)
library(sva)

library(readr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

{  library(Signac)
  library(Seurat)
  library(EnsDb.Hsapiens.v86)
  library(GenomeInfoDb)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  
  library(ggplot2)
  library(stringr)
}

load("../../../covariates_pooled_together_for_all_linextime_5+18+20.RData")
load("./pseudobulk_DE/res_use_combatseq_mat/by_line_pseudobulk_mat_adjusted.RData")
celllines <- str_remove(colnames(GABA_mat_adj), "_[0|1|6]hr$")

covar_18line <- covar_table[covar_table$cell_line %in% lines, ]
covar_18line <- covar_18line[order(covar_18line$cell_line), ]
rm(covar_table)

GABA_fraction <- covar_18line$GABA_fraction
nmglut_fraction <- covar_18line$nmglut_fraction
npglut_fraction <- covar_18line$npglut_fraction

createDGE <- function(count_matrix, isexpr_val){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = celllines)
  A <- rowSums(y$counts)
  #isexpr <- A > 10
  hasant <- rowSums(is.na(y$genes)) == 0
  #y <- y[isexpr & hasant, , keep.lib.size = F]
  y <- y[hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}
filterByCpm <- function(y, cutoff, nsample) {
  cpm <- cpm(y)
  colInd0 <- seq(1, ncol(y), by = 3)
  colInd1 <- seq(2, ncol(y), by = 3)
  colInd6 <- seq(3, ncol(y), by = 3)
  
  cpm_0hr <- cpm[, colInd0]
  cpm_1hr <- cpm[, colInd1]
  cpm_6hr <- cpm[, colInd6]
  #passfilter <- rowSums(zerohr_cpm >= cutoff) >= nsample
  passfilter <- (rowSums(cpm_0hr >= cutoff) >= nsample |
                   rowSums(cpm_1hr >= cutoff) >= nsample |
                   rowSums(cpm_6hr >= cutoff) >= nsample)
  return(passfilter)
}
cnfV <- function(y, design) {
  y <- calcNormFactors(y)
  # voom (no logcpm)
  v <- voom(y, design, plot = T)
  
  return(v)  
}
contrastFit <- function(fit, design) {
  contr.matrix <- makeContrasts(
    onevszero = time1hr-time0hr, 
    sixvszero = time6hr-time0hr, 
    sixvsone = time6hr-time1hr, levels = colnames(design))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit <- eBayes(fit, trend = T)
  print(summary(decideTests(fit)))
  return(fit)
}


# Limma with voom####
design_GABA <- model.matrix(~0 + time + group + age + sex + disease + GABA_fraction, data = covar_18line) #26353 genes
y_GABA <- createDGE(GABA_mat_adj, 10) # plug in adjusted or raw matrix
ind.keep <- filterByCpm(y_GABA, 1, 9)
sum(ind.keep) #17014 13851 #16977 #17547
v_GABA <- cnfV(y_GABA[ind.keep, ], design_GABA)
fit_GABA <- lmFit(v_GABA, design_GABA)
fit_GABA_contr <- contrastFit(fit_GABA, design_GABA)
plotSA(fit_GABA, main = "Final model: Mean-variance trend")
plotSA(fit_GABA_contr, main = "Final model: Mean-variance trend") # same, just has a line

design_nmglut <- model.matrix(~0 + time + group + age + sex + disease + nmglut_fraction, data = covar_18line) #26353 genes
y_nmglut <- createDGE(nmglut_mat_adj, 10)
ind.keep <- filterByCpm(y_nmglut, 1, 9)
sum(ind.keep) #16503 14366 #16474 #17110
v_nmglut <- cnfV(y_nmglut[ind.keep, ], design_nmglut)
fit_nmglut <- lmFit(v_nmglut, design_nmglut)
fit_nmglut_contr <- contrastFit(fit_nmglut, design_nmglut)
plotSA(fit_nmglut_contr, main = "Final model: Mean-variance trend") 

design_npglut <- model.matrix(~0 + time + group + age + sex + disease + npglut_fraction, data = covar_18line) #29400 genes
y_npglut <- createDGE(npglut_mat_adj, 10)
ind.keep <- filterByCpm(y_npglut, 1, 9)
sum(ind.keep) #17019 15416 #17081 #17557
v_npglut <- cnfV(y_npglut[ind.keep, ], design_npglut)
fit_npglut <- lmFit(v_npglut, design_npglut)
fit_npglut_contr <- contrastFit(fit_npglut, design_npglut)
plotSA(fit_npglut_contr, main = "Final model: Mean-variance trend") 

save(v_GABA, v_nmglut, v_npglut, file = "18line_pseudobulk_limma_post_normfactor_voom_dfs.RData")

# output results ####
fit_all <- list(fit_GABA_contr, fit_nmglut_contr, fit_npglut_contr)
setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm_9in18_res/")
cnames <- c("1v0 up", "1v0 down", "6v0 up", "6v0 down")
types <- c("GABA", "nmglut", "npglut")
for (i in 1:3){
  res_1v0 <- topTable(fit_all[[i]], coef = "onevszero", p.value = 0.05, sort.by = "P", number = Inf)
  res_6v0 <- topTable(fit_all[[i]], coef = "sixvszero", p.value = 0.05, sort.by = "P", number = Inf)
  res_6v1 <- topTable(fit_all[[i]], coef = "sixvsone", p.value = 0.05, sort.by = "P", number = Inf)
  res_1v0$gene <- rownames(res_1v0)
  res_6v0$gene <- rownames(res_6v0)
  res_6v1$gene <- rownames(res_6v1)
  res_all <- list(res_1v0, res_6v0, res_6v1)
  filenameElements <- c("res_1v0", "res_6v0", "res_6v1")
  for (j in 1:length(res_all)) {
    filename <- paste0("./upregulated_significant/", types[i], "_", 
                       filenameElements[j], "_upregulated_DEGs.csv")
    print(filename)
    write.table(res_all[[j]][res_all[[j]]$logFC > 0, ], file = filename, 
                quote = F, sep = ",", row.names = F, col.names = T)
    
    filename <- paste0("./downregulated_significant/", types[i], "_", 
                       filenameElements[j], "_downregulated_DEGs.csv")
    print(filename)
    write.table(res_all[[j]][res_all[[j]]$logFC < 0, ], file = filename, 
                quote = F, sep = ",", row.names = F, col.names = T)
  }
  
  res_1v0 <- topTable(fit_all[[i]], coef = "onevszero", p.value = Inf, sort.by = "P", number = Inf)
  res_6v0 <- topTable(fit_all[[i]], coef = "sixvszero", p.value = Inf, sort.by = "P", number = Inf)
  res_6v1 <- topTable(fit_all[[i]], coef = "sixvsone", p.value = Inf, sort.by = "P", number = Inf)
  res_1v0$gene <- rownames(res_1v0)
  res_6v0$gene <- rownames(res_6v0)
  res_6v1$gene <- rownames(res_6v1)
  res_all <- list(res_1v0, res_6v0, res_6v1)
  filenameElements <- c("res_1v0", "res_6v0", "res_6v1")
  for (j in 1:length(res_all)) {
    filename <- paste0("./unfiltered_by_padj/upregulated/", types[i], "_", 
                       filenameElements[j], "_upregulated_DEGs.csv")
    print(filename)
    write.table(res_all[[j]][res_all[[j]]$logFC > 0, ], file = filename, 
                quote = F, sep = ",", row.names = F, col.names = T)
    
    filename <- paste0("./unfiltered_by_padj/downregulated/", types[i], "_", 
                       filenameElements[j], "_downregulated_DEGs.csv")
    print(filename)
    write.table(res_all[[j]][res_all[[j]]$logFC < 0, ], file = filename, 
                quote = F, sep = ",", row.names = F, col.names = T)
    
    filename <- paste0("./unfiltered_by_padj/", types[i], "_", 
                       filenameElements[j], "_all_DEGs.csv")
    print(filename)
    write.table(res_all[[j]], file = filename, 
                quote = F, sep = ",", row.names = F, col.names = T)
  }
} 

# summary table ####
deg_counts <- array(dim = c(3, 4), dimnames = list(types, cnames)) # 3 cell types are rows, 1/6 up/down are cols
for (i in 1:3){                                                               
  res_1v0 <- topTable(fit_all[[i]], coef = "onevszero", p.value = 0.05, 
                      sort.by = "P", number = Inf, adjust.method = "bonferroni")
  res_6v0 <- topTable(fit_all[[i]], coef = "sixvszero", p.value = 0.05, 
                      sort.by = "P", number = Inf, adjust.method = "bonferroni")
  deg_counts[i, 1] <- sum(res_1v0$logFC > 0)
  deg_counts[i, 2] <- sum(res_1v0$logFC < 0)     
  deg_counts[i, 3] <- sum(res_6v0$logFC > 0)
  deg_counts[i, 4] <- sum(res_6v0$logFC < 0)
}
write.table(deg_counts, file = "./deg_summary_combat.csv", quote = F, sep = ",", 
            row.names = T, col.names = T)



gene_ranges <- EnsDb.Hsapiens.v86
ensembldb::seqlevelsStyle(gene_ranges) <- "UCSC"

annotation_ranges <-
  GetGRangesFromEnsDb(ensdb = gene_ranges,
                      verbose = T, 
                      standard.chromosomes = T)
