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


# Limma make DGE and filter ####
design_GABA <- model.matrix(~0 + time + group + age + sex + disease + GABA_fraction, data = covar_18line) #26353 genes
y_GABA <- createDGE(GABA_mat_adj, 10) # plug in adjusted or raw matrix
keep <- filterByCpm(y_GABA, 1, 12)
sum(keep) #16977
y_GABA <- calcNormFactors(y_GABA[keep,, keep.lib.sizes= F])

design_nmglut <- model.matrix(~0 + time + group + age + sex + disease + nmglut_fraction, data = covar_18line) #26353 genes
y_nmglut <- createDGE(nmglut_mat_adj, 10)
keep <- filterByCpm(y_nmglut, 1, 12)
sum(keep) #16977
y_nmglut <- calcNormFactors(y_nmglut[keep,, keep.lib.sizes= F])

design_npglut <- model.matrix(~0 + time + group + age + sex + disease + npglut_fraction, data = covar_18line) #29400 genes
y_npglut <- createDGE(npglut_mat_adj, 10)
keep <- filterByCpm(y_npglut, 1, 12)
sum(keep) #16977
y_npglut <- calcNormFactors(y_npglut[keep,, keep.lib.sizes= F])

# edgeR ####
edgrOps <- function(y, design) {
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf1 <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0, 0, 0, 0, 0))
  qlf6 <- glmQLFTest(fit, contrast = c(-1, 0, 1, 0, 0, 0, 0, 0))
  res1 <- topTags(qlf1, n = Inf, p.value = Inf)
  res6 <- topTags(qlf6, n = Inf, p.value = Inf)
  return(list(y, fit, res1, res6))
}

res_GABA <- edgrOps(y_GABA, design_GABA)
res_nmglut <- edgrOps(y_nmglut, design_nmglut)
res_npglut <- edgrOps(y_npglut, design_npglut)

plotBCV(res_GABA[[1]], main = "Final mean-variance trend", col.common = "transparent")
plotBCV(res_nmglut[[1]], main = "Final mean-variance trend", col.common = "transparent")
plotBCV(res_npglut[[1]], main = "Final mean-variance trend", col.common = "transparent")

#plotSA(res_GABA[[2]])
#plotMD(res_GABA[[2]], main = "GABA 0v6hr")

# output results ####
fit_all <- list(res_GABA[[3]], res_nmglut[[3]], res_npglut[[3]],
                res_GABA[[4]], res_nmglut[[4]], res_npglut[[4]])
setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/pseudobulk_DE/res_use_combatseq_mat/edgeR_res")
cnames <- c("1v0 up", "1v0 down", "6v0 up", "6v0 down")
types <- c("GABA", "nmglut", "npglut")
for (i in 1:3){
  res_1v0 <- fit_all[[i]]$table
  res_6v0 <- fit_all[[i + 3]]$table
  res_all <- list(res_1v0, res_6v0)
  filenameElements <- c("res_1v0", "res_6v0")
  for (j in 1:length(res_all)) {
    filename <- paste0("./upregulated_significant/", types[i], "_", 
                       filenameElements[j], "_upregulated_DEGs.csv")
    print(filename)
    write.table(res_all[[j]][res_all[[j]]$logFC > 0 & res_all[[j]]$FDR < 0.05, ], 
                file = filename, quote = F, sep = ",", row.names = F, col.names = T)
    
    filename <- paste0("./downregulated_significant/", types[i], "_", 
                       filenameElements[j], "_downregulated_DEGs.csv")
    print(filename)
    write.table(res_all[[j]][res_all[[j]]$logFC < 0 & res_all[[j]]$FDR < 0.05, ], 
                file = filename, quote = F, sep = ",", row.names = F, col.names = T)
  }
  
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
  res_1v0 <- fit_all[[i]]
  res_6v0 <- fit_all[[i + 3]]
  deg_counts[i, 1] <- sum(res_1v0$table$logFC > 0 & res_1v0$table$FDR < 0.05)
  deg_counts[i, 2] <- sum(res_1v0$table$logFC < 0 & res_1v0$table$FDR < 0.05)     
  deg_counts[i, 3] <- sum(res_6v0$table$logFC > 0 & res_6v0$table$FDR < 0.05)
  deg_counts[i, 4] <- sum(res_6v0$table$logFC < 0 & res_6v0$table$FDR < 0.05)
}
write.table(deg_counts, file = "./deg_summary_combat.csv", quote = F, sep = ",", 
            row.names = T, col.names = T)

