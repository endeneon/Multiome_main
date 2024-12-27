# Chuxuan Li 07/19/2022
# DE analysis using Limma blocking out all covariates (not cell line, but include
#group) same as those in the data transfer to Lifan in June 2022
# try different filters for genes: 1) total exprs > 20, 50, 100, 
#2) >1cpm in at least 0.5*18=9 samples

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
load("./pseudobulk_DE/res_use_combatseq_mat/by_line_pseudobulk_mat_raw.RData")
load("./pseudobulk_DE/res_use_combatseq_mat/by_line_pseudobulk_mat_adjusted.RData")
celllines <- str_remove(colnames(GABA_mat), "_[0|1|6]hr$")

covar_18line <- covar_table[covar_table$cell_line %in% lines, ]
covar_18line <- covar_18line[order(covar_18line$cell_line), ]
rm(covar_table)

createDGE <- function(count_matrix, isexpr_val){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = celllines)
  A <- rowSums(y$counts)
  isexpr <- A > 10
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[isexpr & hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}

# Limma on uncorrected matrices with trend ####

limmaOps <- function(y, design) {
  y <- calcNormFactors(y)
  # logCPM
  logCPM <- cpm(y, log = T, normalized.lib.sizes = T)
  # fit
  contr.matrix <- makeContrasts(
    onevszero = time1hr-time0hr, 
    sixvszero = time6hr-time0hr, 
    sixvsone = time6hr-time1hr, levels = colnames(design))
  fit <- lmFit(logCPM, design)
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit <- eBayes(fit, trend = T)
  summary(decideTests(fit))
  return(fit)  
}


design_GABA <- model.matrix(~0 + time + group + age + sex + disease + GABA_fraction, data = covar_18line) #26353 genes
y_GABA <- createDGE(GABA_mat, 10)
plotSA(y_GABA, main="Mean-variance trend")
fit_GABA <- limmaOps(y_GABA, design_GABA)
plotSA(fit_GABA, main="Mean-variance trend")

design_nmglut <- model.matrix(~0 + time + group + age + sex + disease + nmglut_fraction, data = covar_18line) #26353 genes
fit_nmglut <- limmaOps(nmglut_mat, 10, design_nmglut)

design_npglut <- model.matrix(~0 + time + group + age + sex + disease + npglut_fraction, data = covar_18line) #29400 genes
fit_npglut <- limmaOps(npglut_mat, 10, design_npglut)

fit_all <- list(fit_GABA, fit_nmglut, fit_npglut)

AT1G60030 <- y$E["AT1G60030",]
plot(AT1G60030 ~ pH, ylim = c(0, 3.5))
intercept <- coef(fit)["AT1G60030", "(Intercept)"]
slope <- coef(fit)["AT1G60030", "pH"]
abline(a = intercept, b = slope)

setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/pseudobulk_DE/res_use_combatseq_mat/limma_mat_raw_res")

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

# summary table
cnames <- c("1v0 up", "1v0 down", "6v0 up", "6v0 down")
types <- c("GABA", "nmglut", "npglut")
deg_counts <- array(dim = c(3, 4), dimnames = list(types, cnames)) # 3 cell types are rows, 1/6 up/down are cols
for (i in 1:3){                                                               
  res_1v0 <- topTable(fit_all[[i]], coef = "onevszero", p.value = 0.05, sort.by = "P", number = Inf)
  res_6v0 <- topTable(fit_all[[i]], coef = "sixvszero", p.value = 0.05, sort.by = "P", number = Inf)
  deg_counts[i, 1] <- sum(res_1v0$logFC > 0)
  deg_counts[i, 2] <- sum(res_1v0$logFC < 0)     
  deg_counts[i, 3] <- sum(res_6v0$logFC > 0)
  deg_counts[i, 4] <- sum(res_6v0$logFC < 0)
}
write.table(deg_counts, file = "./deg_summary_limma_on_mat_raw.csv", quote = F, sep = ",", 
            row.names = T, col.names = T)


# Limma on combat corrected matrices ####
design_GABA <- model.matrix(~0 + time + cell_line + GABA_fraction, data = covar_18line) #27597 genes
fit_GABA <- limmaOps(GABA_mat_adj, 10, design_GABA)

design_nmglut <- model.matrix(~0 + time + cell_line + nmglut_fraction, data = covar_18line) #26353 genes
fit_nmglut <- limmaOps(nmglut_mat_adj, 10, design_nmglut)

design_npglut <- model.matrix(~0 + time + cell_line + npglut_fraction, data = covar_18line) #29400 genes
fit_npglut <- limmaOps(npglut_mat_adj, 10, design_npglut)

fit_all <- list(fit_GABA, fit_nmglut, fit_npglut)
setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/pseudobulk_DE/res_use_combatseq_mat/limma_mat_adj_res")

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

# summary table
cnames <- c("1v0 up", "1v0 down", "6v0 up", "6v0 down")
types <- c("GABA", "nmglut", "npglut")
deg_counts <- array(dim = c(3, 4), dimnames = list(types, cnames)) # 3 cell types are rows, 1/6 up/down are cols
for (i in 1:3){                                                               
  res_1v0 <- topTable(fit_all[[i]], coef = "onevszero", p.value = 0.05, sort.by = "P", number = Inf)
  res_6v0 <- topTable(fit_all[[i]], coef = "sixvszero", p.value = 0.05, sort.by = "P", number = Inf)
  deg_counts[i, 1] <- sum(res_1v0$logFC > 0)
  deg_counts[i, 2] <- sum(res_1v0$logFC < 0)     
  deg_counts[i, 3] <- sum(res_6v0$logFC > 0)
  deg_counts[i, 4] <- sum(res_6v0$logFC < 0)
}

write.table(deg_counts, file = "./deg_summary_limma_on_mat_adj.csv", quote = F, 
            sep = ",", row.names = T, col.names = T)

