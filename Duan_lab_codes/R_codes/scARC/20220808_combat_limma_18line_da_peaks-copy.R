# Chuxuan Li 08/08/2022
# Use combat + Limma to compute pseudobulk differentially accessible peaks
#for 18-line ATACseq data

# init ####
library(limma)
library(edgeR)
library(sva)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(readr)
library(Signac)

# setwd("~/NVME/scARC_Duan_018/Duan_Project_018-ATAC/")

load("~/Data/FASTQ/Duan_Project_024/hybrid_output/R_DE_peaks_pseudobulk/frag_count_by_cell_type_specific_peaks_22Jul2022.RData")

load("../../covariates_pooled_together_for_all_linextime_5+18+20.RData")

# define covariates from covar table
covar_18line <- covar_table[covar_table$cell_line %in% lines, ]
covar_18line <- covar_18line[order(covar_18line$cell_line), ]
rm(covar_table)
GABA_fraction <- covar_18line$GABA_fraction
nmglut_fraction <- covar_18line$nmglut_fraction
npglut_fraction <- covar_18line$npglut_fraction

# make combat dataframes ####
GABA <- Seurat_object_list[[2]]
nmglut <- Seurat_object_list[[3]]
npglut <- Seurat_object_list[[4]]

cellines <- sort(unique(nmglut$cell.line.ident))
times <- sort(unique(nmglut$time.ident))
nmglut$timexline.ident <- NA
for (l in cellines) {
  for (t in times) {
    id <- paste(l, t, sep = "-")
    print(id)
    nmglut$timexline.ident[nmglut$cell.line.ident == l &
                             nmglut$time.ident == t] <- id
  }
}
unique(nmglut$timexline.ident)
lts <- sort(unique(nmglut$timexline.ident))

# make raw count matrix ####
makeMat4CombatseqLine <- function(typeobj, linetimes) {
  rownames <- rep_len(NA, (length(linetimes)))
  for (i in 1:length(linetimes)) {
    print(linetimes[i])
    rownames[i] <- linetimes[i]
    obj <- subset(typeobj, timexline.ident == linetimes[i])
    count_matrix <- obj 
    if (i == 1) {
      mat <- as.array(rowSums(obj@assays$ATAC@counts))
    } else {
      to_bind <- as.array(rowSums(obj@assays$ATAC@counts))
      mat <- rbind(mat, to_bind)
    }
  }
  rownames(mat) <- rownames
  return(mat)
}

GABA_mat <- t(makeMat4CombatseqLine(GABA, lts))
nmglut_mat <- t(makeMat4CombatseqLine(nmglut, lts))
npglut_mat <- t(makeMat4CombatseqLine(npglut, lts))

# combat-seq ####
batch = rep(c("12", "11", "12", "12", "8", "8", "9", "9", "10", "8", "8", "9", 
              "10", "10", "11", "11", "11", "10"), each = 3)
batch
length(batch) #54
ncol(GABA_mat) #54

GABA_mat_adj <- ComBat_seq(GABA_mat, batch = batch, 
                           group = NULL)
nmglut_mat_adj <- ComBat_seq(nmglut_mat, batch = batch,
                             group = NULL)
npglut_mat_adj <- ComBat_seq(npglut_mat, batch = batch,
                             group = NULL)
save(GABA_mat_adj, nmglut_mat_adj, npglut_mat_adj, file = "post-combat_count_matrices.RData")


celllines <- str_extract(lts, "^CD_[0-9][0-9]")


# auxillary functions ####
createDGE <- function(count_matrix){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = celllines)
  A <- rowSums(y$counts)
  #isexpr <- A > 10
  hasant <- rowSums(is.na(y$genes)) == 0
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
y_GABA <- createDGE(GABA_mat_adj)
ind.keep <- filterByCpm(y_GABA, 1, 9)
sum(ind.keep) #311555 -> 172759
v_GABA <- cnfV(y_GABA[ind.keep, ], design_GABA)
fit_GABA <- lmFit(v_GABA, design_GABA)
fit_GABA_contr <- contrastFit(fit_GABA, design_GABA)
plotSA(fit_GABA_contr, main = "Final model: Mean-variance trend") # same, just has a line

design_nmglut <- model.matrix(~0 + time + group + age + sex + disease + nmglut_fraction, data = covar_18line) #26353 genes
y_nmglut <- createDGE(nmglut_mat_adj)
ind.keep <- filterByCpm(y_nmglut, 1, 9)
sum(ind.keep) #325571 -> 196083
v_nmglut <- cnfV(y_nmglut[ind.keep, ], design_nmglut)
fit_nmglut <- lmFit(v_nmglut, design_nmglut)
fit_nmglut_contr <- contrastFit(fit_nmglut, design_nmglut)
plotSA(fit_nmglut_contr, main = "Final model: Mean-variance trend") 

design_npglut <- model.matrix(~0 + time + group + age + sex + disease + npglut_fraction, data = covar_18line) #29400 genes
y_npglut <- createDGE(npglut_mat_adj)
ind.keep <- filterByCpm(y_npglut, 1, 9)
sum(ind.keep) #340779 -> 207450
v_npglut <- cnfV(y_npglut[ind.keep, ], design_npglut)
fit_npglut <- lmFit(v_npglut, design_npglut)
fit_npglut_contr <- contrastFit(fit_npglut, design_npglut)
plotSA(fit_npglut_contr, main = "Final model: Mean-variance trend") 


# output results ####
setwd("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/da_peaks_limma")

fit_all <- list(fit_GABA_contr, fit_nmglut_contr, fit_npglut_contr)
cnames <- c("1v0 up", "1v0 down", "6v0 up", "6v0 down")
types <- c("GABA", "nmglut", "npglut")
for (i in 1:3){
  res_1v0 <- topTable(fit_all[[i]], coef = "onevszero", p.value = 0.05, sort.by = "P", number = Inf)
  res_6v0 <- topTable(fit_all[[i]], coef = "sixvszero", p.value = 0.05, sort.by = "P", number = Inf)
  res_6v1 <- topTable(fit_all[[i]], coef = "sixvsone", p.value = 0.05, sort.by = "P", number = Inf)
  res_1v0$peak <- rownames(res_1v0)
  res_6v0$peak <- rownames(res_6v0)
  res_6v1$peak <- rownames(res_6v1)
  res_all <- list(res_1v0, res_6v0, res_6v1)
  filenameElements <- c("res_1v0", "res_6v0", "res_6v1")
  for (j in 1:length(res_all)) {
    filename <- paste0("./upregulated_significant/", types[i], "_", 
                       filenameElements[j], "_upregulated_peaks.csv")
    print(filename)
    write.table(res_all[[j]][res_all[[j]]$logFC > 0, ], file = filename, 
                quote = F, sep = ",", row.names = F, col.names = T)
    
    filename <- paste0("./downregulated_significant/", types[i], "_", 
                       filenameElements[j], "_downregulated_peaks.csv")
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
                       filenameElements[j], "_upregulated_peaks.csv")
    print(filename)
    write.table(res_all[[j]][res_all[[j]]$logFC > 0, ], file = filename, 
                quote = F, sep = ",", row.names = F, col.names = T)
    
    filename <- paste0("./unfiltered_by_padj/downregulated/", types[i], "_", 
                       filenameElements[j], "_downregulated_peaks.csv")
    print(filename)
    write.table(res_all[[j]][res_all[[j]]$logFC < 0, ], file = filename, 
                quote = F, sep = ",", row.names = F, col.names = T)
    
    filename <- paste0("./unfiltered_by_padj/", types[i], "_", 
                       filenameElements[j], "_all_peaks.csv")
    print(filename)
    write.table(res_all[[j]], file = filename, 
                quote = F, sep = ",", row.names = F, col.names = T)
  }
} 

# summary table ####
da_counts <- array(dim = c(3, 4), dimnames = list(types, cnames)) # 3 cell types are rows, 1/6 up/down are cols
for (i in 1:3){                                                               
  res_1v0 <- topTable(fit_all[[i]], coef = "onevszero", p.value = 0.05, 
                      sort.by = "P", number = Inf)
  res_6v0 <- topTable(fit_all[[i]], coef = "sixvszero", p.value = 0.05, 
                      sort.by = "P", number = Inf)
  da_counts[i, 1] <- sum(res_1v0$logFC > 0)
  da_counts[i, 2] <- sum(res_1v0$logFC < 0)     
  da_counts[i, 3] <- sum(res_6v0$logFC > 0)
  da_counts[i, 4] <- sum(res_6v0$logFC < 0)
}
write.table(da_counts, file = "./da_peaks_limma/da_peaks_count_summary.csv", quote = F, sep = ",", 
            row.names = T, col.names = T)

# save bed files ####
da_peaks_by_time_list_1v0 <- vector("list", 3)
da_peaks_by_time_list_6v0 <- vector("list", 3)
for (i in 1:length(da_peaks_by_time_list_1v0)) {
  da_peaks_by_time_list_1v0[[i]] <- topTable(fit_all[[i]], coef = "onevszero", 
                                             p.value = Inf, sort.by = "P", 
                                             number = Inf)
  da_peaks_by_time_list_6v0[[i]] <- topTable(fit_all[[i]], coef = "sixvszero", 
                                             p.value = Inf, sort.by = "P", 
                                             number = Inf)
}
saveBed <- function(df, file_name){
  df_bed <- data.frame(chr = str_split(string = df$peaks, pattern = "-", simplify = T)[, 1],
                       start = gsub(pattern = "-",
                                    replacement = "",
                                    x = str_extract_all(string = df$peaks,
                                                        pattern = "-[0-9]+-"
                                    )),
                       end = gsub(pattern = "-",
                                  replacement = "",
                                  x = str_extract_all(string = df$peaks,
                                                      pattern = "-[0-9]+$"
                                  )),
                       id = paste0(df$avg_log2FC, '^',
                                   df$p_val, '^',
                                   df$p_val_adj), # log2FC^pval^pval_adj, for annotation
                       p_val = df$p_val,
                       strand = rep("+", length(df$p_val)))
  write.table(df_bed,
              file = file_name,
              quote = F,
              sep = "\t",
              row.names = F,
              col.names = F)
}

lists <- list(da_peaks_by_time_list_1v0, da_peaks_by_time_list_6v0)
times <- c("1", "6")
types
#setwd("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/da_peaks_limma")
for (j in 1:length(lists)) {
  time <- times[j]
  list <- lists[[j]]
  for (i in 1:length(list)){
    type <- types[i]
    ident <- paste(type, time, sep = "_")
    print(ident)
    #list[[i]] <- list[[i]][list[[i]]$adj.P.Val < 0.05, ]
    uplist <- list[[i]][list[[i]]$logFC > 0, ]
    dolist <- list[[i]][list[[i]]$logFC < 0, ]
    # df <- data.frame(peaks = rownames(uplist),
    #                  p_val = uplist$P.Value,
    #                  avg_log2FC = uplist$logFC,
    #                  p_val_adj = uplist$adj.P.Val)
    # file_name <- paste0("./unfiltered_bed/",
    #                     "./bed/",
    #                     ident, "v0hr_upregulated.bed")
    # saveBed(df, file_name)
    # df <- data.frame(peaks = rownames(dolist),
    #                  p_val = dolist$P.Value,
    #                  avg_log2FC = dolist$logFC,
    #                  p_val_adj = dolist$adj.P.Val)
    # file_name <- paste0("./unfiltered_bed/",
    #   "./bed/",
    #   ident, "v0hr_downregulated.bed")
    # saveBed(df, file_name)
    df <- data.frame(peaks = rownames(list[[i]]),
                     p_val = list[[i]]$P.Value,
                     avg_log2FC = list[[i]]$logFC,
                     p_val_adj = list[[i]]$adj.P.Val)
    file_name <- paste0("./unfiltered_bed/",
                        #"./bed/",
                        ident, "v0hr_full.bed")
    saveBed(df, file_name)
  }
}
