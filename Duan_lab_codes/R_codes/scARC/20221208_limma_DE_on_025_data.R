# Chuxuan Li 12/8/2022
# pseudobulk DE with limma on 025 RNAseq data

# init ####
library(Seurat)
library(limma)
library(edgeR)
library(sva)

library(readr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

# make count df ####
load("~/NVME/scARC_Duan_018/Duan_project_025_RNA/025_raw_pseudobulk_matrix_for_varpartition.RData")

GABA_mat <- export_df[, str_detect(colnames(export_df), "GABA")]
nmglut_mat <- export_df[, str_detect(colnames(export_df), "nmglut")]
npglut_mat <- export_df[, str_detect(colnames(export_df), "npglut")]
GABA_mat <- GABA_mat[, sort(colnames(GABA_mat))]
nmglut_mat <- nmglut_mat[, sort(colnames(nmglut_mat))]
npglut_mat <- npglut_mat[, sort(colnames(npglut_mat))]

# combat-seq correct for batch effect ####
lines <- str_extract(colnames(GABA_mat), "CD_[0-9][0-9]")
load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/integrated_labeled.RData")
batch <- c()
for (l in lines) {
  print(l)
  batch <- c(batch, 
             unique(str_split(integrated_labeled$orig.ident[integrated_labeled$cell.line.ident == l],
                              pattern = "-", n = 2, simplify = T)[1]))
}
batch <- as.factor(batch)
group <- str_extract(colnames(GABA_mat), "[0|1|6]hr")
GABA_mat_adj <- ComBat_seq(as.matrix(GABA_mat), batch = batch, group = group)
nmglut_mat_adj <- ComBat_seq(as.matrix(nmglut_mat), batch = batch, group = group)
npglut_mat_adj <- ComBat_seq(as.matrix(npglut_mat), batch = batch, group = group)
GABA_mat_adj <- GABA_mat_adj[, sort(colnames(GABA_mat_adj))]
nmglut_mat_adj <- nmglut_mat_adj[, sort(colnames(nmglut_mat_adj))]
npglut_mat_adj <- npglut_mat_adj[, sort(colnames(npglut_mat_adj))]

setwd("./Analysis_part2_GRCh38/")
save(GABA_mat, nmglut_mat, npglut_mat, file = "three_celltype_pseudobulk_mat_raw.RData")
save(GABA_mat_adj, nmglut_mat_adj, npglut_mat_adj, 
     file = "three_celltype_pseudobulk_mat_corrected.RData")

# make covariate table ####
MGS_table <- read_excel("/nvmefs/scARC_Duan_018/MGS_iPSC lines_60samples_scRNA_ATAC-seq_bulkATAC-seq status_AK_HZ.xlsx")
MGS_table$`new iPSC line ID` <- str_replace(MGS_table$`new iPSC line ID`,
                                            "CD00000",
                                            "CD_")
covar_table <- MGS_table[, c(1, 3, 4, 5, 9)]
covar_table$`sex (1=M, 2=F)` <- str_replace(str_replace(covar_table$`sex (1=M, 2=F)`,
                                            "1", "M"), "2", "F")
colnames(covar_table) <- c("cell_line", "aff", "sex", "age", "batch")
covar_table <- covar_table[!is.na(covar_table$batch), ]
covar_table$batch <- str_remove(covar_table$batch, "[0-9]+ & ")
covar_table$batch <- str_remove(covar_table$batch, "[0-9]+ \\(need redo\\)  & ")
# for (i in unique(covar_table$batch)) {
#   cat(i, ": ", sum(covar_table$batch == i), "\n")
# }
integrated_labeled$linextime.ident <- NA
times <- sort(unique(integrated_labeled$time.ident))
lines <- sort(unique(integrated_labeled$cell.line.ident))
for (t in times) {
  for (l in lines) {
    integrated_labeled$linextime.ident[integrated_labeled$time.ident == t &
                                         integrated_labeled$cell.line.ident == l] <-
      paste(l, t, sep = "_")
  }
}
lts <- unique(integrated_labeled$linextime.ident)
med_df <- data.frame(cell_line = rep_len("", length(lts)),
                     time = rep_len("", length(lts)),
                     GABA_counts = rep_len("", length(lts)),
                     GABA_fraction = rep_len("", length(lts)),
                     nmglut_counts = rep_len("", length(lts)),
                     nmglut_fraction = rep_len("", length(lts)),
                     npglut_counts = rep_len("", length(lts)),
                     npglut_fraction = rep_len("", length(lts)),
                     total_counts = rep_len("", length(lts)))
unique(integrated_labeled$cell.type)
integrated_labeled <- subset(integrated_labeled, cell.type != "unidentified")
for (i in 1:length(lts)) {
  obj <- subset(integrated_labeled, linextime.ident == lts[[i]])
  line <- unique(obj$cell.line.ident)
  time <- unique(obj$time.ident)
  types <- unique(obj$cell.type)
  print(line)
  print(time)
  print(types)
  med_df$cell_line[i] <- line
  med_df$time[i] <- time
  med_df$GABA_counts[i] <- sum(obj$cell.type == "GABA")
  med_df$nmglut_counts[i] <- sum(obj$cell.type == "nmglut")
  med_df$npglut_counts[i] <- sum(obj$cell.type == "npglut")
  med_df$total_counts[i] <- ncol(obj)
  med_df$GABA_fraction[i] <- sum(obj$cell.type == "GABA") / ncol(obj)
  med_df$nmglut_fraction[i] <- sum(obj$cell.type == "nmglut") / ncol(obj)
  med_df$npglut_fraction[i] <- sum(obj$cell.type == "npglut") / ncol(obj)
}

# combine into covariate table
covar_df <- med_df
covar_df$aff <- rep_len("", nrow(covar_df))
covar_df$age <- rep_len("", nrow(covar_df))
covar_df$sex <- rep_len("", nrow(covar_df))
covar_df$batch <- rep_len("", nrow(covar_df))
for (i in 1:length(covar_df$cell_line)) {
  line <- covar_df$cell_line[i]
  print(line)
  covar_df$aff[i] <- covar_table$aff[covar_table$cell_line == line]
  covar_df$age[i] <- covar_table$age[covar_table$cell_line == line]
  covar_df$sex[i] <- covar_table$sex[covar_table$cell_line == line]
  covar_df$batch[i] <- covar_table$batch[covar_table$cell_line == line]
}
covar_df <- covar_df[order(covar_df$cell_line), ]
covar_df <- covar_df[order(covar_df$time), ]
covar_df$GABA_fraction <- as.numeric(covar_df$GABA_fraction)
covar_df$nmglut_fraction <- as.numeric(covar_df$nmglut_fraction)
covar_df$npglut_fraction <- as.numeric(covar_df$npglut_fraction)
covar_df$age <- as.numeric(covar_df$age)
save(covar_df, file = "025_covariates_full_df.RData")

# limma ####
createDGE <- function(count_matrix){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix))
  A <- rowSums(y$counts)
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}
filterByCpm <- function(y, cutoff, nsample) {
  cpm <- cpm(y)
  colInd0 <- str_detect(colnames(y), "0hr")
  colInd1 <- str_detect(colnames(y), "1hr")
  colInd6 <- str_detect(colnames(y), "6hr")
  
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
design_GABA <- model.matrix(~0 + time + batch + age + sex + aff + GABA_fraction, 
                            data = covar_df)
y_GABA <- createDGE(GABA_mat_adj) # plug in adjusted or raw matrix
ind.keep <- filterByCpm(y_GABA, 1, 10)
sum(ind.keep) #18865
v_GABA <- cnfV(y_GABA[ind.keep, ], design_GABA)
fit_GABA <- lmFit(v_GABA, design_GABA)
fit_GABA_contr <- contrastFit(fit_GABA, design_GABA)
plotSA(fit_GABA_contr, main = "Final model: Mean-variance trend") # same, just has a line

design_nmglut <- model.matrix(~0 + time + batch + age + sex + aff + nmglut_fraction, 
                              data = covar_df)
y_nmglut <- createDGE(nmglut_mat_adj)
ind.keep <- filterByCpm(y_nmglut, 1, 10)
sum(ind.keep) #18869
v_nmglut <- cnfV(y_nmglut[ind.keep, ], design_nmglut)
fit_nmglut <- lmFit(v_nmglut, design_nmglut)
fit_nmglut_contr <- contrastFit(fit_nmglut, design_nmglut)
plotSA(fit_nmglut_contr, main = "Final model: Mean-variance trend") 

design_npglut <- model.matrix(~0 + time + group + age + sex + aff + npglut_fraction, 
                              data = covar_df) 
y_npglut <- createDGE(npglut_mat_adj)
ind.keep <- filterByCpm(y_npglut, 1, 10)
sum(ind.keep) #17996
v_npglut <- cnfV(y_npglut[ind.keep, ], design_npglut)
fit_npglut <- lmFit(v_npglut, design_npglut)
fit_npglut_contr <- contrastFit(fit_npglut, design_npglut)
plotSA(fit_npglut_contr, main = "Final model: Mean-variance trend") 

fit_all <- list(fit_GABA_contr, fit_nmglut_contr, fit_npglut_contr)
setwd("./limma_DE_results/")
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

# summary table
cnames <- c("1v0 up", "1v0 down", "6v0 up", "6v0 down", "total")
deg_counts <- array(dim = c(3, 5), dimnames = list(types, cnames)) # 3 cell types are rows, 1/6 up/down are cols
for (i in 1:3){                                                               
  res_1v0 <- topTable(fit_all[[i]], coef = "onevszero", p.value = Inf, sort.by = "P", number = Inf)
  res_6v0 <- topTable(fit_all[[i]], coef = "sixvszero", p.value = Inf, sort.by = "P", number = Inf)
  deg_counts[i, 1] <- sum(res_1v0$logFC > 0 & res_1v0$adj.P.Val < 0.05)
  deg_counts[i, 2] <- sum(res_1v0$logFC < 0 & res_1v0$adj.P.Val < 0.05)     
  deg_counts[i, 3] <- sum(res_6v0$logFC > 0 & res_6v0$adj.P.Val < 0.05)
  deg_counts[i, 4] <- sum(res_6v0$logFC < 0 & res_6v0$adj.P.Val < 0.05)
  deg_counts[i, 5] <- nrow(res_1v0)
}
write.table(deg_counts, file = "./deg_summary_limma.csv", quote = F, sep = ",", 
            row.names = T, col.names = T)

