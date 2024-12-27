# Chuxuan Li 07/11/2022
# DE analysis using Limma blocking out cell line that covers other covariates
#(age, sex, disease, group)

# init ####
library(Seurat)
library(limma)
library(edgeR)
library(sva)

library(readr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

load("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")
RNA_18line <- integrated_labeled
rm(integrated_labeled)

unique(RNA_18line$cell.type)

GABA <- subset(RNA_18line, cell.type %in% c("GABA", "SST_pos_GABA", "SEMA3E_pos_GABA"))
npglut <- subset(RNA_18line, cell.type == "NEFM_pos_glut")
nmglut <- subset(RNA_18line, cell.type == "NEFM_neg_glut")

# make count df ####
makePseudobulkMat <- function(typeobj, linetimes) {
  rownames <- rep_len(NA, (length(linetimes)))
  for (i in 1:length(linetimes)) {
    print(linetimes[i])
    rownames[i] <- linetimes[i]
    obj <- subset(typeobj, timexline.ident == linetimes[i])
    count_matrix <- obj 
    if (i == 1) {
      mat <- as.array(rowSums(obj@assays$RNA@counts))
    } else {
      to_bind <- as.array(rowSums(obj@assays$RNA@counts))
      mat <- rbind(mat, to_bind)
    }
  }
  rownames(mat) <- rownames
  return(mat)
}

lts <- sort(unique(RNA_18line$timexline.ident))
lines <- sort(unique(RNA_18line$cell.line.ident))
rm(RNA_18line)
GABA_mat <- t(makePseudobulkMat(GABA, lts))
nmglut_mat <- t(makePseudobulkMat(nmglut, lts))
npglut_mat <- t(makePseudobulkMat(npglut, lts))
GABA_mat <- GABA_mat[, sort(colnames(GABA_mat))]
nmglut_mat <- nmglut_mat[, sort(colnames(nmglut_mat))]
npglut_mat <- npglut_mat[, sort(colnames(npglut_mat))]

# combat-seq correct for batch effect - by group
batch <- str_split(string = lts, pattern = "_", n = 3, simplify = T)[,2]

GABA_mat_adj <- ComBat_seq(GABA_mat, batch = batch, 
                           #group = group
                           group = NULL)
nmglut_mat_adj <- ComBat_seq(nmglut_mat, batch = batch,
                             #group = group
                             group = NULL)
npglut_mat_adj <- ComBat_seq(npglut_mat, batch = batch,
                             #group = group
                             group = NULL)
GABA_mat_adj <- GABA_mat_adj[, sort(colnames(GABA_mat_adj))]
nmglut_mat_adj <- nmglut_mat_adj[, sort(colnames(nmglut_mat_adj))]
npglut_mat_adj <- npglut_mat_adj[, sort(colnames(npglut_mat_adj))]

setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/pseudobulk_DE/res_use_combatseq_mat")
save(GABA_mat, nmglut_mat, npglut_mat, file = "by_line_pseudobulk_mat_raw.RData")
save(GABA_mat_adj, nmglut_mat_adj, npglut_mat_adj, file = "by_line_pseudobulk_mat_adjusted.RData")

# Limma on uncorrected matrices ####
covar_add <- read_csv("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/covariates_with_matching_cell_line_id.csv")
covar_table <- read_csv("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/cell_type_proportions_by_cellline_time.csv")
covar_table$sex <- rep("", length(covar_table$cell_line))
covar_table$age <- rep("", length(covar_table$cell_line))
covar_table$disease <- rep("", length(covar_table$cell_line))
covar_table$group <- rep("", length(covar_table$cell_line))
for (i in 1:length(covar_table$cell_line)) {
  line <- covar_table$cell_line[i]
  covar_table$disease[i] <- covar_add$disease_status[covar_add$cell_line == line]
  covar_table$sex[i] <- covar_add$sex[covar_add$cell_line == line]
  covar_table$age[i] <- covar_add$age[covar_add$cell_line == line]
  covar_table$group[i] <- covar_add$batch[covar_add$cell_line == line]
}
save(covar_table, lines, file = "../../../covariates_pooled_together_for_all_linextime_5+18+20.RData")

covar_18line <- covar_table[covar_table$cell_line %in% lines, ]
covar_18line <- covar_18line[order(covar_18line$cell_line), ]

# try with GABA ####
# prefiltering and normalization

celllines <- str_remove(colnames(GABA_mat), "_[0|1|6]hr$")
# times <- str_split(colnames(GABA_mat), pattern = "_", n = 3, simplify = T)[,3]
# sexes <- rep_len("", length(lines))
# ages <- rep_len("", length(lines))
# GABA_prop <- rep_len("", length(lines))
# nmglut_prop <- rep_len("", length(lines))
# npglut_prop <- rep_len("", length(lines))
# dis <- rep_len("", length(lines))
# groups <- rep_len("", length(lines))
# for (i in 1:length(lines)){
#   sexes[i] <- covar_table$sex[covar_table$cell_line == lines[i] & covar_table$time == times[i]]
#   ages[i] <- covar_table$age[covar_table$cell_line == lines[i] & covar_table$time == times[i]]
#   GABA_prop[i] <- covar_table$GABA_fraction[covar_table$cell_line == lines[i] & covar_table$time == times[i]]
#   nmglut_prop[i] <- covar_table$nmglut_fraction[covar_table$cell_line == lines[i] & covar_table$time == times[i]]
#   npglut_prop[i] <- covar_table$npglut_fraction[covar_table$cell_line == lines[i] & covar_table$time == times[i]]
#   dis[i] <- covar_table$disease[covar_table$cell_line == lines[i] & covar_table$time == times[i]]
#   groups[i] <- covar_table$group[covar_table$cell_line == lines[i] & covar_table$time == times[i]]
# }
# cov_GABA <- as.data.frame(cbind(lines, times, ages, sexes, dis, GABA_prop))
y <- DGEList(counts = GABA_mat, genes = rownames(GABA_mat), group = celllines)
A <- rowSums(y$counts)
isexpr <- A > 10
hasant <- rowSums(is.na(y$genes)) == 0
y <- y[isexpr & hasant, , keep.lib.size = F]
dim(y) #27597 genes 
y <- calcNormFactors(y)
plotMDS(y, label = sprintf(covar_18line$GABA_fraction,
                           fmt = '%#.2f'))

# logCPM
logCPM <- cpm(y, log = T, normalized.lib.sizes = T)
col <- brewer.pal(ncol(y), "Paired")
boxplot(logCPM, las=2, col=col, main="") 

# fit
design <- model.matrix(~0 + time + cell_line + GABA_fraction, data = covar_18line)
contr.matrix <- makeContrasts(
  onevszero = time1hr-time0hr, 
  sixvszero = time6hr-time0hr, 
  sixvsone = time6hr-time1hr, levels = colnames(design))
fit <- lmFit(logCPM, design)
fit <- contrasts.fit(fit, contrasts = contr.matrix)
fit <- eBayes(fit, trend = T)
summary(decideTests(fit))
plotSA(fit, main="Mean-variance trend")

#topTable(fit, coef = "time6hr", p.value = 0.05, sort.by = "P", number = 20)


# function for limma operations ####
limmaOps <- function(count_matrix, isexpr_val, design){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = celllines)
  A <- rowSums(y$counts)
  isexpr <- A > 10
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[isexpr & hasant, , keep.lib.size = F]
  print(dim(y)) 
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

design_GABA <- model.matrix(~0 + time + cell_line + GABA_fraction, data = covar_18line) #26353 genes
fit_GABA <- limmaOps(GABA_mat, 10, design_GABA)

design_nmglut <- model.matrix(~0 + time + cell_line + nmglut_fraction, data = covar_18line) #26353 genes
fit_nmglut <- limmaOps(nmglut_mat, 10, design_nmglut)

design_npglut <- model.matrix(~0 + time + cell_line + npglut_fraction, data = covar_18line) #29400 genes
fit_npglut <- limmaOps(npglut_mat, 10, design_npglut)

fit_all <- list(fit_GABA, fit_nmglut, fit_npglut)
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

