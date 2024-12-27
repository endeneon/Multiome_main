# Chuxuan Li 02/21/2023
# Differential accessibility with limma on 029 data, based on 08Aug2022_combat_limma_18line_da_peaks.R

# init ####
library(limma)
library(edgeR)
library(sva)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(readr)
library(readxl)
library(Signac)
library(Seurat)

load("/nvmefs/scARC_Duan_018/Duan_project_029_ATAC/multiomic_obj_with_new_peaks_labeled.RData")
multiomic_obj_new <- subset(multiomic_obj_new, RNA.cell.type != "unidentified")

# make covariate table ####
CIRM_table <- read_excel("~/NVME/scARC_Duan_018/CIRM control iPSC lines_40_Duan (003).xlsx")
CIRM_table <- CIRM_table[, c(4, 6, 13, 16, 21)]
CIRM_table <- CIRM_table[!is.na(CIRM_table$`Catalog ID`), ]
colnames(CIRM_table) <- c("group", "cell_line", "age", "sex", "aff")
CIRM_table$sex[CIRM_table$sex == "Female"] <- "F"
CIRM_table$sex[CIRM_table$sex == "Male"] <- "M"
CIRM_table$aff <- "control"

multiomic_obj_new$linextime.ident <- NA
lines <- unique(multiomic_obj_new$cell.line.ident)
times <- unique(multiomic_obj_new$time.ident)
for (l in lines) {
  for (t in times) {
    multiomic_obj_new$linextime.ident[multiomic_obj_new$cell.line.ident == l &
                                        multiomic_obj_new$time.ident == t] <- 
      paste0(l, "_", t)
  }
}
lts <- unique(multiomic_obj_new$linextime.ident)
med_df <- data.frame(cell_line = rep_len("", length(lts)),
                     time = rep_len("", length(lts)),
                     GABA_counts = rep_len("", length(lts)),
                     GABA_fraction = rep_len("", length(lts)),
                     nmglut_counts = rep_len("", length(lts)),
                     nmglut_fraction = rep_len("", length(lts)),
                     npglut_counts = rep_len("", length(lts)),
                     npglut_fraction = rep_len("", length(lts)),
                     total_counts = rep_len("", length(lts)))
unique(multiomic_obj_new$RNA.cell.type)
for (i in 1:length(lts)) {
  obj <- subset(multiomic_obj_new, linextime.ident == lts[[i]])
  line <- unique(obj$cell.line.ident)
  time <- unique(obj$time.ident)
  types <- unique(obj$RNA.cell.type)
  print(line)
  print(time)
  print(types)
  med_df$cell_line[i] <- line
  med_df$time[i] <- time
  med_df$GABA_counts[i] <- sum(obj$RNA.cell.type == "GABA")
  med_df$nmglut_counts[i] <- sum(obj$RNA.cell.type == "nmglut")
  med_df$npglut_counts[i] <- sum(obj$RNA.cell.type == "npglut")
  med_df$total_counts[i] <- ncol(obj)
  med_df$GABA_fraction[i] <- sum(obj$RNA.cell.type == "GABA") / ncol(obj)
  med_df$nmglut_fraction[i] <- sum(obj$RNA.cell.type == "nmglut") / ncol(obj)
  med_df$npglut_fraction[i] <- sum(obj$RNA.cell.type == "npglut") / ncol(obj)
}

# combine into covariate table
covar_df <- med_df
covar_df$age <- rep_len("", nrow(covar_df))
covar_df$sex <- rep_len("", nrow(covar_df))
covar_df$batch <- rep_len("", nrow(covar_df))
for (i in 1:length(covar_df$cell_line)) {
  line <- covar_df$cell_line[i]
  print(line)
  covar_df$age[i] <- CIRM_table$age[CIRM_table$cell_line == line]
  covar_df$sex[i] <- CIRM_table$sex[CIRM_table$cell_line == line]
  covar_df$batch[i] <- CIRM_table$group[CIRM_table$cell_line == line]
}
covar_df <- covar_df[order(covar_df$cell_line), ]
covar_df <- covar_df[order(covar_df$time), ]
covar_df$GABA_fraction <- as.numeric(covar_df$GABA_fraction)
covar_df$nmglut_fraction <- as.numeric(covar_df$nmglut_fraction)
covar_df$npglut_fraction <- as.numeric(covar_df$npglut_fraction)
covar_df$age <- as.numeric(covar_df$age)
save(covar_df, file = "029_covariates_full_df_ATAC.RData")

CIRM_table$time <- "0hr"
CIRM_table <- rbind(CIRM_table, CIRM_table, CIRM_table)
CIRM_table$time[53:104] <- "1hr"
CIRM_table$time[105:156] <- "6hr"

# define covariates from covar table
GABA_fraction <- covar_df$GABA_fraction
nmglut_fraction <- covar_df$nmglut_fraction
npglut_fraction <- covar_df$npglut_fraction

# make combat dataframes ####
GABA <- subset(multiomic_obj_new, RNA.cell.type == "GABA")
nmglut <- subset(multiomic_obj_new, RNA.cell.type == "nmglut")
npglut <- subset(multiomic_obj_new, RNA.cell.type == "npglut")

unique(nmglut$linextime.ident)
lts <- sort(unique(nmglut$linextime.ident))

# make raw count matrix ####
makeMat4CombatseqLine <- function(typeobj, linetimes) {
  colnames <- rep_len(NA, (length(linetimes)))
  for (i in 1:length(linetimes)) {
    print(linetimes[i])
    colnames[i] <- linetimes[i]
    obj <- subset(typeobj, linextime.ident == linetimes[i])
    count_matrix <- obj 
    if (i == 1) {
      mat <- as.array(rowSums(obj@assays$peaks@counts))
    } else {
      to_bind <- as.array(rowSums(obj@assays$peaks@counts))
      mat <- cbind(mat, to_bind)
    }
  }
  colnames(mat) <- colnames
  return(mat)
}

GABA_mat <- makeMat4CombatseqLine(GABA, lts)
nmglut_mat <- makeMat4CombatseqLine(nmglut, lts)
npglut_mat <- makeMat4CombatseqLine(npglut, lts)

# combat-seq ####
batch = c()
for (i in lts) {
  batch <- c(batch, unique(covar_df$batch[covar_df$cell_line == str_extract(i, "[A-Z]+[0-9]+")]))
}
batch
length(batch) #60
ncol(GABA_mat) #60
time <- str_extract(lts, "[0|1|6]hr")
GABA_mat_adj <- ComBat_seq(GABA_mat, batch = batch, group = time)
nmglut_mat_adj <- ComBat_seq(nmglut_mat, batch = batch, group = time)
npglut_mat_adj <- ComBat_seq(npglut_mat, batch = batch, group = time)
GABA_mat_adj <- GABA_mat_adj[, order(colnames(GABA_mat_adj))]
nmglut_mat_adj <- nmglut_mat_adj[, order(colnames(nmglut_mat_adj))]
npglut_mat_adj <- npglut_mat_adj[, order(colnames(npglut_mat_adj))]
save(GABA_mat_adj, nmglut_mat_adj, npglut_mat_adj, file = "post-combat_020_ATAC_count_matrices.RData")


celllines <- str_extract(lts, "[A-Z]+[0-9]+")

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
covar_df <- covar_df[order(covar_df$cell_line), ]
design_GABA <- model.matrix(~0 + time + batch + age + sex + GABA_fraction, data = covar_df) 
y_GABA <- createDGE(GABA_mat_adj)
ind.keep <- filterByCpm(y_GABA, 1, 10)
sum(ind.keep) #433918 -> 159613
v_GABA <- cnfV(y_GABA[ind.keep, ], design_GABA)
fit_GABA <- lmFit(v_GABA, design_GABA)
fit_GABA_contr <- contrastFit(fit_GABA, design_GABA)
plotSA(fit_GABA_contr, main = "Final model: Mean-variance trend") # same, just has a line

design_nmglut <- model.matrix(~0 + time + batch + age + sex + nmglut_fraction, data = covar_df) 
y_nmglut <- createDGE(nmglut_mat_adj)
ind.keep <- filterByCpm(y_nmglut, 1, 10)
sum(ind.keep) #433918 -> 163732
v_nmglut <- cnfV(y_nmglut[ind.keep, ], design_nmglut)
fit_nmglut <- lmFit(v_nmglut, design_nmglut)
fit_nmglut_contr <- contrastFit(fit_nmglut, design_nmglut)
plotSA(fit_nmglut_contr, main = "Final model: Mean-variance trend") 

design_npglut <- model.matrix(~0 + time + batch + age + sex + npglut_fraction, data = covar_df)
y_npglut <- createDGE(npglut_mat_adj)
ind.keep <- filterByCpm(y_npglut, 1, 10)
sum(ind.keep) #433918 -> 191354
v_npglut <- cnfV(y_npglut[ind.keep, ], design_npglut)
fit_npglut <- lmFit(v_npglut, design_npglut)
fit_npglut_contr <- contrastFit(fit_npglut, design_npglut)
plotSA(fit_npglut_contr, main = "Final model: Mean-variance trend") 


# output results ####
setwd("/nvmefs/scARC_Duan_018/Duan_project_029_ATAC/dapeaks_limma/")

fit_all <- list(fit_GABA_contr, fit_nmglut_contr, fit_npglut_contr)
save(fit_all, file = "029_limma_da_peaks_list_of_three_celltypes_contrast_fit_results.RData")
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
da_counts <- array(dim = c(3, 5), dimnames = list(types, c(cnames, "total"))) # 3 cell types are rows, 1/6 up/down are cols
for (i in 1:3){                                                               
  res_1v0 <- topTable(fit_all[[i]], coef = "onevszero", p.value = 0.05, 
                      sort.by = "P", number = Inf)
  res_6v0 <- topTable(fit_all[[i]], coef = "sixvszero", p.value = 0.05, 
                      sort.by = "P", number = Inf)
  da_counts[i, 1] <- sum(res_1v0$logFC  > 0)
  da_counts[i, 2] <- sum(res_1v0$logFC < 0)     
  da_counts[i, 3] <- sum(res_6v0$logFC > 0)
  da_counts[i, 4] <- sum(res_6v0$logFC < 0)
  da_counts[i, 5] <- nrow(topTable(fit_all[[i]], coef = "onevszero", p.value = Inf, 
                                   sort.by = "P", number = Inf))
}
write.table(da_counts, file = "./da_peaks_count_summary.csv", quote = F, sep = ",", 
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
    df <- data.frame(peaks = rownames(uplist),
                     p_val = uplist$P.Value,
                     avg_log2FC = uplist$logFC,
                     p_val_adj = uplist$adj.P.Val)
    file_name <- paste0("./unfiltered_bed/",
      #"./filtered_bed/",
      ident, "v0hr_upregulated.bed")
    saveBed(df, file_name)
    df <- data.frame(peaks = rownames(dolist),
                     p_val = dolist$P.Value,
                     avg_log2FC = dolist$logFC,
                     p_val_adj = dolist$adj.P.Val)
    file_name <- paste0("./unfiltered_bed/",
      #"./filtered_bed/",
      ident, "v0hr_downregulated.bed")
    saveBed(df, file_name)
    df <- data.frame(peaks = rownames(list[[i]]),
                     p_val = list[[i]]$P.Value,
                     avg_log2FC = list[[i]]$logFC,
                     p_val_adj = list[[i]]$adj.P.Val)
    file_name <- paste0("./unfiltered_bed/",
      #"./filtered_bed/",
      ident, "v0hr_full.bed")
    saveBed(df, file_name)
  }
}
