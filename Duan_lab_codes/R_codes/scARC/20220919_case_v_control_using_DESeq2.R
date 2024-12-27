# Chuxuan Li 09/19/2022
# Test case v control difference at 0,1,6 and response difference, using DESeq2

# init ####
library(DESeq2)
library(stringr)
library(readr)
library(ggplot2)
library(readxl)

load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/case_v_control_DE/5+18+20line_sepby_type_colby_linxetime_adj_mat_list.RData")
load("../../../covariates_pooled_together_for_all_linextime_5+18+20.RData")
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/5+18+20line_residuals_matrices_w_agesqd_notime.RData")


# 1. test case v control at 0, 1, 6 ####
# separate 3 time points
GABA_3times <- vector("list", 3)
nmglut_3times <- vector("list", 3)
npglut_3times <- vector("list", 3)
times = sort(unique(str_extract(string = colnames(adj_mat_lst$GABA),
                   pattern = "[0|1|6]hr")))
for (i in 1:3) {
  GABA_3times[[i]] <- adj_mat_lst$GABA[, str_detect(colnames(adj_mat_lst$GABA), times[i])]
  nmglut_3times[[i]] <- adj_mat_lst$GABA[, str_detect(colnames(adj_mat_lst$NEFM_neg_glut), times[i])]
  npglut_3times[[i]] <- adj_mat_lst$GABA[, str_detect(colnames(adj_mat_lst$NEFM_pos_glut), times[i])]
}
names(GABA_3times) <- times
names(nmglut_3times) <- times
names(npglut_3times) <- times

# create coldata df
celllines <- sort(str_remove(colnames(GABA_3times[[1]]), "_[0|1|6]hr$"))

for (i in 1:length(celllines)) {
  print(celllines[i])
  if (i == 1) {
    batch = unique(covar_table$group[covar_table$cell_line == celllines[i]])
  } else {
    batch = c(batch, unique(covar_table$group[covar_table$cell_line == celllines[i]]))
  }
}
batch <- as.factor(batch)

for (i in 1:length(celllines)) {
  print(celllines[i])
  if (i == 1) {
    age = unique(covar_table$age[covar_table$cell_line == celllines[i]])
    sex = unique(covar_table$sex[covar_table$cell_line == celllines[i]])
    disease = unique(covar_table$disease[covar_table$cell_line == celllines[i]])
  } else {
    age = c(age, unique(covar_table$age[covar_table$cell_line == celllines[i]]))
    sex = c(sex, unique(covar_table$sex[covar_table$cell_line == celllines[i]]))
    disease = c(disease, unique(covar_table$disease[covar_table$cell_line == celllines[i]]))
  }
}
# test GABA
reslist_GABA <- vector("list", 3)
for (i in 1:3){
  GABA_fraction <- covar_table$GABA_fraction[covar_table$time == names(GABA_3times)[i]]
  coldata_GABA <- data.frame(batch = batch,
                             age = age,
                             aff = disease,
                             sex = sex,
                             GABA_fraction = GABA_fraction)
  rownames(coldata_GABA) <- colnames(GABA_3times[[i]])
  dds_GABA <- DESeqDataSetFromMatrix(countData = GABA_3times[[i]],
                                     colData = coldata_GABA,
                                     design = ~ aff + batch + age + sex + GABA_fraction)  
  keep <- rowSums(counts(dds_GABA)) > 500
  dds_GABA <- dds_GABA[keep,]
  
  # de analysis
  dds_GABA <- DESeq(dds_GABA)
  res_GABA <- results(dds_GABA, contrast = c("aff", "case", "control"))
  reslist_GABA[[i]] <- res_GABA[res_GABA$padj < 0.1, ]
  write.table(x = rownames(res_GABA[res_GABA$pvalue < 0.05, ]), 
              file = paste0("./case_v_control_DE/DESeq2_res/GABA_genelist_", times[i], ".txt"),
              quote = F, sep = ",", row.names = F, col.names = F)
}

# test nmglut
reslist_nmglut <- vector("list", 3)
for (i in 1:3){
  nmglut_fraction <- covar_table$nmglut_fraction[covar_table$time == names(nmglut_3times)[i]]
  coldata_nmglut <- data.frame(batch = batch,
                             age = age,
                             aff = disease,
                             sex = sex,
                             nmglut_fraction = nmglut_fraction)
  rownames(coldata_nmglut) <- colnames(nmglut_3times[[i]])
  dds_nmglut <- DESeqDataSetFromMatrix(countData = nmglut_3times[[i]],
                                     colData = coldata_nmglut,
                                     design = ~ aff + batch + age + sex + nmglut_fraction)  
  keep <- rowSums(counts(dds_nmglut)) > 500
  dds_nmglut <- dds_nmglut[keep,]
  
  # de analysis
  dds_nmglut <- DESeq(dds_nmglut)
  res_nmglut <- results(dds_nmglut, contrast = c("aff", "case", "control"))
  reslist_nmglut[[i]] <- res_nmglut[res_nmglut$padj < 0.1, ]
  write.table(x = rownames(res_nmglut[res_nmglut$pvalue < 0.05, ]), 
              file = paste0("./case_v_control_DE/DESeq2_res/nmglut_genelist_", times[i], ".txt"),
              quote = F, sep = ",", row.names = F, col.names = F)
  
}

# test npglut
reslist_npglut <- vector("list", 3)
for (i in 1:3){
  npglut_fraction <- covar_table$npglut_fraction[covar_table$time == names(npglut_3times)[i]]
  coldata_npglut <- data.frame(batch = batch,
                               age = age,
                               aff = disease,
                               sex = sex,
                               npglut_fraction = npglut_fraction)
  rownames(coldata_npglut) <- colnames(npglut_3times[[i]])
  dds_npglut <- DESeqDataSetFromMatrix(countData = npglut_3times[[i]],
                                       colData = coldata_npglut,
                                       design = ~ aff + batch + age + sex + npglut_fraction)  
  keep <- rowSums(counts(dds_npglut)) > 500
  dds_npglut <- dds_npglut[keep,]
  
  # de analysis
  dds_npglut <- DESeq(dds_npglut)
  res_npglut <- results(dds_npglut, contrast = c("aff", "case", "control"))
  reslist_npglut[[i]] <- res_npglut[res_npglut$padj < 0.1, ]
  write.table(x = rownames(res_npglut[res_npglut$pvalue < 0.05, ]), 
              file = paste0("./case_v_control_DE/DESeq2_res/npglut_genelist_", times[i], ".txt"),
              quote = F, sep = ",", row.names = F, col.names = F)
}

# 2. test on response ####
for (i in 1:length(fitted_1v0)){
  genes <- rownames(fitted_1v0[[i]])
  fitted_1v0[[i]] <- apply(apply(fitted_1v0[[i]], 2, exp), 2, as.integer)
  rownames(fitted_1v0[[i]]) <- genes
  genes <- rownames(fitted_6v0[[i]])
  fitted_6v0[[i]] <- apply(apply(fitted_6v0[[i]], 2, exp), 2, as.integer)
  rownames(fitted_6v0[[i]]) <- genes
  genes <- rownames(fitted_6v1[[i]])
  fitted_6v1[[i]] <- apply(apply(fitted_6v1[[i]], 2, exp), 2, as.integer)
  rownames(fitted_6v1[[i]]) <- genes
}
responses <- list(fitted_1v0, fitted_6v0, fitted_6v1)
reslist_GABA <- vector("list", 3)
reslist_nmglut <- vector("list", 3)
reslist_npglut <- vector("list", 3)

for (i in 1:length(responses)) {
  # same coldata for 3 cell types
  coldata <- data.frame(aff = disease)
  rownames(coldata) <- colnames(responses[[i]]$GABA)
  # GABA
  dds_GABA <- DESeqDataSetFromMatrix(countData = responses[[i]]$GABA,
                                     colData = coldata,
                                     design = ~ aff)  
  keep <- rowSums(counts(dds_GABA)) > 10
  dds_GABA <- dds_GABA[keep,]
  dds_GABA <- DESeq(dds_GABA)
  res_GABA <- results(dds_GABA, contrast = c("aff", "case", "control"))
  res_GABA <- res_GABA[!is.na(res_GABA$padj), ]
  #reslist_GABA[[i]] <- res_GABA[res_GABA$padj < 0.05, ]
  reslist_GABA[[i]] <- res_GABA
  # nmglut
  dds_nmglut <- DESeqDataSetFromMatrix(countData = responses[[i]]$NEFM_neg_glut,
                                     colData = coldata,
                                     design = ~ aff)  
  keep <- rowSums(counts(dds_nmglut)) > 10
  dds_nmglut <- dds_nmglut[keep,]
  dds_nmglut <- DESeq(dds_nmglut)
  res_nmglut <- results(dds_nmglut, contrast = c("aff", "case", "control"))
  res_nmglut <- res_nmglut[!is.na(res_nmglut$padj), ]
  #reslist_nmglut[[i]] <- res_nmglut[res_nmglut$padj < 0.05, ]
  reslist_nmglut[[i]] <- res_nmglut
  # npglut
  dds_npglut <- DESeqDataSetFromMatrix(countData = responses[[i]]$NEFM_pos_glut,
                                       colData = coldata,
                                       design = ~ aff)  
  keep <- rowSums(counts(dds_npglut)) > 10
  dds_npglut <- dds_npglut[keep,]
  dds_npglut <- DESeq(dds_npglut)
  res_npglut <- results(dds_npglut, contrast = c("aff", "case", "control"))
  res_npglut <- res_npglut[!is.na(res_npglut$padj), ]
  #reslist_npglut[[i]] <- res_npglut[res_npglut$padj < 0.05, ]
  reslist_npglut[[i]] <- res_npglut
}
names(reslist_GABA) <- c("1v0", "6v0", "6v1")
names(reslist_nmglut) <- c("1v0", "6v0", "6v1")
names(reslist_npglut) <- c("1v0", "6v0", "6v1")
writeGenelistTxt <- function(dflist, celltype) {
  tnames <- names(dflist)
  for (i in 1:length(dflist)) {
    print(tnames[i])
    df <- dflist[[i]]
    write.table(x = rownames(df[df$pvalue < 0.05 & df$log2FoldChange > 0, ]), 
                file = paste0("./case_v_control_DE/DESeq2_res/", celltype, 
                              "_", tnames[i], "_genelist_up.txt"),
                quote = F, sep = ",", row.names = F, col.names = F)
    write.table(x = rownames(df[df$pvalue < 0.05 & df$log2FoldChange < 0, ]), 
                file = paste0("./case_v_control_DE/DESeq2_res/", celltype, 
                              "_", tnames[i], "_genelist_down.txt"),
                quote = F, sep = ",", row.names = F, col.names = F)
  }
}
writeGenelistTxt(reslist_GABA, "GABA")
writeGenelistTxt(reslist_nmglut, "nmglut")
writeGenelistTxt(reslist_npglut, "npglut")


# test enrichment on differential response genes ####

# SZ genes
SZ_1 <- read_excel("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                   sheet = "prioritized SZ single genes")
SZ_1 <- SZ_1$gene.symbol
SZ_2 <- read_excel("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                   sheet = "SZ_SCHEMA")
SZ_2 <- SZ_2$...17
SZgenes <- union(SZ_1, SZ_2)
rm(SZ_1, SZ_2)

fishersExactTest <- function(dds, dflist) {
  reslist <- vector("list", length(dflist))
  for (i in 1:length(dflist)) {
    df <- dflist[[i]]
    # upregulated 
    up_sig_genes <- rownames(df[df$padj < 0.05 & df$log2FoldChange > 0, ])
    df4test <- data.frame(all = c(nrow(dds), 
                                  length(intersect(SZgenes, rownames(dds)))),
                          sig = c(length(up_sig_genes),
                                  length(intersect(SZgenes, up_sig_genes))),
                          row.names = c("nonSZ", "SZ"))
    up_res <- fisher.test(df4test, alternative = "g")
    # downregulated 
    down_sig_genes <- rownames(df[df$padj < 0.05 & df$log2FoldChange < 0, ])
    df4test <- data.frame(all = c(nrow(dds), 
                                  length(intersect(SZgenes, rownames(dds)))),
                          sig = c(length(down_sig_genes),
                                  length(intersect(SZgenes, down_sig_genes))),
                          row.names = c("nonSZ", "SZ"))
    down_res <- fisher.test(df4test, alternative = "g")
    reslist[[i]] <- list(up_res, down_res)
  }
  return(reslist)
}
GABA_fisher <- fishersExactTest(dds_GABA, reslist_GABA)
nmglut_fisher <- fishersExactTest(dds_nmglut, reslist_nmglut)
npglut_fisher <- fishersExactTest(dds_npglut, reslist_npglut)


# plot dotplots for significant genes ####
for (i in 1:length(reslist_GABA)) {
  reslist_GABA[[i]]$gene <- rownames(reslist_GABA[[i]])
  reslist_nmglut[[i]]$gene <- rownames(reslist_nmglut[[i]])
  reslist_npglut[[i]]$gene <- rownames(reslist_npglut[[i]])
  df4plot <- reslist_GABA[[i]][reslist_GABA[[i]]$padj < 0.05, ]
  df4plot$cell.type <- "GABA"
  df1 <- reslist_nmglut[[i]][reslist_nmglut[[i]]$padj < 0.05, ]
  df1$cell.type <- "NEFM- glut"
  df2 <- reslist_npglut[[i]][reslist_npglut[[i]]$padj < 0.05, ]
  df2$cell.type <- "NEFM+ glut"
  df4plot <- rbind(df4plot, df1, df2)
  df4plot <- as.data.frame(df4plot)
  df4plot$neglogp <- (-1) * log2(df4plot$padj)
  title <- paste0("Differential response genes - ", names(reslist_GABA)[i], "hr")
  print(title)
  filename <- paste0("./case_v_control_DE/DESeq2_res/DRG_", names(reslist_GABA)[i], "hr.pdf")
  pdf(file = filename, width = 12, height = 5)
  p <- ggplot(df4plot, aes(x = gene, y = cell.type,
                      color = log2FoldChange, size = neglogp)) +
    geom_point() +
    scale_color_gradientn(colours = c("royalblue4", "white", "darkred")) +
    theme_light() +
    theme(text = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(size = "-log2(p-value)", x = "Genes", y = "Cell Type") +
    ggtitle(title)
  print(p)
  dev.off()
}

# output gene lists ####
for (i in 1:length(reslist_GABA)) {
  reslist_GABA[[i]]$gene <- rownames(reslist_GABA[[i]])
  reslist_nmglut[[i]]$gene <- rownames(reslist_nmglut[[i]])
  reslist_npglut[[i]]$gene <- rownames(reslist_npglut[[i]])
  df4plot <- reslist_GABA[[i]][reslist_GABA[[i]]$padj < 0.05, ]
  df4plot$cell.type <- "GABA"
  df1 <- reslist_nmglut[[i]][reslist_nmglut[[i]]$padj < 0.05, ]
  df1$cell.type <- "NEFM- glut"
  df2 <- reslist_npglut[[i]][reslist_npglut[[i]]$padj < 0.05, ]
  df2$cell.type <- "NEFM+ glut"
  df4plot <- rbind(df4plot, df1, df2)
  write.table(df4plot, 
              file = paste0("./case_v_control_DE/DESeq2_res/padj_sig_genes_df/DESeq2_DRGs_", 
                            names(reslist_GABA)[i], "hr.csv"),
              sep = ",", quote = F, row.names = T, col.names = T)
}
