# Chuxuan Li 05/05/2023
# compare gene expression in 2 conditions in rabbit cells using limma

# init ####
# init ####
library(Seurat)
library(limma)
library(edgeR)
library(sva)

library(readr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# make count df ####
load("./oc2_mm10_analysis/integrated_labeled_after_filtering2.Rdata")

# make pseudobulk samples ####
DefaultAssay(integrated_labeled) <- "RNA"
libs <- unique(integrated_labeled$orig.ident)
Idents(integrated_labeled)
times <- unique(integrated_labeled$time.ident)

# GABA 
GABA <- subset(integrated_labeled, cell.type == "GABA")
# count number of cells in each line, type, and time
for (i in 1:length(libs)) {
  print(libs[i])
  temp_count <- rowSums(subset(GABA, orig.ident == libs[i]))
  if (i == 1) {
    df_GABA <- data.frame(temp_count,
                            stringsAsFactors = F)
    colnames(df_GABA) <- libs[i]
  } else {
    temp_colnames <- colnames(df_GABA)
    df_GABA <- data.frame(cbind(df_GABA,
                                  temp_count),
                            stringsAsFactors = F)
    colnames(df_GABA) <- c(temp_colnames, libs[i])
  }
}
colnames(df_GABA) 
rownames(df_GABA)

# glut 
glut <- subset(integrated_labeled, cell.type == "glut")
# count number of cells in each line, type, and time
for (i in 1:length(libs)) {
  print(libs[i])
  temp_count <- rowSums(subset(glut, orig.ident == libs[i]))
  if (i == 1) {
    df_glut <- data.frame(temp_count,
                          stringsAsFactors = F)
    colnames(df_glut) <- libs[i]
  } else {
    temp_colnames <- colnames(df_glut)
    df_glut <- data.frame(cbind(df_glut,
                                temp_count),
                          stringsAsFactors = F)
    colnames(df_glut) <- c(temp_colnames, libs[i])
  }
}
colnames(df_glut) 
rownames(df_glut)

# corneu 
corneu <- subset(integrated_labeled, cell.type %in% c("GABA", "glut"))
# count number of cells in each line, type, and time
for (i in 1:length(libs)) {
  print(libs[i])
  temp_count <- rowSums(subset(corneu, orig.ident == libs[i]))
  if (i == 1) {
    df_corneu <- data.frame(temp_count,
                          stringsAsFactors = F)
    colnames(df_corneu) <- libs[i]
  } else {
    temp_colnames <- colnames(df_corneu)
    df_corneu <- data.frame(cbind(df_corneu,
                                temp_count),
                          stringsAsFactors = F)
    colnames(df_corneu) <- c(temp_colnames, libs[i])
  }
}
colnames(df_corneu) 
rownames(df_corneu)

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
  
  cpm_1 <- cpm[, 1:3]
  cpm_2 <- cpm[, 4:6]
  passfilter <- (rowSums(cpm_1 >= cutoff) >= nsample |
                   rowSums(cpm_2 >= cutoff) >= nsample)
  return(passfilter)
}
cnfV <- function(y, design) {
  y <- calcNormFactors(y)
  # voom (no logcpm)
  v <- voom(y, design, plot = T)
  
  return(v)  
}
contrastFit <- function(fit, design) {
  contr.matrix <- makeContrasts(cond1vs2 = condcondition1-condcondition2,
                                levels = colnames(design))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit <- eBayes(fit, trend = T)
  print(summary(decideTests(fit)))
  return(fit)
}

# run limma on GABA ####
covar_df <- data.frame(lib = libs,
                       cond = c("condition1", "condition1", "condition1",
                                "condition2", "condition2", "condition2"))
design <- model.matrix(~0 + cond, data = covar_df)
y <- createDGE(df_GABA) # plug in adjusted or raw matrix
ind.keep <- filterByCpm(y, 1, 3)
sum(ind.keep) #19645
v <- cnfV(y[ind.keep, ], design)
fit <- lmFit(v, design)
fit_contr <- contrastFit(fit, design)
plotSA(fit_contr, main = "Final model: Mean-variance trend") # same, just has a line

# output results
setwd("./oc2_mm10_analysis/limma_DE_results/GABA/")
res <- topTable(fit_contr, coef = "cond1vs2", p.value = Inf, sort.by = "P", number = Inf)
write.table(res[res$logFC > 0 & res$P.Value < 0.05, ], 
            file = "GABA_upregulated_genes_in_condition1_pvalue_significant.csv",
            quote = F, sep = ",", row.names = F)
write.table(res[res$logFC < 0 & res$P.Value < 0.05, ], 
            file = "GABA_downregulated_genes_in_condition1_pvalue_significant.csv",
            quote = F, sep = ",", row.names = F)
write.table(res[res$logFC > 0, ], 
            file = "GABA_upregulated_genes_in_condition1_full.csv",
            quote = F, sep = ",", row.names = F)
write.table(res[res$logFC < 0, ], 
            file = "GABA_downregulated_genes_in_condition1_full.csv",
            quote = F, sep = ",", row.names = F)

# run limma on glut ####
covar_df <- data.frame(lib = libs,
                       cond = c("condition1", "condition1", "condition1",
                                "condition2", "condition2", "condition2"))
design <- model.matrix(~0 + cond, data = covar_df)
y <- createDGE(df_glut) # plug in adjusted or raw matrix
ind.keep <- filterByCpm(y, 1, 3)
sum(ind.keep) #19221
v <- cnfV(y[ind.keep, ], design)
fit <- lmFit(v, design)
fit_contr <- contrastFit(fit, design)
plotSA(fit_contr, main = "Final model: Mean-variance trend") # same, just has a line

# output results
setwd("../glut/")
res <- topTable(fit_contr, coef = "cond1vs2", p.value = Inf, sort.by = "P", number = Inf)
write.table(res[res$logFC > 0 & res$P.Value < 0.05, ], 
            file = "glut_upregulated_genes_in_condition1_pvalue_significant.csv",
            quote = F, sep = ",", row.names = F)
write.table(res[res$logFC < 0 & res$P.Value < 0.05, ], 
            file = "glut_downregulated_genes_in_condition1_pvalue_significant.csv",
            quote = F, sep = ",", row.names = F)
write.table(res[res$logFC > 0, ], 
            file = "glut_upregulated_genes_in_condition1_full.csv",
            quote = F, sep = ",", row.names = F)
write.table(res[res$logFC < 0, ], 
            file = "glut_downregulated_genes_in_condition1_full.csv",
            quote = F, sep = ",", row.names = F)

# run limma on cortical neurons ####
covar_df <- data.frame(lib = libs,
                       cond = c("condition1", "condition1", "condition1",
                                "condition2", "condition2", "condition2"))
design <- model.matrix(~0 + cond, data = covar_df)
y <- createDGE(df_corneu) # plug in adjusted or raw matrix
ind.keep <- filterByCpm(y, 1, 3)
sum(ind.keep) #20227
v <- cnfV(y[ind.keep, ], design)
fit <- lmFit(v, design)
fit_contr <- contrastFit(fit, design)
plotSA(fit_contr, main = "Final model: Mean-variance trend") # same, just has a line

# output results
setwd("./oc2_mm10_analysis/limma_DE_results/GABA_glut_combined/")
res <- topTable(fit_contr, coef = "cond1vs2", p.value = Inf, sort.by = "P", number = Inf)
write.table(res[res$logFC > 0 & res$P.Value < 0.05, ], 
            file = "corneu_upregulated_genes_in_condition1_pvalue_significant.csv",
            quote = F, sep = ",", row.names = F)
write.table(res[res$logFC < 0 & res$P.Value < 0.05, ], 
            file = "corneu_downregulated_genes_in_condition1_pvalue_significant.csv",
            quote = F, sep = ",", row.names = F)
write.table(res[res$logFC > 0, ], 
            file = "corneu_upregulated_genes_in_condition1_full.csv",
            quote = F, sep = ",", row.names = F)
write.table(res[res$logFC < 0, ], 
            file = "corneu_downregulated_genes_in_condition1_full.csv",
            quote = F, sep = ",", row.names = F)

# summary table 
cnames <- c("up in condition 1 (20088)", "down in condition 1 (20088)", "total passed filter")
deg_counts <- array(dim = c(2, 3), 
                    dimnames = list(c("FDR < 0.05", "p-value < 0.05"), cnames)) # 3 cell types are rows, 1/6 up/down are cols
deg_counts[1, 1] <- sum(res$logFC > 0 & res$adj.P.Val < 0.05)
deg_counts[1, 2] <- sum(res$logFC < 0 & res$adj.P.Val < 0.05)     
deg_counts[1, 3] <- length(res$logFC)
deg_counts[2, 1] <- sum(res$logFC > 0 & res$P.Value < 0.05)
deg_counts[2, 2] <- sum(res$logFC < 0 & res$P.Value < 0.05)     
deg_counts[2, 3] <- length(res$logFC)
write.table(deg_counts, file = "DEG_counts_rabbit.csv", quote = F, sep = ",")

# volcano plots
res$significance <- "nonsignificant"
res$significance[res$P.Value < 0.05 & res$logFC > 0] <- "up (p < 0.05)"
res$significance[res$P.Value < 0.05 & res$logFC < 0] <- "down (p < 0.05)"
unique(res$significance)
res$significance <- factor(res$significance, levels = c("up (p < 0.05)", 
                                                        "nonsignificant", "down (p < 0.05)"))
res$neg_log_pval <- (0 - log2(res$P.Value))
png("rabbit_experiment_volcano_plot.png", width = 400, height = 400)
ggplot(data = as.data.frame(res),
       aes(x = logFC, 
           y = neg_log_pval, 
           color = significance)) + 
  geom_point(size = 0.2) +
  scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
  theme_minimal() +
  xlim(c(-5, 5)) +
  ggtitle("Differential expression in rabbit experiment") +
  ylab("-log10(p)")
dev.off()
