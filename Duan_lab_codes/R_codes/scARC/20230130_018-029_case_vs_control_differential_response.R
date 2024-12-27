# Chuxuan Li 01/30/2023
# differential response analysis on 018-029 data
# pseudobulk count matrix by cell line and time point -> combat -> limma -> 
#residual subtraction -> linear fit again

# init ####
library(limma)
library(edgeR)
library(ggplot2)
library(stringr)
library(readr)
library(readxl)
load("./018-029_covar_table_final.RData")
load("/nvmefs/scARC_Duan_018/case_control_DE_analysis/018-029_by_type_linextime_combat_adj_mat_lst.RData")


# auxillary functions ####
createDGE <- function(count_matrix){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = celllines)
  A <- rowSums(y$counts)
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}

cnfV <- function(y, design) {
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = T)
  return(v)  
}

# make design matrices ####
covar_table_final <- covar_table_final[order(covar_table_final$cell_line), ]
covar_table_0hr <- covar_table_final[covar_table_final$time == "0hr", ]
covar_table_1hr <- covar_table_final[covar_table_final$time == "1hr", ]
covar_table_6hr <- covar_table_final[covar_table_final$time == "6hr", ]
covar_table_final$age_sq <- covar_table_final$age * covar_table_final$age

celllines <- str_remove(colnames(adj_mat_lst[[1]]), "_[0|1|6]hr$")
designs <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  type <- names(adj_mat_lst)[i]
  if (type == "GABA") {
    f = covar_table_final$GABA_fraction
    print(f)
  } else if (type == "NEFM_neg_glut") {
    f = covar_table_final$nmglut_fraction
  } else {
    f = covar_table_final$npglut_fraction
  }
  designs[[i]] <- model.matrix(~0 + batch + age_sq + sex + f, 
                               data = covar_table_final)
}

DGELists <- vector("list", length(adj_mat_lst))
voomedObjs <- vector("list", length(adj_mat_lst))
fits <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  DGELists[[i]] <- createDGE(adj_mat_lst[[i]])
  voomedObjs[[i]] <- cnfV(DGELists[[i]], designs[[i]])
  fits[[i]] <- lmFit(voomedObjs[[i]], designs[[i]])
}
names(fits) <- names(adj_mat_lst)


# extract residuals ####
resid_0hr <- vector("list", length(adj_mat_lst))
resid_1hr <- vector("list", length(adj_mat_lst))
resid_6hr <- vector("list", length(adj_mat_lst))
for (i in 1:length(fits)) {
  resid <- residuals.MArrayLM(object = fits[[i]], voomedObjs[[i]])
  resid_0hr[[i]] <- resid[, seq(1, ncol(resid), by = 3)]
  resid_1hr[[i]] <- resid[, seq(2, ncol(resid), by = 3)]
  resid_6hr[[i]] <- resid[, seq(3, ncol(resid), by = 3)]
}
names(resid_0hr) <- names(adj_mat_lst)
names(resid_1hr) <- names(adj_mat_lst)
names(resid_6hr) <- names(adj_mat_lst)
resid_1v0 <- vector("list", length(adj_mat_lst))
resid_6v0 <- vector("list", length(adj_mat_lst))
resid_6v1 <- vector("list", length(adj_mat_lst))
for (i in 1:length(resid_0hr)) {
  resid_1v0[[i]] <- resid_1hr[[i]] - resid_0hr[[i]]
  resid_6v0[[i]] <- resid_6hr[[i]] - resid_0hr[[i]]
  resid_6v1[[i]] <- resid_6hr[[i]] - resid_1hr[[i]]
}
names(resid_1v0) <- names(resid_1hr)
names(resid_6v0) <- names(resid_6hr)
names(resid_6v1) <- names(resid_6hr)

save(resid_1v0, resid_6v0, resid_6v1, file = "residual_differences_1v0_6v0_6v1.RData")

# re-fit linear model exp = b0 + b1 * aff with limma ####
designs_1v0 <- vector("list", length(resid_1v0))
designs_6v0 <- vector("list", length(resid_6v0))
designs_6v1 <- vector("list", length(resid_6v1))
for (i in 1:length(designs_1v0)) {
  designs_1v0[[i]] <- model.matrix(~0 + aff, 
                                   data = covar_table_0hr)
  designs_6v0[[i]] <- model.matrix(~0 + aff, 
                                   data = covar_table_0hr)
  designs_6v1[[i]] <- model.matrix(~0 + aff, 
                                   data = covar_table_0hr)
}
names(designs_1v0) <- names(resid_1v0)
names(designs_6v0) <- names(resid_6v0)
names(designs_6v1) <- names(resid_6v1)

contrastOnSubtractedVal <- function(df, design) {
  fit <- lmFit(df, design)
  contr.matrix <- makeContrasts(
    casevscontrol = affcase-affcontrol,
    levels = colnames(design))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit <- eBayes(fit, trend = T)
  print(summary(decideTests(fit)))
  return(fit)
}
contrasts_1v0 <- vector("list", 3)
contrasts_6v0 <- vector("list", 3)
contrasts_6v1 <- vector("list", 3)

for (i in 1:length(resid_1v0)) {
  contrasts_1v0[[i]] <- contrastOnSubtractedVal(resid_1v0[[i]], designs_1v0[[i]])
  contrasts_6v0[[i]] <- contrastOnSubtractedVal(resid_6v0[[i]], designs_6v0[[i]])
  contrasts_6v1[[i]] <- contrastOnSubtractedVal(resid_6v1[[i]], designs_6v1[[i]])
}


# loop fit linear model exp = b0 + b1 * aff with lm ####
types <- names(resid_1v0)
aff <- factor(covar_table_final$aff[covar_table_final$time == "0hr"], levels = c("control", "case"))
times <- sort(unique(covar_table_final$time))

lmfeeder <- function(row) {
  return(lm(row ~ aff))
}
cvc_1v0 <- vector("list", 3)
cvc_6v0 <- vector("list", 3)
cvc_6v1 <- vector("list", 3)
betas_1v0 <- vector("list", 3)
betas_6v0 <- vector("list", 3)
betas_6v1 <- vector("list", 3)
for (i in 1:length(resid_1v0)) {
  cvc_1v0[[i]] <- apply(X = resid_1v0[[i]], MARGIN = 1, FUN = lmfeeder)
  cvc_6v0[[i]] <- apply(X = resid_6v0[[i]], MARGIN = 1, FUN = lmfeeder)
  cvc_6v1[[i]] <- apply(X = resid_6v1[[i]], MARGIN = 1, FUN = lmfeeder)
}
for (i in 1:length(cvc_1v0)) {
  for(j in 1:length(cvc_1v0[[i]])) {
    if (j == 1) {
      mat1 <- summary(cvc_1v0[[i]][[j]])$coefficients["affcase", ]
      mat6 <- summary(cvc_6v0[[i]][[j]])$coefficients["affcase", ]
      mat61 <- summary(cvc_6v1[[i]][[j]])$coefficients["affcase", ]
    } else {
      mat1 <- rbind(mat1, summary(cvc_1v0[[i]][[j]])$coefficients["affcase", ])
      mat6 <- rbind(mat6, summary(cvc_6v0[[i]][[j]])$coefficients["affcase", ])
      mat61 <- rbind(mat61, summary(cvc_6v1[[i]][[j]])$coefficients["affcase", ])
    }
  }
  rownames(mat1) <- names(cvc_1v0[[i]])
  betas_1v0[[i]] <- as.data.frame(mat1)
  rownames(mat6) <- names(cvc_6v0[[i]])
  betas_6v0[[i]] <- as.data.frame(mat6)
  rownames(mat61) <- names(cvc_6v1[[i]])
  betas_6v1[[i]] <- as.data.frame(mat61)
}
names(betas_1v0) <- names(resid_1v0)
names(betas_6v0) <- names(resid_6v0)
names(betas_6v1) <- names(resid_6v1)

# adjust p
for (i in 1:length(betas_1v0)) {
  betas_1v0[[i]]$p.adj <- p.adjust(betas_1v0[[i]]$`Pr(>|t|)`, method = "BH")
  betas_6v0[[i]]$p.adj <- p.adjust(betas_6v0[[i]]$`Pr(>|t|)`, method = "BH")
  betas_6v1[[i]]$p.adj <- p.adjust(betas_6v1[[i]]$`Pr(>|t|)`, method = "BH")
}


# output results ####
for (i in 1:length(betas_1v0)) {
  colnames(betas_1v0[[i]]) <- c("beta", "Std. Error", "t value", "pvalue", "p.adj")
  colnames(betas_6v0[[i]]) <- c("beta", "Std. Error", "t value", "pvalue", "p.adj")
  colnames(betas_6v1[[i]]) <- c("beta", "Std. Error", "t value", "pvalue", "p.adj")
}

# summary table 
deg_counts1 <- array(dim = c(3, 4), 
                     dimnames = list(paste(types, "1v0", sep = "-"),
                                     c("positive", "negative", "nonsig", "total"))) 
for (i in 1:length(betas_1v0)) {
  deg_counts1[i,1] <- sum(betas_1v0[[i]]$pvalue < 0.05 & betas_1v0[[i]]$beta > 0)
  deg_counts1[i,2] <- sum(betas_1v0[[i]]$pvalue < 0.05 & betas_1v0[[i]]$beta < 0)
  deg_counts1[i,3] <- sum(betas_1v0[[i]]$pvalue > 0.05)
  deg_counts1[i,4] <- nrow(betas_1v0[[i]])
} 
deg_counts6 <- array(dim = c(3, 4), 
                     dimnames = list(paste(types, "6v0", sep = "-"),
                                     c("positive", "negative", "nonsig", "total"))) 
for (i in 1:length(betas_6v0)) {
  deg_counts6[i,1] <- sum(betas_6v0[[i]]$pvalue < 0.05 & betas_6v0[[i]]$beta > 0)
  deg_counts6[i,2] <- sum(betas_6v0[[i]]$pvalue < 0.05 & betas_6v0[[i]]$beta < 0)
  deg_counts6[i,3] <- sum(betas_6v0[[i]]$pvalue > 0.05)
  deg_counts6[i,4] <- nrow(betas_6v0[[i]])
} 
deg_counts <- rbind(deg_counts1, deg_counts6)
write.table(deg_counts, 
            file = "./018-029_response_as_residual_difference_casevcontrol_lm_results/DRG_counts_filter_by_pval.csv", 
            sep = ",", quote = F)

for (i in 1:length(betas_1v0)) {
  df <- betas_1v0[[i]][betas_1v0[[i]]$p.adj < 0.05 & betas_1v0[[i]]$beta > 0, ]
  filename <- 
    paste0("./018-029_response_as_residual_difference_casevcontrol_lm_results/filtered_by_padj/",
           names(betas_1v0)[i], "_1v0response_upregulated_padj.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_1v0[[i]][betas_1v0[[i]]$p.adj < 0.05 & betas_1v0[[i]]$beta < 0, ]
  filename <- 
    paste0("./018-029_response_as_residual_difference_casevcontrol_lm_results/filtered_by_padj/",
           names(betas_1v0)[i], "_1v0response_downregulated_padj.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_1v0[[i]][betas_1v0[[i]]$pvalue < 0.05 & betas_1v0[[i]]$beta > 0, ]
  filename <- 
    paste0("./018-029_response_as_residual_difference_casevcontrol_lm_results/filtered_by_pvalue/",
           names(betas_1v0)[i], "_1v0response_upregulated_pvalue.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_1v0[[i]][betas_1v0[[i]]$pvalue < 0.05 & betas_1v0[[i]]$beta < 0, ]
  filename <- 
    paste0("./018-029_response_as_residual_difference_casevcontrol_lm_results/filtered_by_pvalue/",
           names(betas_1v0)[i], "_1v0response_downregulated_pvalue.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_1v0[[i]]
  filename <-
    paste0("./018-029_response_as_residual_difference_casevcontrol_lm_results/unfiltered_full_results/",
           names(betas_1v0)[i], "_1v0response_all_genes.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
} 


# volcano plot using lm results ####
SZ_1 <- read_excel("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                   sheet = "prioritized SZ single genes")
SZ_1 <- SZ_1$gene.symbol
SZ_2 <- read_excel("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                   sheet = "SZ_SCHEMA")
SZ_2 <- SZ_2$...17
SZgenes <- union(SZ_1, SZ_2)
rm(SZ_1, SZ_2)
library(ggrepel)
for (i in 1:length(betas_1v0)){
  res <- betas_1v0[[i]]
  #res$gene.symbol <- rownames(res)
  res$significance <- "nonsignificant"
  res$significance[res$pvalue < 0.05 & res$beta > 0] <- "up (p-value < 0.05)"
  res$significance[res$pvalue < 0.05 & res$beta < 0] <- "down (p-value < 0.05)"
  unique(res$significance)
  res$significance <- factor(res$significance, levels = c("up (p-value < 0.05)",
                                                          "nonsignificant", 
                                                          "down (p-value < 0.05)"))
  res$neg_log_pval <- (0 - log2(res$pvalue))
  res$labelling <- ""
  for (j in SZgenes){
    res$labelling[res$gene %in% j] <- j
  }
  is.allnonsig <- length(unique(res$significance)) == 1
  if (is.allnonsig) {
    dotColors <- c("grey50")
  } else {
    dotColors <- c("red3", "grey50", "steelblue3")
  }
  png(paste0("./018-029_response_as_residual_difference_casevcontrol_lm_results/volcano_plots/", 
             names(betas_1v0)[i], "_casevcontrol_in_1v0_response_volcano_plot.png"), 
      width = 750, height = 750)
  p <- ggplot(data = as.data.frame(res),
              aes(x = beta, 
                  y = neg_log_pval, 
                  color = significance,
                  label = labelling)) + 
    geom_point(size = 0.2) +
    xlab("Beta") +
    ylab("-log10(p-value)") +
    scale_color_manual(values = dotColors) +
    theme_minimal() +
    geom_text_repel(box.padding = unit(0.05, 'lines'),
                    min.segment.length = 0,
                    force = 0.2,
                    max.overlaps = Inf,
                    force_pull = 0.1,
                    show.legend = F) +
    xlim(c(-2, 2)) +
    ggtitle(paste0(str_replace(str_replace(str_replace_all(names(betas_1v0)[i], 
                                                           "_", " "), " pos", "+"), 
                               " neg", "-")), "Case vs control in 1v0 response genes")
  print(p)
  dev.off()
}
