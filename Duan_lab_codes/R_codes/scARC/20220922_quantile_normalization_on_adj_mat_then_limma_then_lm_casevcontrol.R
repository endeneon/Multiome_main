# Chuxuan Li 09/22/2022
# Use quantile normalization on the combat-normalized matrix, doing it separately
#for the case and control groups, then use this result to fit linear model,
#retrieve residuals, and fit model again.

# init ####
library(limma)
library(sva)
library(edgeR)
library(readr)
library(readxl)
library(stringr)

load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/case_v_control_DE/5+18+20line_sepby_type_colby_linxetime_adj_mat_list.RData")
load("../../../covariates_pooled_together_for_all_linextime_5+18+20.RData")

# quantile normalization ####
# separate case from control
for (i in 1:length(adj_mat_lst)) {
  cols <- colnames(adj_mat_lst[[i]])
  for (j in 1:length(cols)) {
    l <- str_extract(cols[j], "CD_[0-9][0-9]")
    print(l)
    aff <- unique(covar_table$disease[covar_table$cell_line == l])
    if (j == 1) {
      aff_lst <- aff
    } else {
      aff_lst <- c(aff_lst, aff)
    }
  }
}
case_lst <- vector("list", 3)
ctrl_lst <- vector("list", 3)
for (i in 1:length(adj_mat_lst)) {
  case_lst[[i]] <- adj_mat_lst[[i]][, aff_lst == "case"]
  ctrl_lst[[i]] <- adj_mat_lst[[i]][, aff_lst == "control"]
}
# write a function
findMeanByRank <- function(ranks, means) {
  return(means[ranks])
}
quantileNormalize <- function(mat) {
  ranks <- apply(mat, 2, rank, ties.method = "min")
  mat.sorted <- apply(mat, 2, sort)
  means <- apply(mat.sorted, 1, mean)
  mat.final <- apply(ranks, 2, findMeanByRank, means = means)
  rownames(mat.final) <- rownames(mat)
  return(mat.final)
}
# normalize case and control separately
for (i in 1:length(case_lst)) {
  case_lst[[i]] <- quantileNormalize(case_lst[[i]])
  ctrl_lst[[i]] <- quantileNormalize(ctrl_lst[[i]])
}
# join case and control back together
normed_mat_lst <- vector("list", length(case_lst))
for (i in 1:length(case_lst)) {
  normed_mat_lst[[i]] <- cbind(case_lst[[i]], ctrl_lst[[i]])
}
names(normed_mat_lst) <- names(adj_mat_lst)


# make design matrices ####
covar_table <- covar_table[order(covar_table$cell_line), ]
covar_table_0hr <- covar_table[covar_table$time == "0hr", ]
covar_table_1hr <- covar_table[covar_table$time == "1hr", ]
covar_table_6hr <- covar_table[covar_table$time == "6hr", ]
celllines <- str_remove(colnames(adj_mat_lst[[1]]), "_[0|1|6]hr$")
for (i in 1:length(celllines)) {
  print(celllines[i])
  if (i == 1) {
    batch = unique(covar_table$group[covar_table$cell_line == celllines[i]])
  } else {
    batch = c(batch, unique(covar_table$group[covar_table$cell_line == celllines[i]]))
  }
}
designs <- vector("list", length(adj_mat_lst))
covar_table$age2 <- covar_table$age * covar_table$age
for (i in 1:length(adj_mat_lst)) {
  type <- names(adj_mat_lst)[i]
  if (type == "GABA") {
    f = covar_table$GABA_fraction
    print(f)
  } else if (type == "NEFM_neg_glut") {
    f = covar_table$nmglut_fraction
  } else {
    f = covar_table$npglut_fraction
  }
  designs[[i]] <- model.matrix(~0 + group + age2 + sex + f, 
                               data = covar_table)
}


# do limma voom ####
createDGE <- function(count_matrix){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = celllines)
  A <- rowSums(y$counts)
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}
filterByCpm <- function(y, cutoff, nsample) {
  cpm <- cpm(y)
  passfilter <- (rowSums(cpm >= cutoff) >= nsample)
  return(passfilter)
}
cnfV <- function(y, design) {
  y <- calcNormFactors(y)
  # voom (no logcpm needed, voom does it within the function)
  v <- voom(y, design, plot = T)
  
  return(v)  
}
contrastFit <- function(fit, design) {
  contr.matrix <- makeContrasts(
    casevscontrol = diseasecase-diseasecontrol,
    levels = colnames(design))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit <- eBayes(fit, trend = T)
  print(summary(decideTests(fit)))
  return(fit)
}
DGELists <- vector("list", length(adj_mat_lst))
voomedObjs <- vector("list", length(adj_mat_lst))
fits <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  DGELists[[i]] <- createDGE(adj_mat_lst[[i]]) # plug in adjusted or raw matrix
  ind.keep <- filterByCpm(DGELists[[i]], 1, 21)
  print(sum(ind.keep)) # 17496-17667-17006
  voomedObjs[[i]] <- cnfV(DGELists[[i]][ind.keep, ], designs[[i]])
  fits[[i]] <- lmFit(voomedObjs[[i]], designs[[i]])
}
names(fits) <- names(adj_mat_lst)

# extract residuals ####
fits_0hr <- vector("list", length(adj_mat_lst))
fits_1hr <- vector("list", length(adj_mat_lst))
fits_6hr <- vector("list", length(adj_mat_lst))
for (i in 1:length(fits)) {
  ftd <- residuals.MArrayLM(object = fits[[i]], voomedObjs[[i]])
  fits_0hr[[i]] <- ftd[, seq(1, ncol(ftd), by = 3)]
  fits_1hr[[i]] <- ftd[, seq(2, ncol(ftd), by = 3)]
  fits_6hr[[i]] <- ftd[, seq(3, ncol(ftd), by = 3)]
}
names(fits_0hr) <- names(adj_mat_lst)
names(fits_1hr) <- names(adj_mat_lst)
names(fits_6hr) <- names(adj_mat_lst)
fitted_1v0 <- vector("list", length(adj_mat_lst))
fitted_6v0 <- vector("list", length(adj_mat_lst))
fitted_6v1 <- vector("list", length(adj_mat_lst))
for (i in 1:length(fits_0hr)) {
  fitted_1v0[[i]] <- fits_1hr[[i]] - fits_0hr[[i]]
  fitted_6v0[[i]] <- fits_6hr[[i]] - fits_0hr[[i]]
  fitted_6v1[[i]] <- fits_6hr[[i]] - fits_1hr[[i]]
}
names(fitted_1v0) <- names(fits_1hr)
names(fitted_6v0) <- names(fits_6hr)
names(fitted_6v1) <- names(fits_6hr)

# test the effect of case and control directly on residuals (not response) ####
cvc_0hr <- vector("list", length(fits_0hr))
cvc_1hr <- vector("list", length(fits_1hr))
cvc_6hr <- vector("list", length(fits_6hr))
for (i in 1:length(fits_0hr)){
  d <- model.matrix(~0 + disease, data = covar_table[covar_table$time == "0hr", ])
  fit <- lmFit(fits_0hr[[i]], d)
  contr.matrix <- makeContrasts(
    casevscontrol = diseasecase-diseasecontrol,
    levels = colnames(d))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  cvc_0hr[[i]] <- eBayes(fit, trend = F)
  print(summary(decideTests(cvc_0hr[[i]])))
}

# re-fit linear model ####
disease <- factor(covar_table$disease[covar_table$time == "0hr"], levels = c("control", "case"))
lmfeeder <- function(row, celltype) {
  if (celltype == "GABA") {
    return(lm(row ~ disease))
  } else if (celltype == "NEFM_neg_glut") {
    return(lm(row ~ disease))
  } else {
    return(lm(row ~ disease))
  }
}
fits_1v0 <- vector("list", 3)
fits_6v0 <- vector("list", 3)
fits_6v1 <- vector("list", 3)
betas_1v0 <- vector("list", 3)
betas_6v0 <- vector("list", 3)
betas_6v1 <- vector("list", 3)
for (i in 1:length(fitted_1v0)) {
  t <- names(fitted_1v0)[i]
  fits_1v0[[i]] <- apply(X = fitted_1v0[[i]], MARGIN = 1, FUN = lmfeeder, celltype = t)
  t <- names(fitted_6v0)[i]
  fits_6v0[[i]] <- apply(X = fitted_6v0[[i]], MARGIN = 1, FUN = lmfeeder, celltype = t)
  t <- names(fitted_6v1)[i]
  fits_6v1[[i]] <- apply(X = fitted_6v1[[i]], MARGIN = 1, FUN = lmfeeder, celltype = t)
}
for (i in 1:length(fits_1v0)) {
  for(j in 1:length(fits_1v0[[i]])) {
    if (j == 1) {
      mat10 <- summary(fits_1v0[[i]][[j]])$coefficients["diseasecase", ]
      mat60 <- summary(fits_6v0[[i]][[j]])$coefficients["diseasecase", ]
      mat61 <- summary(fits_6v1[[i]][[j]])$coefficients["diseasecase", ]
    } else {
      mat10 <- rbind(mat10, summary(fits_1v0[[i]][[j]])$coefficients["diseasecase", ])
      mat60 <- rbind(mat60, summary(fits_6v0[[i]][[j]])$coefficients["diseasecase", ])
      mat61 <- rbind(mat61, summary(fits_6v1[[i]][[j]])$coefficients["diseasecase", ])
    }
  }
  rownames(mat10) <- names(fits_1v0[[i]])
  betas_1v0[[i]] <- as.data.frame(mat10)
  rownames(mat60) <- names(fits_6v0[[i]])
  betas_6v0[[i]] <- as.data.frame(mat60)
  rownames(mat61) <- names(fits_6v1[[i]])
  betas_6v1[[i]] <- as.data.frame(mat61)
}
names(betas_1v0) <- names(fitted_1v0)
names(betas_6v0) <- names(fitted_6v0)
names(betas_6v1) <- names(fitted_6v1)

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
for (i in 1:length(betas_6v1)) {
  df <- betas_6v1[[i]][betas_6v1[[i]]$p.adj < 0.05 & betas_6v1[[i]]$beta > 0, ]
  filename <- 
    paste0("./case_v_control_DE/qnres/positive_beta_filtered_by_padj/",
           names(betas_6v1)[i], "_differential_6v1response_genes_positiveb_sig_by_padj.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_6v1[[i]][betas_6v1[[i]]$p.adj < 0.05 & betas_6v1[[i]]$beta < 0, ]
  filename <- 
    paste0("./case_v_control_DE/qnres/negative_beta_filtered_by_padj/",
           names(betas_6v1)[i], "_differential_6v1response_genes_negativeb_sig_by_padj.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_6v1[[i]][betas_6v1[[i]]$pvalue < 0.05 & betas_6v1[[i]]$beta > 0, ]
  filename <- 
    paste0("./case_v_control_DE/qnres/positive_beta_filtered_by_pvalue/",
           names(betas_6v1)[i], "_differential_6v1response_genes_positiveb_sig_by_pvalue.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_6v1[[i]][betas_6v1[[i]]$pvalue < 0.05 & betas_6v1[[i]]$beta < 0, ]
  filename <- 
    paste0("./case_v_control_DE/qnres/negative_beta_filtered_by_pvalue/",
           names(betas_6v1)[i], "_differential_6v1response_genes_negativeb_sig_by_pvalue.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_6v1[[i]]
  filename <-
    paste0("./case_v_control_DE/qnres/unfiltered_full_results/",
           names(betas_6v1)[i], "_differential_6v1response_genes_full.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
} 

# summary table ####
types <- names(fitted_1v0)
deg_counts1 <- array(dim = c(3, 4), 
                     dimnames = list(paste(types, "1v0", sep = "-"),
                                     c("positive", "negative", "nonsig", "total"))) 
deg_counts6 <- array(dim = c(3, 4), 
                     dimnames = list(paste(types, "6v0", sep = "-"),
                                     c("positive", "negative", "nonsig", "total"))) 
deg_counts61 <- array(dim = c(3, 4), 
                      dimnames = list(paste(types, "6v1", sep = "-"),
                                      c("positive", "negative", "nonsig", "total"))) 

for (i in 1:length(betas_1v0)) {
  deg_counts1[i,1] <- sum(betas_1v0[[i]]$pvalue < 0.05 & betas_1v0[[i]]$beta > 0)
  deg_counts1[i,2] <- sum(betas_1v0[[i]]$pvalue < 0.05 & betas_1v0[[i]]$beta < 0)
  deg_counts1[i,3] <- sum(betas_1v0[[i]]$pvalue > 0.05)
  deg_counts1[i,4] <- nrow(betas_1v0[[i]])
  deg_counts6[i,1] <- sum(betas_6v0[[i]]$pvalue < 0.05 & betas_6v0[[i]]$beta > 0)
  deg_counts6[i,2] <- sum(betas_6v0[[i]]$pvalue < 0.05 & betas_6v0[[i]]$beta < 0)
  deg_counts6[i,3] <- sum(betas_6v0[[i]]$pvalue > 0.05)
  deg_counts6[i,4] <- nrow(betas_6v0[[i]])
  deg_counts61[i,1] <- sum(betas_6v1[[i]]$pvalue < 0.05 & betas_6v1[[i]]$beta > 0)
  deg_counts61[i,2] <- sum(betas_6v1[[i]]$pvalue < 0.05 & betas_6v1[[i]]$beta < 0)
  deg_counts61[i,3] <- sum(betas_6v1[[i]]$pvalue > 0.05)
  deg_counts61[i,4] <- nrow(betas_6v1[[i]])
} 
deg_counts <- rbind(deg_counts1, deg_counts6, deg_counts61)
write.table(deg_counts, 
            file = "./case_v_control_DE/qnres/casevcontrol_lm_beta_filtered_by_pvalue_summary.csv", 
            quote = F, sep = ",", row.names = T, col.names = T)

# volcano plot ####
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
  res$significance[res$pvalue < 0.05 & res$beta > 0] <- "up"
  res$significance[res$pvalue < 0.05 & res$beta < 0] <- "down"
  unique(res$significance)
  res$significance <- factor(res$significance, levels = c("up", "nonsignificant", "down"))
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
  png(paste0("./case_v_control_DE/qnres/volcano_plots/by_pvalue/", 
             names(betas_1v0)[i], "_1v0hr_volcano_plot.png"))
  p <- ggplot(data = as.data.frame(res),
              aes(x = beta, 
                  y = neg_log_pval, 
                  color = significance,
                  label = labelling)) + 
    geom_point(size = 0.2) +
    scale_color_manual(values = dotColors) +
    theme_minimal() +
    geom_text_repel(box.padding = unit(0.05, 'lines'),
                    min.segment.length = 0,
                    force = 0.2,
                    max.overlaps = Inf,
                    force_pull = 0.1,
                    show.legend = F) +
    xlim(c(-2, 2)) +
    ggtitle(str_replace(str_replace(str_replace_all(names(betas_1v0)[i], "_", " "), " pos", "+"), " neg", "-"))
  print(p)
  dev.off()
}

