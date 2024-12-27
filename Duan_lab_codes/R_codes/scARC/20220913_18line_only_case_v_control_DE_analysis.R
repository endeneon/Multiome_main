# Chuxuan li 09/12/2022
# Use 18-line data only to compare case and control against residual of linear
#model

# init ####
library(stringr)
library(readr)
library(Seurat)
library(limma)
library(edgeR)
library(DESeq2)
library(sva)

load("./labeled_nfeat3000.RData")
load("./pseudobulk_DE/res_use_combatseq_mat/by_line_pseudobulk_mat_adjusted.RData")
load("../covariates_pooled_together_for_all_linextime_5+18+20.RData")
adj_mat_lst <- list(GABA_mat_adj, nmglut_mat_adj, npglut_mat_adj)
names(adj_mat_lst) <- c("GABA", "nmglut", "npglut")
celllines <- sort(unique(integrated_labeled$cell.line.ident))
covar_table <- covar_table[covar_table$cell_line %in% celllines, ]

# auxillary functions ####
createDGE <- function(count_matrix){
  y <- DGEList(counts = count_matrix, 
               genes = rownames(count_matrix), group = colnames(count_matrix))
  A <- rowSums(y$counts)
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

# Limma with voom####
covar_table <- covar_table[order(covar_table$cell_line), ]
designs <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  type <- str_split(names(adj_mat_lst)[i], pattern = "-", n = 2, simplify = T)[,1]
  if (type == "GABA") {
    f = covar_table$GABA_fraction
  } else if (type == "NEFM_neg_glut") {
    f = covar_table$nmglut_fraction
  } else {
    f = covar_table$npglut_fraction
  }
  designs[[i]] <- model.matrix(~0 + group + age^2 + sex + f, 
                               data = covar_table)
}

DGELists <- vector("list", length(adj_mat_lst))
voomedObjs <- vector("list", length(adj_mat_lst))
fits <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  DGELists[[i]] <- createDGE(adj_mat_lst[[i]]) # plug in adjusted or raw matrix
  ind.keep <- filterByCpm(DGELists[[i]], 1, 9)
  print(sum(ind.keep)) # 17547-17110-17557
  voomedObjs[[i]] <- cnfV(DGELists[[i]][ind.keep, ], designs[[i]])
  fits[[i]] <- lmFit(voomedObjs[[i]], designs[[i]])
}
names(fits) <- names(adj_mat_lst)

# extract fitted values ####
fits_0hr <- vector("list", length(adj_mat_lst))
fits_1hr <- vector("list", length(adj_mat_lst))
fits_6hr <- vector("list", length(adj_mat_lst))
for (i in 1:length(fits)) {
  ftd <- residuals.MArrayLM(fits[[i]], voomedObjs[[i]])
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

# re-fit linear model ####
designs_1v0 <- vector("list", length(fitted_1v0))
designs_6v0 <- vector("list", length(fitted_6v0))

for (i in 1:length(designs_1v0)) {
 designs_1v0[[i]] <- model.matrix(~0 + disease, 
                                   data = covar_table[covar_table$time == "0hr", ])
designs_6v0[[i]] <- model.matrix(~0 + disease, 
                                   data = covar_table[covar_table$time == "0hr", ])
}
names(designs_1v0) <- names(fitted_1v0)
names(designs_6v0) <- names(fitted_6v0)

contrastOnSubtractedVal <- function(df, design) {
  fit <- lmFit(df, design)
  contr.matrix <- makeContrasts(
    casevscontrol = diseasecase-diseasecontrol,
    levels = colnames(design))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit <- eBayes(fit, trend = T)
  print(summary(decideTests(fit)))
  return(fit)
}
contrasts_1v0 <- vector("list", 3)
contrasts_6v0 <- vector("list", 3)

for (i in 1:length(fitted_1v0)) {
  contrasts_1v0[[i]] <- contrastOnSubtractedVal(fitted_1v0[[i]], designs_1v0[[i]])
  contrasts_6v0[[i]] <- contrastOnSubtractedVal(fitted_6v0[[i]], designs_6v0[[i]])
}

# loop fit ####
types <- names(fitted_1v0)
disease <- factor(covar_table$disease[covar_table$time == "0hr"], levels = c("control", "case"))
age <- covar_table$age[covar_table$time == "0hr"]
sex <- covar_table$sex[covar_table$time == "0hr"]
group <- covar_table$group[covar_table$time == "0hr"]
fractions_GABA <- vector("list", 3)
fractions_nmglut <- vector("list", 3)
fractions_npglut <- vector("list", 3)
times <- sort(unique(covar_table$time))
for (i in 1:length(fractions_GABA)) {
  fractions_GABA[[i]] <- as.array(covar_table$GABA_fraction[covar_table$time == times[i]])
  fractions_nmglut[[i]] <- as.array(covar_table$nmglut_fraction[covar_table$time == times[i]])
  fractions_npglut[[i]] <- as.array(covar_table$npglut_fraction[covar_table$time == times[i]])
}

lmfeeder <- function(row) {
  return(lm(row ~ disease))
}
fits_1v0 <- vector("list", 3)
fits_6v0 <- vector("list", 3)
fits_6v1 <- vector("list", 3)
betas_1v0 <- vector("list", 3)
betas_6v0 <- vector("list", 3)
betas_6v1 <- vector("list", 3)
for (i in 1:length(fitted_1v0)) {
  fits_1v0[[i]] <- apply(X = fitted_1v0[[i]], MARGIN = 1, FUN = lmfeeder)
  fits_6v0[[i]] <- apply(X = fitted_6v0[[i]], MARGIN = 1, FUN = lmfeeder)
  fits_6v1[[i]] <- apply(X = fitted_6v1[[i]], MARGIN = 1, FUN = lmfeeder)
}
for (i in 1:length(fits_1v0)) {
  for(j in 1:length(fits_1v0[[i]])) {
    if (j == 1) {
      mat1 <- summary(fits_1v0[[i]][[j]])$coefficients["diseasecase", ]
      mat6 <- summary(fits_6v0[[i]][[j]])$coefficients["diseasecase", ]
      mat61 <- summary(fits_6v1[[i]][[j]])$coefficients["diseasecase", ]
    } else {
      mat1 <- rbind(mat1, summary(fits_1v0[[i]][[j]])$coefficients["diseasecase", ])
      mat6 <- rbind(mat6, summary(fits_6v0[[i]][[j]])$coefficients["diseasecase", ])
      mat61 <- rbind(mat61, summary(fits_6v1[[i]][[j]])$coefficients["diseasecase", ])
    }
  }
  rownames(mat1) <- names(fits_1v0[[i]])
  betas_1v0[[i]] <- as.data.frame(mat1)
  rownames(mat6) <- names(fits_6v0[[i]])
  betas_6v0[[i]] <- as.data.frame(mat6)
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

# summary table 
deg_counts1 <- array(dim = c(3, 4), 
                     dimnames = list(paste(types, "1v0", sep = "-"),
                                     c("positive", "negative", "nonsig", "total"))) 
for (i in 1:length(betas_1v0)) {
  deg_counts1[i,1] <- sum(betas_1v0[[i]]$`Pr(>|t|)` < 0.05 & betas_1v0[[i]]$Estimate > 0)
  deg_counts1[i,2] <- sum(betas_1v0[[i]]$`Pr(>|t|)` < 0.05 & betas_1v0[[i]]$Estimate < 0)
  deg_counts1[i,3] <- sum(betas_1v0[[i]]$`Pr(>|t|)` > 0.05)
  deg_counts1[i,4] <- nrow(betas_1v0[[i]])
} 
deg_counts6 <- array(dim = c(3, 4), 
                     dimnames = list(paste(types, "6v0", sep = "-"),
                                     c("positive", "negative", "nonsig", "total"))) 
for (i in 1:length(betas_6v0)) {
  deg_counts6[i,1] <- sum(betas_6v0[[i]]$`Pr(>|t|)` < 0.05 & betas_6v0[[i]]$Estimate > 0)
  deg_counts6[i,2] <- sum(betas_6v0[[i]]$`Pr(>|t|)` < 0.05 & betas_6v0[[i]]$Estimate < 0)
  deg_counts6[i,3] <- sum(betas_6v0[[i]]$`Pr(>|t|)` > 0.05)
  deg_counts6[i,4] <- nrow(betas_6v0[[i]])
} 
deg_counts <- rbind(deg_counts1, deg_counts6)
