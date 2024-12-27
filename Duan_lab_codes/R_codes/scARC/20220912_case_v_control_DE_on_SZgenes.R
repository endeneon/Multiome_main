# Chuxuan Li 08/26/2022
# Do 06Aug2022_use_residuals_case_control analysis on SZ genes only


# init ####
library(limma)
library(edgeR)
library(ggplot2)
library(stringr)
library(readr)
library(readxl)

load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/case_v_control_DE/5+18+20line_sepby_type_colby_linxetime_adj_mat_list.RData")
load("../../../covariates_pooled_together_for_all_linextime_5+18+20.RData")

SZ_1 <- read_excel("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                   sheet = "prioritized SZ single genes")
SZ_1 <- SZ_1$gene.symbol
SZ_2 <- read_excel("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/SZ_BP_ASD_MD_PTSD_genes_updated.xlsx", 
                   sheet = "SZ_SCHEMA")
SZ_2 <- SZ_2$...17
SZgenes <- union(SZ_1, SZ_2)
rm(SZ_1, SZ_2)

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
celllines <- str_remove(colnames(adj_mat_lst[[1]]), "_[0|1|6]hr$")
covar_table <- covar_table[order(covar_table$cell_line), ]
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

DGELists <- vector("list", length(adj_mat_lst))
voomedObjs <- vector("list", length(adj_mat_lst))
fits <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  DGELists[[i]] <- createDGE(adj_mat_lst[[i]][rownames(adj_mat_lst[[i]]) %in% SZgenes, ])
  voomedObjs[[i]] <- cnfV(DGELists[[i]], designs[[i]])
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
for (i in 1:length(fits_1hr)){
  d <- model.matrix(~0 + disease, data = covar_table[covar_table$time == "1hr", ])
  fit <- lmFit(fits_1hr[[i]], d)
  contr.matrix <- makeContrasts(
    casevscontrol = diseasecase-diseasecontrol,
    levels = colnames(d))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  cvc_1hr[[i]] <- eBayes(fit, trend = F)
  print(summary(decideTests(cvc_1hr[[i]])))
}

# re-fit linear model ####
designs_1v0 <- vector("list", length(fitted_1v0))
designs_6v0 <- vector("list", length(fitted_6v0))
designs_6v1 <- vector("list", length(fitted_6v1))
for (i in 1:length(designs_1v0)) {
  designs_1v0[[i]] <- model.matrix(~0 + disease, 
                                   data = covar_table[covar_table$time == "0hr", ])
  designs_6v0[[i]] <- model.matrix(~0 + disease, 
                                   data = covar_table[covar_table$time == "0hr", ])
  designs_6v1[[i]] <- model.matrix(~0 + disease, 
                                   data = covar_table[covar_table$time == "0hr", ])
}
names(designs_1v0) <- names(fitted_1v0)
names(designs_6v0) <- names(fitted_6v0)
names(designs_6v1) <- names(fitted_6v1)

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
contrasts_6v1 <- vector("list", 3)

for (i in 1:length(fitted_1v0)) {
  contrasts_1v0[[i]] <- contrastOnSubtractedVal(fitted_1v0[[i]], designs_1v0[[i]])
  contrasts_6v0[[i]] <- contrastOnSubtractedVal(fitted_6v0[[i]], designs_6v0[[i]])
  contrasts_6v1[[i]] <- contrastOnSubtractedVal(fitted_6v1[[i]], designs_6v1[[i]])
}


# loop fit ####
types <- names(fitted_1v0)
disease <- factor(covar_table$disease[covar_table$time == "0hr"], levels = c("control", "case"))
age <- covar_table$age2[covar_table$time == "0hr"]
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
