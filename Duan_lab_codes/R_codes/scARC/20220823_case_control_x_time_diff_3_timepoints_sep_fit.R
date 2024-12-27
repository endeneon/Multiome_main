# Chuxuan Li 08/22/2022
# Examine the relationship between logFC and case v control in 5, 18, 20 line 
#combined data using limma, using already generated combated count matrices in 
#17Aug2022_case_control_DE_w_eexps_limma_combat.R


# init ####
library(limma)
library(edgeR)
library(sva)
library(Seurat)

library(ggplot2)
library(stringr)
library(readr)

load("../../../covariates_pooled_together_for_all_linextime_5+18+20.RData")
load("5+18+20line_sepby_type_time_colby_line_adj_mat_list.RData")

# prepare variables ####
lines <- colnames(adj_mat_lst[[1]])
for (i in 1:length(lines)) {
  print(lines[i])
  if (i == 1) {
    batch = unique(covar_table$group[covar_table$cell_line == lines[i]])
  } else {
    batch = c(batch, unique(covar_table$group[covar_table$cell_line == lines[i]]))
  }
}
batch <- as.factor(batch)

covar_table <- covar_table[order(covar_table$cell_line), ]
covar_table_0hr <- covar_table[covar_table$time == "0hr", ]
covar_table_1hr <- covar_table[covar_table$time == "1hr", ]
covar_table_6hr <- covar_table[covar_table$time == "6hr", ]

designs <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  type <- str_split(names(adj_mat_lst)[i], pattern = "-", n = 2, simplify = T)[,1]
  time <- str_split(names(adj_mat_lst)[i], pattern = "-", n = 2, simplify = T)[,2]
  if (time == "0hr") {
    d = covar_table_0hr
  } else if (time == "1hr") {
    d = covar_table_1hr
  } else {
    d = covar_table_6hr
  }
  if (type == "GABA") {
    f = d$GABA_fraction
  } else if (type == "NEFM_neg_glut") {
    f = d$nmglut_fraction
  } else {
    f = d$npglut_fraction
  }
  designs[[i]] <- model.matrix(~0 + disease + age + sex + group + f, 
                               data = d)
}
names(designs) <- names(adj_mat_lst)

# auxillary functions ####
createDGE <- function(count_matrix){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = lines)
  A <- rowSums(y$counts)
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}
filterByCpm <- function(typelists, cutoff, nsample) {
  filters <- vector("list", length(typelists))
  for (i in 1:length(typelists)) {
    print(i)
    cpm <- cpm(typelists[[i]])
    filters[[i]] <- (rowSums(cpm >= cutoff) >= nsample)
  }
  print(length(filters))
  if (length(filters) > 3) {
    stop("list does not contain 3 elements")
  } else {
    passfilter <- filters[[1]] + filters[[2]] + filters[[3]]
    passfilter[passfilter >= 1] <- TRUE
    passfilter[passfilter == 0] <- FALSE
  }
  return(as.logical(passfilter))
}
cnfV <- function(y, design) {
  y <- calcNormFactors(y)
  # voom (no logcpm needed, voom does it within the function)
  v <- voom(y, design, plot = T)
  
  return(v)  
}

# obtain gene list for each cell type across 3 time points ####
GABA_DGEs <- vector("list", length(adj_mat_lst) / 3)
nmglut_DGEs <- vector("list", length(adj_mat_lst) / 3)
npglut_DGEs <- vector("list", length(adj_mat_lst) / 3)
names(adj_mat_lst)
for (i in 1:length(adj_mat_lst)) {
  ind <- i %% 3
  if (ind == 0) {
    ind <- 3
  }
  if (str_detect(names(adj_mat_lst)[i], "GABA")) {
    GABA_DGEs[[ind]] <- createDGE(adj_mat_lst[[i]])
    names(GABA_DGEs)[ind] <- names(adj_mat_lst)[i]
  } else if (str_detect(names(adj_mat_lst)[i], "NEFM_neg_glut")) {
    nmglut_DGEs[[ind]] <- createDGE(adj_mat_lst[[i]])
    names(nmglut_DGEs)[ind] <- names(adj_mat_lst)[i]
  } else {
    npglut_DGEs[[ind]] <- createDGE(adj_mat_lst[[i]])
    names(npglut_DGEs)[ind] <- names(adj_mat_lst)[i]
  }
}  

GABA_filter <- filterByCpm(GABA_DGEs, 1, 21)
sum(GABA_filter) #17712
nmglut_filter <- filterByCpm(nmglut_DGEs, 1, 21)
sum(nmglut_filter) #17022
npglut_filter <- filterByCpm(npglut_DGEs, 1, 21)
sum(npglut_filter) #17566

# limma voom ####
voomedDGEs <- vector("list", length(adj_mat_lst))
fits <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  ind <- i %% 3
  if (ind == 0) {
    ind <- 3
  }
  print(paste(names(adj_mat_lst)[i], ind, i))
  if (str_detect(names(adj_mat_lst)[i], "GABA")) {
    voomedDGEs[[i]] <- cnfV(GABA_DGEs[[ind]][GABA_filter, ], designs[[i]])
    fit <- lmFit(voomedDGEs[[i]], designs[[i]])
    fits[[i]] <- eBayes(fit)
  } else if (str_detect(names(adj_mat_lst)[i], "NEFM_neg_glut")) {
    voomedDGEs[[i]] <- cnfV(nmglut_DGEs[[ind]][nmglut_filter, ], designs[[i]])
    fit <- lmFit(voomedDGEs[[i]], designs[[i]])
    fits[[i]] <- eBayes(fit)
  } else {
    voomedDGEs[[i]] <- cnfV(npglut_DGEs[[ind]][npglut_filter, ], designs[[i]])
    fit <- lmFit(voomedDGEs[[i]], designs[[i]])
    fits[[i]] <- eBayes(fit)
  }
}
names(fits) <- names(adj_mat_lst)

# extract fitted value differences ####
fitted_1v0 <- vector("list", 3)
fitted_6v0 <- vector("list", 3)
for (i in 1:length(fitted_1v0)) {
  ind_0hr <- 3 * i - 2
  ind_1hr <- 3 * i - 1
  ind_6hr <- 3 * i
  fitted_1v0[[i]] <- fitted.MArrayLM(fits[[ind_1hr]]) - fitted.MArrayLM(fits[[ind_0hr]])
  names(fitted_1v0)[i] <- names(fits)[ind_1hr]
  fitted_6v0[[i]] <- fitted.MArrayLM(fits[[ind_6hr]]) - fitted.MArrayLM(fits[[ind_0hr]])
  names(fitted_6v0)[i] <- names(fits)[ind_6hr]
}

# re-fit linear model ####
designs_1v0 <- vector("list", length(fitted_1v0))
designs_6v0 <- vector("list", length(fitted_6v0))

for (i in 1:length(designs_1v0)) {
  type <- str_split(names(fitted_1v0)[i], pattern = "-", n = 2, simplify = T)[,1]
  if (type == "GABA") {
    fraction_0hr = covar_table_0hr$GABA_fraction
    fraction_1hr = covar_table_1hr$GABA_fraction
  } else if (type == "NEFM_neg_glut") {
    fraction_0hr = covar_table_0hr$nmglut_fraction
    fraction_1hr = covar_table_1hr$nmglut_fraction
  } else {
    fraction_0hr = covar_table_0hr$npglut_fraction
    fraction_1hr = covar_table_1hr$npglut_fraction
  }
  designs_1v0[[i]] <- model.matrix(~0 + disease + age + sex + group + fraction_0hr + fraction_1hr, 
                               data = d)
}
for (i in 1:length(designs_6v0)) {
  type <- str_split(names(fitted_6v0)[i], pattern = "-", n = 2, simplify = T)[,1]
  if (type == "GABA") {
    fraction_0hr = covar_table_0hr$GABA_fraction
    fraction_6hr = covar_table_6hr$GABA_fraction
  } else if (type == "NEFM_neg_glut") {
    fraction_0hr = covar_table_0hr$nmglut_fraction
    fraction_6hr = covar_table_6hr$nmglut_fraction
  } else {
    fraction_0hr = covar_table_0hr$npglut_fraction
    fraction_6hr = covar_table_6hr$npglut_fraction
  }
  designs_6v0[[i]] <- model.matrix(~0 + disease + age + sex + group + fraction_0hr + fraction_6hr, 
                               data = d)
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

