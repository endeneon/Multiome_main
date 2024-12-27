# Chuxuan Li 09/15/2022
# Repeat case v control using residual values from first pass linear model fitting,
#this time do t-test directly, comparing the average residuals from case/control

# init ####
library(limma)
library(edgeR)
library(ggplot2)
library(stringr)
library(readr)
library(readxl)

load("../../../covariates_pooled_together_for_all_linextime_5+18+20.RData")
load("5+18+20line_residuals_matrices_w_agesqd_notime.RData")
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/case_v_control_DE/5+18+20line_sepby_type_colby_linxetime_adj_mat_list.RData")

# auxillary functions ####
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

celllines <- str_remove(colnames(adj_mat_lst[[1]]), "_[0|1|6]hr$")
covar_table <- covar_table[order(covar_table$cell_line), ]
covar_table$age2 <- covar_table$age * covar_table$age

# do t-test on the residuals
fitted_list <- list(fitted_1v0, fitted_6v0, fitted_6v1)
case_lines <- unique(covar_table$cell_line[covar_table$disease == "case"])
ctrl_lines <- unique(covar_table$cell_line[covar_table$disease == "control"])


tres_1v0 <- vector("list", length(fitted_1v0))
tres_6v0 <- vector("list", length(fitted_6v0))
tres_6v1 <- vector("list", length(fitted_6v1))
tres <- list(tres_1v0, tres_6v0, tres_6v1)
for (i in 1:length(fitted_list)) {
  for (j in 1:length(fitted_list[[i]])) {
    case <- fitted_list[[i]][[j]][, str_extract(colnames(fitted_list[[i]][[j]]), "CD_[0-9][0-9]") %in% case_lines]
    ctrl <- fitted_list[[i]][[j]][, str_extract(colnames(fitted_list[[i]][[j]]), "CD_[0-9][0-9]") %in% ctrl_lines]
    tres[[i]][[j]] <- data.frame(matrix(ncol = 3, nrow = nrow(case), 
                                        dimnames = list(rownames(case),
                                                        c("mean_x", "mean_y", "p.value"))))
    for (k in 1:nrow(case)){
      row1 <- case[k, ]
      row2 <- ctrl[rownames(ctrl) == rownames(case)[k], ]
      res <- t.test(row1, row2, paired = F)
      tres[[i]][[j]]$mean_x[k] <- res$estimate[1]
      tres[[i]][[j]]$mean_y[k] <- res$estimate[2]
      tres[[i]][[j]]$p.value[k] <- res$p.value
    }
    tres[[i]][[j]]$p.adj <- p.adjust(tres[[i]][[j]]$p.value, method = "BY")
    print(sum(tres[[i]][[j]]$p.value < 0.05))
  }
}
