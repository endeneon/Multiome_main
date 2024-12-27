# Siwei 30 Jan 2024
# test differential response by sex

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
createDGE <- 
  function(count_matrix) {
  y <- 
    DGEList(counts = count_matrix, 
            genes = rownames(count_matrix), 
            group = celllines)
  A <- rowSums(y$counts)
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}

filterByCpm <- 
  function(df1, df2, cutoff, proportion = 1) {
  cpm1 <- cpm(df1)
  cpm2 <- cpm(df2)
  passfilter <- (rowSums(cpm1 >= cutoff) >= ncol(cpm1) * proportion |
                   rowSums(cpm2 >= cutoff) >= ncol(cpm2) * proportion) # number of samples from either group > ns
  return(passfilter)
}

cnfV <- 
  function(y, design) {
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = T)
  return(v)  
}

# make design matrices ####
covar_table_final <- 
  covar_table_final[order(covar_table_final$cell_line), ]
covar_table_final$age_sq <- 
  covar_table_final$age * covar_table_final$age
covar_table_0hr <- 
  covar_table_final[covar_table_final$time == "0hr", ]
covar_table_1hr <- 
  covar_table_final[covar_table_final$time == "1hr", ]
covar_table_6hr <- 
  covar_table_final[covar_table_final$time == "6hr", ]
tables <- 
  list(covar_table_0hr, 
       covar_table_1hr, 
       covar_table_6hr)
names(tables) <- 
  c("covar_table_0hr", 
    "covar_table_1hr", 
    "covar_table_6hr")

designs_0hr <- vector("list", length(adj_mat_lst))
designs_1hr <- vector("list", length(adj_mat_lst))
designs_6hr <- vector("list", length(adj_mat_lst))
designs <- list(designs_0hr, 
                designs_1hr, 
                designs_6hr)
names(designs) <- c("designs_0hr", 
                    "designs_1hr", 
                    "designs_6hr")

celllines <- 
  unique(str_remove(colnames(adj_mat_lst[[1]]), 
                    pattern = "_[0|1|6]hr$"))

# Here keep the $sex factor, regress out the $aff
for (i in 1:length(tables)) {
  table <- tables[[i]]
  for (j in 1:length(adj_mat_lst)) {
    type <- names(adj_mat_lst)[j]
    print(paste(i, j))
    if (type == "GABA") {
      f = table$GABA_fraction
      # print(f)
    } else if (type == "NEFM_neg_glut") {
      f = table$nmglut_fraction
    } else {
      f = table$npglut_fraction
    }
    designs[[i]][[j]] <- model.matrix(~ 0 + 
                                        batch + 
                                        age_sq + 
                                        aff + 
                                        f, 
                                 data = table)
  }
  names(designs[[i]]) <- names(adj_mat_lst)
}

# split count matrices ####
df_0hr <- list(adj_mat_lst[[1]][, str_detect(colnames(adj_mat_lst[[1]]), "0hr")],
               adj_mat_lst[[2]][, str_detect(colnames(adj_mat_lst[[2]]), "0hr")],
               adj_mat_lst[[3]][, str_detect(colnames(adj_mat_lst[[3]]), "0hr")])
df_1hr <- list(adj_mat_lst[[1]][, str_detect(colnames(adj_mat_lst[[1]]), "1hr")],
               adj_mat_lst[[2]][, str_detect(colnames(adj_mat_lst[[2]]), "1hr")],
               adj_mat_lst[[3]][, str_detect(colnames(adj_mat_lst[[3]]), "1hr")])
df_6hr <- list(adj_mat_lst[[1]][, str_detect(colnames(adj_mat_lst[[1]]), "6hr")],
               adj_mat_lst[[2]][, str_detect(colnames(adj_mat_lst[[2]]), "6hr")],
               adj_mat_lst[[3]][, str_detect(colnames(adj_mat_lst[[3]]), "6hr")])
df_list <- list(df_0hr, 
                df_1hr, 
                df_6hr)
names(df_list) <- 
  c("df_0hr", 
    "df_1hr", 
    "df_6hr")

# filter genes by cpm in either 1 of 2 groups ####
voomed_1v0 <- vector("list", length(adj_mat_lst))
voomed_6v0 <- vector("list", length(adj_mat_lst))
fitted_1v0 <- vector("list", length(adj_mat_lst))
fitted_6v0 <- vector("list", length(adj_mat_lst))

for (i in 1:length(df_list[[1]])) {
  passfilter_1v0 <- filterByCpm(df_list[[1]][[i]], df_list[[2]][[i]], 1)
  passfilter_6v0 <- filterByCpm(df_list[[1]][[i]], df_list[[3]][[i]], 1)
  
  voomed_1v0[[i]] <- list(cnfV(createDGE(df_list[[1]][[i]][passfilter_1v0, ]), 
                                          designs[[1]][[i]]),
                          cnfV(createDGE(df_list[[2]][[i]][passfilter_1v0, ]), 
                               designs[[2]][[i]]))
  voomed_6v0[[i]] <- list(cnfV(createDGE(df_list[[1]][[i]][passfilter_6v0, ]), 
                               designs[[1]][[i]]),
                          cnfV(createDGE(df_list[[3]][[i]][passfilter_6v0, ]), 
                               designs[[3]][[i]]))
  
  fitted_1v0[[i]] <- list(lmFit(voomed_1v0[[i]][[1]], designs[[1]][[i]]),
                          lmFit(voomed_1v0[[i]][[2]], designs[[2]][[i]]))
  fitted_6v0[[i]] <- list(lmFit(voomed_6v0[[i]][[1]], designs[[1]][[i]]),
                          lmFit(voomed_6v0[[i]][[2]], designs[[3]][[i]]))
}
names(fitted_1v0) <- names(adj_mat_lst)
names(fitted_6v0) <- names(adj_mat_lst)


# extract residuals ####
resid_1v0 <- vector("list", length(adj_mat_lst))
resid_6v0 <- vector("list", length(adj_mat_lst))
for (i in 1:length(fitted_1v0)) {
  resid_1v0[[i]] <- 
    list(residuals.MArrayLM(object = fitted_1v0[[i]][[1]], 
                            voomed_1v0[[i]][[1]]),
         residuals.MArrayLM(object = fitted_1v0[[i]][[2]], 
                            voomed_1v0[[i]][[2]]))
  resid_6v0[[i]] <- 
    list(residuals.MArrayLM(object = fitted_6v0[[i]][[1]], 
                            voomed_6v0[[i]][[1]]),
         residuals.MArrayLM(object = fitted_6v0[[i]][[2]], 
                            voomed_6v0[[i]][[2]]))
}
names(resid_1v0) <- names(adj_mat_lst)
names(resid_1v0) <- names(adj_mat_lst)

resid_diff_1v0 <- vector("list", length(adj_mat_lst))
resid_diff_6v0 <- vector("list", length(adj_mat_lst))

for (i in 1:length(resid_1v0)) {
  resid_diff_1v0[[i]] <- resid_1v0[[i]][[1]] - resid_1v0[[i]][[2]]
  resid_diff_6v0[[i]] <- resid_6v0[[i]][[1]] - resid_6v0[[i]][[2]]
}
names(resid_diff_1v0) <- names(adj_mat_lst)
names(resid_diff_6v0) <- names(adj_mat_lst)

save(resid_diff_1v0, resid_diff_6v0, file = "cpm_filtered_residual_differences_1v0_6v0_dif_sex.RData")

# re-fit linear model exp = b0 + b1 * aff with limma ####
designs_1v0 <- vector("list", length(resid_diff_1v0))
designs_6v0 <- vector("list", length(resid_diff_6v0))
for (i in 1:length(designs_1v0)) {
  designs_1v0[[i]] <- 
    model.matrix(~ 0 + sex, 
                 data = covar_table_1hr)
  designs_6v0[[i]] <- 
    model.matrix(~ 0 + sex, 
                 data = covar_table_0hr)
}
names(designs_1v0) <- names(resid_diff_1v0)
names(designs_6v0) <- names(resid_diff_6v0)

# ! contrast was made by M vs F !! ####
contrastOnSubtractedVal <- 
  function(df, design) {
  fit <- lmFit(df, design)
  contr.matrix <- makeContrasts(
    MvsF = sexM - sexF,
    levels = colnames(design))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit <- eBayes(fit, trend = T)
  print(summary(decideTests(fit)))
  return(fit)
}
contrasts_1v0 <- vector("list", 3)
contrasts_6v0 <- vector("list", 3)

for (i in 1:length(resid_diff_1v0)) {
  contrasts_1v0[[i]] <- contrastOnSubtractedVal(resid_diff_1v0[[i]], designs_1v0[[i]])
  contrasts_6v0[[i]] <- contrastOnSubtractedVal(resid_diff_6v0[[i]], designs_6v0[[i]])
}


# loop fit linear model exp = b0 + b1 * sex with lm ####
types <- names(resid_diff_1v0)
## since M vs F, F should be the first level !
factor_sex <- 
  factor(covar_table_final$sex[covar_table_final$time == "0hr"], 
         levels = c("F", "M"))
times <- sort(unique(covar_table_final$time))

lmfeeder <- function(row) {
  return(lm(row ~ factor_sex))
}

cvc_1v0 <- vector("list", 3)
cvc_6v0 <- vector("list", 3)
betas_1v0 <- vector("list", 3)
betas_6v0 <- vector("list", 3)

for (i in 1:length(resid_1v0)) {
  cvc_1v0[[i]] <- apply(X = resid_diff_1v0[[i]], 
                        MARGIN = 1, 
                        FUN = lmfeeder)
  cvc_6v0[[i]] <- apply(X = resid_diff_6v0[[i]], 
                        MARGIN = 1, 
                        FUN = lmfeeder)
}

for (i in 1:length(cvc_1v0)) {
  for (j in 1:length(cvc_1v0[[i]])) {
    if (j == 1) {
      mat1 <- summary(cvc_1v0[[i]][[j]])$coefficients["factor_sexM", ]
    } else {
      mat1 <- rbind(mat1, summary(cvc_1v0[[i]][[j]])$coefficients["factor_sexM", ])
    }
  }
  rownames(mat1) <- names(cvc_1v0[[i]])
  betas_1v0[[i]] <- as.data.frame(mat1)
}

for (i in 1:length(cvc_6v0)) {
  for (j in 1:length(cvc_6v0[[i]])) {
    if (j == 1) {
      mat1 <- summary(cvc_6v0[[i]][[j]])$coefficients["factor_sexM", ]
    } else {
      mat1 <- rbind(mat1, summary(cvc_6v0[[i]][[j]])$coefficients["factor_sexM", ])
    }
  }
  rownames(mat1) <- names(cvc_6v0[[i]])
  betas_6v0[[i]] <- as.data.frame(mat1)
}
names(betas_1v0) <- names(resid_diff_1v0)
names(betas_6v0) <- names(resid_diff_6v0)

# adjust p
for (i in 1:length(betas_1v0)) {
  betas_1v0[[i]]$p.adj <- p.adjust(betas_1v0[[i]]$`Pr(>|t|)`, method = "BH")
  betas_6v0[[i]]$p.adj <- p.adjust(betas_6v0[[i]]$`Pr(>|t|)`, method = "BH")
}


# output results ####
for (i in 1:length(betas_1v0)) {
  colnames(betas_1v0[[i]]) <- c("beta", "SEM", "T", "pvalue", "p.adj")
  colnames(betas_6v0[[i]]) <- c("beta", "SEM", "T", "pvalue", "p.adj")
}

# summary table 
deg_counts1 <- array(dim = c(3, 4), 
                     dimnames = list(paste(types, 
                                           "1v0", 
                                           sep = "-"),
                                     c("positive", 
                                       "negative", 
                                       "nonsig", 
                                       "total"))) 
for (i in 1:length(betas_1v0)) {
  deg_counts1[i, 1] <- sum(betas_1v0[[i]]$pvalue < 0.05 & betas_1v0[[i]]$beta > 0)
  deg_counts1[i, 2] <- sum(betas_1v0[[i]]$pvalue < 0.05 & betas_1v0[[i]]$beta < 0)
  deg_counts1[i, 3] <- sum(betas_1v0[[i]]$pvalue > 0.05)
  deg_counts1[i, 4] <- nrow(betas_1v0[[i]])
} 

deg_counts6 <- array(dim = c(3, 4), 
                     dimnames = list(paste(types, 
                                           "6v0", 
                                           sep = "-"),
                                     c("positive", 
                                       "negative", 
                                       "nonsig", 
                                       "total"))) 
for (i in 1:length(betas_6v0)) {
  deg_counts6[i,1] <- sum(betas_6v0[[i]]$pvalue < 0.05 & betas_6v0[[i]]$beta > 0)
  deg_counts6[i,2] <- sum(betas_6v0[[i]]$pvalue < 0.05 & betas_6v0[[i]]$beta < 0)
  deg_counts6[i,3] <- sum(betas_6v0[[i]]$pvalue > 0.05)
  deg_counts6[i,4] <- nrow(betas_6v0[[i]])
} 

output_dir <- "Siwei_MvsF_lm_results_30Jan2024"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

deg_counts <- rbind(deg_counts1, 
                    deg_counts6)

write.table(deg_counts, 
            file = paste(".",
                         output_dir,
                         "Siwei_MvsF_lm_018-029_30Jan2024.txt",
                         sep = "/"), 
            sep = "\t", quote = F)

for (i in 1:length(betas_1v0)) {
  print(i)
  df <- betas_1v0[[i]]
  df <- df[order(df$pvalue), ]
  filename <-
    paste0("./",
           output_dir,
           "/",
           names(betas_1v0)[i], "_1v0response_all_genes_MvsF.txt")
  write.table(df, 
              file = filename, 
              quote = F, sep = "\t", 
              row.names = T, 
              col.names = T)
  
  df <- betas_6v0[[i]]
  df <- df[order(df$pvalue), ]
  filename <-
    paste0("./",
           output_dir,
           "/",
           names(betas_6v0)[i], "_6v0response_all_genes_MvsF.txt")
  write.table(df, 
              file = filename, 
              quote = F, sep = "\t", 
              row.names = T, 
              col.names = T)
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
for (i in 1:length(betas_6v0)){
  res <- betas_6v0[[i]]
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
  png(paste0("./018-029_diff_response_filter_by_cpm_casevcontrol_lm_results/volcano_plots/", 
             names(betas_6v0)[i], "_casevcontrol_in_6v0_response_volcano_plot.png"), 
      width = 750, height = 750)
  p <- ggplot(data = as.data.frame(res),
              aes(x = beta, 
                  y = neg_log_pval, 
                  color = significance,
                  label = labelling)) + 
    geom_point(size = 0.6) +
    xlab("Beta") +
    ylab("-log10(p-value)") +
    scale_color_manual(values = dotColors) +
    theme_minimal(base_size = 15) +
    geom_text_repel(box.padding = unit(0.05, 'lines'),
                    min.segment.length = 0,
                    force = 0.2,
                    max.overlaps = Inf,
                    force_pull = 0.1,
                    show.legend = F) +
    xlim(c(-1, 1)) +
    ggtitle(paste0(str_replace(str_replace(str_replace_all(names(betas_6v0)[i], 
                                                           "_", " "), " pos", "+"), 
                               " neg", "-")), "Case vs control in 6v0 response genes")
  print(p)
  dev.off()
  pdf(paste0("./018-029_diff_response_filter_by_cpm_casevcontrol_lm_results/volcano_plots/", 
             names(betas_6v0)[i], "_casevcontrol_in_6v0_response_volcano_plot.pdf"))
  print(p)
  dev.off()
}
