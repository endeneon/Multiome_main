# Chuxuan Li 08/26/2022
# Use fitted values obtained from 24Aug2022_case_control_x_time_fit_3timepoints_together.R
# Redo fit linear model using lm() instead of limma, then extract the beta
#for disease, the p-value for beta, and adjust p-values using BH; then manually
#compute fold change and do t-test on the FC.

# init ####
library(limma)
library(edgeR)
library(ggplot2)
library(stringr)
library(readxl)

load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/5+18+20line_residuals_matrices_w_agesqd_notime.RData")
load("../../../../covariates_pooled_together_for_all_linextime_5+18+20.RData")
covar_table <- covar_table[order(covar_table$cell_line), ]

# loop fit ####
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
  betas_6v1[[i]]$p.adj <- p.adjust(betas_6v1[[i]]$, method = "BH")
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
    paste0("./case_v_control_DE/lmres/positive_beta_filtered_by_padj/",
          names(betas_6v1)[i], "_differential_6v1response_genes_positiveb_sig_by_padj.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_6v1[[i]][betas_6v1[[i]]$p.adj < 0.05 & betas_6v1[[i]]$beta < 0, ]
  filename <- 
    paste0("./case_v_control_DE/lmres/negative_beta_filtered_by_padj/",
           names(betas_6v1)[i], "_differential_6v1response_genes_negativeb_sig_by_padj.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_6v1[[i]][betas_6v1[[i]]$pvalue < 0.05 & betas_6v1[[i]]$beta > 0, ]
  filename <- 
    paste0("./case_v_control_DE/lmres/positive_beta_filtered_by_pvalue/",
           names(betas_6v1)[i], "_differential_6v1response_genes_positiveb_sig_by_pvalue.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_6v1[[i]][betas_6v1[[i]]$pvalue < 0.05 & betas_6v1[[i]]$beta < 0, ]
  filename <- 
    paste0("./case_v_control_DE/lmres/negative_beta_filtered_by_pvalue/",
           names(betas_6v1)[i], "_differential_6v1response_genes_negativeb_sig_by_pvalue.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_6v1[[i]]
  filename <-
    paste0("./case_v_control_DE/lmres/unfiltered_full_results/",
           names(betas_6v1)[i], "_differential_6v1response_genes_full.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
} 

# summary table ####
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
  deg_counts1[i,1] <- sum(betas_1v0[[i]]$p.adj < 0.05 & betas_1v0[[i]]$beta > 0)
  deg_counts1[i,2] <- sum(betas_1v0[[i]]$p.adj < 0.05 & betas_1v0[[i]]$beta < 0)
  deg_counts1[i,3] <- sum(betas_1v0[[i]]$p.adj > 0.05)
  deg_counts1[i,4] <- nrow(betas_1v0[[i]])
  deg_counts6[i,1] <- sum(betas_6v0[[i]]$p.adj < 0.05 & betas_6v0[[i]]$beta > 0)
  deg_counts6[i,2] <- sum(betas_6v0[[i]]$p.adj < 0.05 & betas_6v0[[i]]$beta < 0)
  deg_counts6[i,3] <- sum(betas_6v0[[i]]$p.adj > 0.05)
  deg_counts6[i,4] <- nrow(betas_6v0[[i]])
  deg_counts61[i,1] <- sum(betas_6v1[[i]]$p.adj < 0.05 & betas_6v1[[i]]$beta > 0)
  deg_counts61[i,2] <- sum(betas_6v1[[i]]$p.adj < 0.05 & betas_6v1[[i]]$beta < 0)
  deg_counts61[i,3] <- sum(betas_6v1[[i]]$p.adj > 0.05)
  deg_counts61[i,4] <- nrow(betas_6v1[[i]])
} 
deg_counts <- rbind(deg_counts1, deg_counts6, deg_counts61)
write.table(deg_counts, 
            file = "./case_v_control_DE/lmres/casevcontrol_lm_beta_filtered_by_padj_summary.csv", 
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
  pdf(paste0("./case_v_control_DE/lmres/volcano_plots/by_pvalue/", 
             names(betas_1v0)[i], "_1v0hr_volcano_plot.pdf"))
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

# dotplot ####
hist(betas_1v0$NEFM_neg_glut$Estimate)
max(betas_1v0$NEFM_neg_glut$Estimate)
min(betas_1v0$NEFM_neg_glut$Estimate)
betas_combined <- c(betas_1v0, betas_6v0)
names(betas_combined)[1:3] <- paste(names(betas_combined)[1:3], "1v0", sep = "-")
names(betas_combined)[4:6] <- paste(names(betas_combined)[4:6], "6v0", sep = "-")

deg_df <- NULL
for (i in 1:length(betas_1v0)) {
  print(names(betas_1v0)[i])
  df <- betas_1v0[[i]][betas_1v0[[i]]$p.adj < 0.05, ]
  if (nrow(df) != 0) {
    print("nonzero")
    df$typextime <- names(betas_1v0)[i]
    if (is.null(deg_df)) {
      deg_df <- df
    } else (
      deg_df <- rbind(deg_df, df)
    )
  } 
}
deg_df$neg_log_p <- (-1) * log2(deg_df$p.adj)
deg_df$genes <- rownames(deg_df)
ggplot(deg_df, aes(x = genes, y = typextime, color = Estimate, size = neg_log_p)) +
  geom_point() +
  scale_color_gradient2(low = "steelblue", high = "darkred") +
  ylab("Cell type - Time point") +
  labs(size = "-log(adjusted p-val)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

deg_df <- NULL
for (i in 1:length(betas_6v0)) {
  print(names(betas_6v0)[i])
  df <- betas_6v0[[i]][betas_6v0[[i]]$p.adj < 0.05, ]
  if (nrow(df) != 0) {
    print("nonzero")
    df$typextime <- names(betas_6v0)[i]
    df$genes <- rownames(df)
    if (is.null(deg_df)) {
      deg_df <- df
    } else {
      if (i == 3) {
        df <- df[df$genes %in% deg_df$genes, ]
      }
      deg_df <- rbind(deg_df, df)
    }
  } 
}
deg_df$neg_log_p <- (-1) * log2(deg_df$p.adj)
ggplot(deg_df, aes(x = genes, y = typextime, color = Estimate, size = neg_log_p)) +
  geom_point() +
  scale_color_gradient2(low = "steelblue", high = "darkred") +
  ylab("Cell type - Time point") +
  labs(size = "-log(adjusted p-val)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# calculate fold change and test ####
ttestfeeder <- function(row) {
  ttres = t.test(row$fitted.values[disease == "case"], row$fitted.values[disease == "control"])
  pval = ttres$p.value
  diff = sum(row$fitted.values[disease == "case"]/sum(disease == "case")) - 
    sum(row$fitted.values[disease == "control"]/sum(disease == "control"))
  return(c(diff, pval))
}
fc1 <- vector("list", 3)
for (i in 1:length(fits_1v0)) {
  fc1[[i]] <-lapply(fits_1v0[[i]], ttestfeeder)
}
fc6 <- vector("list", 3)
for (i in 1:length(fits_1v0)) {
  fc6[[i]] <-lapply(fits_6v0[[i]], ttestfeeder)
}
save(fc1, fc6, betas_1v0, betas_6v0, file = "FC_beta_results_0826.RData")

