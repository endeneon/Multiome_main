# Chuxuan Li 02/06/2023
# differential response analysis on 018-029 data, filtering genes with previous
#limma DE results

# init ####
library(limma)
library(edgeR)
library(ggplot2)
library(stringr)
library(readr)
library(readxl)
load("./018-029_covar_table_final.RData")
load("./residual_differences_1v0_6v0_6v1.RData")

paths_024 <- list.files("../Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/pseudobulk_DE/res_use_combatseq_mat/limma_voom_res/adj_cpm_9in18_res/",
          pattern = "[1|6]v0.*regulated_DEGs.csv", recursive = T, full.names = T, include.dirs = F)
paths_024 <- sort(paths_024[str_detect(paths_024, "unfiltered", T)])

typextime <- str_extract(paths_024, "[A-Za-z]+_res_[1|6]v0_[a-z]+")
path_lst <- vector("list", length(typextime))
for (i in 1:length(typextime)) {
  path_lst[[i]] <- paths_024[str_detect(paths_024, typextime[i])]
}
names(path_lst) <- typextime

paths_025 <- list.files("../Duan_project_025_RNA/Analysis_part2_GRCh38/limma_DE_results/",
                        pattern = "[1|6]v0.*regulated_DEGs.csv", recursive = T, 
                        full.names = T, include.dirs = F)
paths_025 <- paths_025[str_detect(paths_025, "unfiltered", T)]
for (i in 1:length(typextime)) {
  path_lst[[i]] <- c(path_lst[[i]],
                     paths_025[str_detect(paths_025, typextime[i])])
}

paths_022 <- list.files("../R_scARC_Duan_018/Analysis_RNAseq_pseudobulk/20line_limma/",
                        pattern = "[1|6]v0.*regulated_DEGs.csv", recursive = T, 
                        full.names = T, include.dirs = F)
paths_022 <- paths_022[str_detect(paths_022, "unfiltered", T)]
for (i in 1:length(typextime)) {
  path_lst[[i]] <- c(path_lst[[i]],
                     paths_022[str_detect(paths_022, typextime[i])])
}

paths_029 <- list.files("../Duan_project_029_RNA/limma_DE_results",
                        pattern = "[1|6]v0.*regulated_DEGs.csv", recursive = T, 
                        full.names = T, include.dirs = F)
paths_029 <- paths_029[str_detect(paths_029, "unfiltered", T)]
for (i in 1:length(typextime)) {
  path_lst[[i]] <- c(path_lst[[i]],
                     paths_029[str_detect(paths_029, typextime[i])])
}

# intersect the upregulated/downregulated genes from each batch
genes_lst <- vector("list", length(path_lst))
names(genes_lst) <- typextime
for (i in 1:length(path_lst)) {
  for (j in 1:length(path_lst[[i]])) {
    DEGs <- read_csv(path_lst[[i]][j])
    DEGs <- DEGs$genes
    if (j == 1) {
      genes_lst[[i]] <- DEGs
    } else {
      genes_lst[[i]] <- intersect(DEGs, genes_lst[[i]])
    }
  }
}
# take union of the up/downregulated genes to get final gene lists
typextime <- unique(str_extract(paths_024, "[A-Za-z]+_res_[1|6]v0"))
final_genelist <- vector("list", length(typextime))
for (i in 1:length(typextime)) {
  final_genelist[[i]] <- genes_lst[str_detect(names(genes_lst), typextime[i])]
}
names(final_genelist) <- typextime
final_genelist_1v0 <- final_genelist[str_detect(names(final_genelist), "1v0")]
final_genelist_6v0 <- final_genelist[str_detect(names(final_genelist), "6v0")]

# loop fit linear model exp = b0 + b1 * aff with lm ####
types <- names(resid_1v0)
aff <- factor(covar_table_final$aff[covar_table_final$time == "0hr"], levels = c("control", "case"))
times <- sort(unique(covar_table_final$time))

lmfeeder <- function(row) {
  return(lm(row ~ aff))
}
cvc_1v0 <- vector("list", 3)
cvc_6v0 <- vector("list", 3)
betas_1v0 <- vector("list", 3)
betas_6v0 <- vector("list", 3)
for (i in 1:length(resid_1v0)) {
  cvc_1v0[[i]] <- apply(X = resid_1v0[[i]][rownames(resid_1v0[[i]]) %in% final_genelist_1v0[[i]], ], 
                        MARGIN = 1, FUN = lmfeeder)
  cvc_6v0[[i]] <- apply(X = resid_6v0[[i]][rownames(resid_6v0[[i]]) %in% final_genelist_6v0[[i]], ], 
                        MARGIN = 1, FUN = lmfeeder)
}
for (i in 1:length(cvc_1v0)) {
  for(j in 1:length(cvc_1v0[[i]])) {
    if (j == 1) {
      mat <- summary(cvc_1v0[[i]][[j]])$coefficients["affcase", ]
    } else {
      mat <- rbind(mat, summary(cvc_1v0[[i]][[j]])$coefficients["affcase", ])
    }
  }
  rownames(mat) <- names(cvc_1v0[[i]])
  betas_1v0[[i]] <- as.data.frame(mat)
}
for (i in 1:length(cvc_6v0)) {
  for(j in 1:length(cvc_6v0[[i]])) {
    if (j == 1) {
      mat <- summary(cvc_6v0[[i]][[j]])$coefficients["affcase", ]
    } else {
      mat <- rbind(mat, summary(cvc_6v0[[i]][[j]])$coefficients["affcase", ])
    }
  }
  rownames(mat) <- names(cvc_6v0[[i]])
  betas_6v0[[i]] <- as.data.frame(mat)
}
names(betas_1v0) <- names(resid_1v0)
names(betas_6v0) <- names(resid_6v0)

# adjust p
for (i in 1:length(betas_1v0)) {
  betas_1v0[[i]]$p.adj <- p.adjust(betas_1v0[[i]]$`Pr(>|t|)`, method = "BH")
  betas_6v0[[i]]$p.adj <- p.adjust(betas_6v0[[i]]$`Pr(>|t|)`, method = "BH")
}


# output results ####
for (i in 1:length(betas_1v0)) {
  colnames(betas_1v0[[i]]) <- c("beta", "Std. Error", "t value", "pvalue", "p.adj")
  colnames(betas_6v0[[i]]) <- c("beta", "Std. Error", "t value", "pvalue", "p.adj")
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

for (i in 1:length(betas_6v0)) {
  df <- betas_6v0[[i]][betas_6v0[[i]]$p.adj < 0.05 & betas_6v0[[i]]$beta > 0, ]
  filename <- 
    paste0("./018-029_response_as_residual_difference_casevcontrol_lm_results/filtered_by_padj/",
           names(betas_6v0)[i], "_6v0response_upregulated_padj.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_6v0[[i]][betas_6v0[[i]]$p.adj < 0.05 & betas_6v0[[i]]$beta < 0, ]
  filename <- 
    paste0("./018-029_response_as_residual_difference_casevcontrol_lm_results/filtered_by_padj/",
           names(betas_6v0)[i], "_6v0response_downregulated_padj.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_6v0[[i]][betas_6v0[[i]]$pvalue < 0.05 & betas_6v0[[i]]$beta > 0, ]
  filename <- 
    paste0("./018-029_response_as_residual_difference_casevcontrol_lm_results/filtered_by_pvalue/",
           names(betas_6v0)[i], "_6v0response_upregulated_pvalue.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_6v0[[i]][betas_6v0[[i]]$pvalue < 0.05 & betas_6v0[[i]]$beta < 0, ]
  filename <- 
    paste0("./018-029_response_as_residual_difference_casevcontrol_lm_results/filtered_by_pvalue/",
           names(betas_6v0)[i], "_6v0response_downregulated_pvalue.csv")
  write.table(df, file = filename, quote = F, sep = ",", row.names = T, col.names = T)
  df <- betas_6v0[[i]]
  filename <-
    paste0("./018-029_response_as_residual_difference_casevcontrol_lm_results/unfiltered_full_results/",
           names(betas_6v0)[i], "_6v0response_all_genes.csv")
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
