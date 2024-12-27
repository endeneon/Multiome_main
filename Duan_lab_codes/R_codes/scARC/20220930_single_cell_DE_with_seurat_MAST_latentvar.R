# Chuxuan Li 09/30/2022
# compute single-cell based differentially expressed genes with findmarkers()
#and MAST, plus co-cultured batch as latent.vars

# init ####
library(Seurat)
library(SeuratWrappers)

library(dplyr)
library(tidyr)
library(stringr)

library(ggplot2)
library(RColorBrewer)

load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")
load("../covariates_pooled_together_for_all_linextime_5+18+20.RData")

# assign covariates based on cell line ident ####
obj <- subset(integrated_labeled, cell.line.ident != "unmatched")
obj$cell.type[obj$cell.type %in% c("GABA", "SST_pos_GABA", "SEMA3E_pos_GABA")] <- "GABA"
rm(integrated_labeled)

unique(obj$cell.type)
lines <- sort(unique(obj$cell.line.ident))
covar_table <- covar_table[covar_table$cell_line %in% lines, ]

obj$age <- NA
obj$sex <- NA
obj$aff <- NA
obj$cell.type.fraction <- NA

for (l in lines) {
  print(l)
  age <- unique(covar_table$age[covar_table$cell_line == l])
  cat("\nage: ", age)
  obj$age[obj$cell.line.ident == l] <- age
  sex <- unique(covar_table$sex[covar_table$cell_line == l])
  cat("\nsex: ", sex)
  obj$sex[obj$cell.line.ident == l] <- sex
  aff <- unique(covar_table$disease[covar_table$cell_line == l])
  cat("\naff: ", aff, "\n")
  obj$aff[obj$cell.line.ident == l] <- aff
}

times <- sort(unique(obj$time.ident))
types <- sort(unique(obj$cell.type))

for (y in types) {
  for (i in times) {
    for (l in lines) {
      if (y == "GABA") {
        obj$cell.type.fraction[obj$cell.type == y & obj$time.ident == i & obj$cell.line.ident == l] <-
          covar_table$GABA_fraction[covar_table$time == i & covar_table$cell_line == l]
      } else if (y == "NEFM_neg_glut") {
        obj$cell.type.fraction[obj$cell.type == y & obj$time.ident == i & obj$cell.line.ident == l] <-
          covar_table$nmglut_fraction[covar_table$time == i & covar_table$cell_line == l]
      } else if (y == "NEFM_pos_glut") {
        obj$cell.type.fraction[obj$cell.type == y & obj$time.ident == i & obj$cell.line.ident == l] <-
          covar_table$npglut_fraction[covar_table$time == i & covar_table$cell_line == l]
      }
    }
  }
}
unique(obj$cell.type.fraction)
save(obj, file = "labeled_nfeat3000_with_covars_from_30Sep2022_single_cell_DE.RData")

GABA <- subset(obj, cell.type == "GABA")
nmglut <- subset(obj, cell.type == "NEFM_neg_glut")
npglut <- subset(obj, cell.type == "NEFM_pos_glut")

# compute deg ####
find_time_DEG <- function(obj, min.pct, type, time){
  DefaultAssay(obj) <- "RNA"
  obj <- ScaleData(obj)
  Idents(obj) <- "time.ident"
  df <- FindMarkers(object = obj, assay = "RNA", slot = "scale.data",
                         ident.1 = time, group.by = 'time.ident', 
                         ident.2 = "0hr",
                         min.pct = min.pct, logfc.threshold = 0.0,
                         test.use = "MAST", latent.vars = "orig.ident", 
                         random.seed = 42)
  df$gene_symbol <- rownames(df)
  df$cell_type <- type
  return(df)
}

GABA_res_1v0 <- find_time_DEG(GABA, 0.01, "GABA", "1hr")
GABA_res_6v0 <- find_time_DEG(GABA, 0.01, "GABA", "6hr")

nmglut_res_1v0 <- find_time_DEG(nmglut, 0.01, "nmglut", "1hr")
nmglut_res_6v0 <- find_time_DEG(nmglut, 0.01, "nmglut", "6hr")

npglut_res_1v0 <- find_time_DEG(npglut, 0.01, "npglut", "1hr")
npglut_res_6v0 <- find_time_DEG(npglut, 0.01, "npglut", "6hr")

save(npglut_res_1v0, file = "seurat_MAST_w_latentvars_18line_DE_results_npglut1v0.RData")
save(npglut_res_6v0, file = "seurat_MAST_w_latentvars_18line_DE_results_npglut6v0.RData")


# summarize counts ####
reslist_1v0 <- list(GABA_res_1v0, nmglut_res_1v0, npglut_res_1v0)
reslist_6v0 <- list(GABA_res_6v0, nmglut_res_6v0, npglut_res_6v0)
types <- c("GABA", "nmglut", "npglut")
deg_counts1 <- array(dim = c(3, 4), 
                     dimnames = list(paste(types, "1v0", sep = "-"),
                                     c("positive", "negative", "nonsig", "total"))) 
deg_counts6 <- array(dim = c(3, 4), 
                     dimnames = list(paste(types, "6v0", sep = "-"),
                                     c("positive", "negative", "nonsig", "total")))  

for (i in 1:length(reslist_1v0)) {
  deg_counts1[i,1] <- sum(reslist_1v0[[i]]$p_val_adj < 1e-8 & reslist_1v0[[i]]$avg_diff > 0)
  deg_counts1[i,2] <- sum(reslist_1v0[[i]]$p_val_adj < 1e-8 & reslist_1v0[[i]]$avg_diff < 0)
  deg_counts1[i,3] <- sum(reslist_1v0[[i]]$p_val_adj > 0.05)
  deg_counts1[i,4] <- nrow(reslist_1v0[[i]])
  deg_counts6[i,1] <- sum(reslist_6v0[[i]]$p_val_adj < 1e-8 & reslist_6v0[[i]]$avg_diff > 0)
  deg_counts6[i,2] <- sum(reslist_6v0[[i]]$p_val_adj < 1e-8 & reslist_6v0[[i]]$avg_diff < 0)
  deg_counts6[i,3] <- sum(reslist_6v0[[i]]$p_val_adj > 0.05)
  deg_counts6[i,4] <- nrow(reslist_6v0[[i]])
  
} 
deg_counts <- rbind(deg_counts1, deg_counts6)
write.table(deg_counts, 
            file = "./signac_MAST_DE/DEG_counts_summary.csv", 
            quote = F, sep = ",", row.names = T, col.names = T)


# volcano plots ####
response_genes <- c("FOS", "NPAS4", "BDNF", "VGF")
library(ggrepel)
for (i in 1:length(reslist_6v0)){
  res <- reslist_6v0[[i]]
  #res$gene.symbol <- rownames(res)
  res$significance <- "nonsignificant"
  res$significance[res$p_val_adj < 0.05 & res$avg_diff > 0] <- "up"
  res$significance[res$p_val_adj < 0.05 & res$avg_diff < 0] <- "down"
  unique(res$significance)
  res$significance <- factor(res$significance, levels = c("up", "nonsignificant", "down"))
  res$neg_log_pval <- (0 - log2(res$p_val_adj))
  res$labelling <- ""
  for (j in response_genes){
    res$labelling[res$gene %in% j] <- j
  }
  is.allnonsig <- length(unique(res$significance)) == 1
  if (is.allnonsig) {
    dotColors <- c("grey50")
  } else {
    dotColors <- c("red3", "grey50", "steelblue3")
  }
  png(paste0("./signac_MAST_DE/volcano_plots/", 
             types[i], "_6v0hr_volcano_plot.png"))
  p <- ggplot(data = as.data.frame(res),
              aes(x = avg_diff, 
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
    ggtitle(paste0(types[i]))
  print(p)
  dev.off()
}
