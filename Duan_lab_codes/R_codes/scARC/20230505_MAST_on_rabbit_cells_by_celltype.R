# Chuxuan Li 05/05/2023
# Use MAST to do single cell DE analysis on rabbit cells in 029 RNAseq data

# init ####
library(ggplot2)
library(ggrepel)
library(GGally)
library(GSEABase)
library(limma)
library(reshape2)
library(data.table)
library(knitr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(NMF)
library(rsvd)
library(RColorBrewer)
library(MAST)
library(Seurat)

setwd("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/oc2_mm10_analysis")
load("./integrated_labeled_after_filtering2.RData")

# make object ####
unique(integrated_labeled$cell.type)
integrated_labeled$cond <- "anesthesia"
integrated_labeled$cond[integrated_labeled$orig.ident %in% 
                          c("20088-1", "60060-1", "60060-6")] <- "control"
integrated_labeled$cond <- factor(integrated_labeled$cond, levels = c("control", "anesthesia"))
unique(integrated_labeled$cond)
GABA <- subset(integrated_labeled, cell.type == "GABA")
glut <- subset(integrated_labeled, cell.type == "glut")
corneu <- subset(integrated_labeled, cell.type %in% c("GABA", "glut"))
subobj_list <- list(GABA, glut, corneu)

# loop MAST ####
types <- c("GABA", "glut", "GABA and glut combined")
deg_counts <- array(dim = c(3, 4), dimnames = list(types, c("total genes", "passed filter", 
                                      "upregulated (in anesthesia)", "downregulated (in anesthesia)")))
for (i in 1:length(types)) {
  type <- types[i]
  subobj <- subobj_list[[i]]
  deg_counts[i, 1] <- nrow(subobj)
  print(type)

  feat_mat <- data.frame(geneid = subobj@assays$SCT@data@Dimnames[[1]], 
                         row.names = subobj@assays$SCT@data@Dimnames[[1]])
  cell_mat <- subobj@meta.data # has to be a data.frame, contains metadata for each cell
  stopifnot(ncol(subobj@assays$SCT@data) == nrow(cell_mat), 
            nrow(subobj@assays$SCT@data) == nrow(feat_mat))
  scaRaw <- FromMatrix(exprsArray = as.array(subobj@assays$SCT@data),
                       cData = cell_mat, fData = feat_mat)
  save(scaRaw, 
       file = paste0("./MAST_scDE/scaRAW_obj_", type, ".RData"))
  expressed_genes <- freq(scaRaw) > 0.05 # proportion of nonzero values across all cells
  sca <- scaRaw[expressed_genes, ]
  rm(scaRaw)
  gc()
  deg_counts[i, 2] <- sum(expressed_genes)
  zlmCond <- zlm(~cond, sca)
  save(zlmCond, file = paste0("./MAST_scDE/MAST_zlmCond_", type, ".RData"))
  
  summaryCond <- summary(zlmCond, doLRT = 'condanesthesia') 
  rm(zlmCond)
  gc()
  save(summaryCond, 
       file = paste0("./MAST_scDE/summaryCond_", type, ".RData"))
  
  summaryDt <- summaryCond$datatable
  rm(summaryCond)
  gc()
  fcHurdle <- merge(summaryDt[contrast=='condanesthesia' & component=='H',
                              .(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='condanesthesia' & component=='logFC', 
                              .(primerid, coef, ci.hi, ci.lo)], by='primerid')
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
  setorder(fcHurdle, fdr)
  save(fcHurdle, file = paste0("./MAST_scDE/results_", type, ".RData"))
  deg_counts[i, 3] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef > 0)
  deg_counts[i, 4] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef < 0)
  
  # plot
  fcHurdle$significance <- "nonsignificant"
  fcHurdle$significance[fcHurdle$fdr < 0.05 & fcHurdle$coef > 0] <- "up"
  fcHurdle$significance[fcHurdle$fdr < 0.05 & fcHurdle$coef < 0] <- "down"
  unique(fcHurdle$significance)
  fcHurdle$significance <- factor(fcHurdle$significance, levels = c("up", "nonsignificant", "down"))
  fcHurdle$neg_log_pval <- (0 - log2(fcHurdle$`Pr(>Chisq)`))
  png(paste0("./MAST_scDE/volcano_plots/", type, "_volcano_plot.png"), width = 400, height = 400)
  p <- ggplot(data = as.data.frame(fcHurdle),
              aes(x = coef, 
                  y = neg_log_pval, 
                  color = significance)) + 
    geom_point(size = 0.2) +
    scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
    theme_minimal() +
    xlim(c(-1, 1)) +
    ggtitle(paste0("Differential Expression between 2 conditions in rabbit\n", type, " cells"))
  print(p)
  dev.off()
  
  pdf(paste0("./MAST_scDE/volcano_plots/", type, "_volcano_plot.pdf"))
  p <- ggplot(data = as.data.frame(fcHurdle),
              aes(x = coef, 
                  y = neg_log_pval, 
                  color = significance)) + 
    geom_point(size = 0.2) +
    scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
    theme_minimal() +
    xlim(c(-2, 2)) +
    ggtitle(paste0("Differential Expression between 2 conditions in rabbit ", type, " cells"))
  print(p)
  dev.off()
}


# output DEG counts table ####
write.table(deg_counts, file = "./MAST_scDE/deg_counts.csv", 
            quote = F, sep = ",", row.names = T, col.names = T)

