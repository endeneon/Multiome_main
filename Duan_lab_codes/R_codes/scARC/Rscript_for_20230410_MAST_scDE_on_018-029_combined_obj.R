# Chuxuan Li 04/14/2023
# script to run MAST on all cell type x time point combinations

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

load("/nvmefs/scARC_Duan_018/case_control_DE_analysis/018-029_covar_table_final.RData")
load("./integrated_labeled_018-029_RNAseq_obj_with_scDE_covars.RData")

# variables
types <- sort(unique(integrated_labeled$cell.type))
lines <- sort(unique(integrated_labeled$cell.line.ident))
times <- sort(unique(integrated_labeled$time.ident))

# loop MAST ####
deg_counts <- array(dim = c(3, 10), 
                    dimnames = list(types,
                                    c("total genes", "passed filter - 0hr", 
                                      "upregulated - 0hr", "downregulated - 0hr",
                                      "passed filter - 1hr", 
                                      "upregulated - 1hr", "downregulated - 1hr",
                                      "passed filter - 6hr", 
                                      "upregulated - 6hr", "downregulated - 6hr")))
for (i in 1:length(types)) {
  type <- types[i]
  subobj <- subset(integrated_labeled, cell.type == type)
  deg_counts[i, 1] <- nrow(subobj)
  
  zero <- subset(subobj, time.ident == "0hr")
  one <- subset(subobj, time.ident == "1hr")
  six <- subset(subobj, time.ident == "6hr")
  timeobj_list <- list(zero, one, six)
  rm(zero, one, six)
  
  print(type)
  
  for (j in 1:length(timeobj_list)) {
    if (i == 3 & j == 1) {
      print("npglut_0hr met")
    } else {
      time <- times[j]
      print(time)
      timeobj <- timeobj_list[[j]]
      feat_mat <- data.frame(geneid = timeobj@assays$SCT@data@Dimnames[[1]], 
                             row.names = timeobj@assays$SCT@data@Dimnames[[1]])
      cell_mat <- timeobj@meta.data # has to be a data.frame, contains metadata for each cell
      stopifnot(ncol(timeobj@assays$SCT@data) == nrow(cell_mat), 
                nrow(timeobj@assays$SCT@data) == nrow(feat_mat))
      scaRaw <- FromMatrix(exprsArray = as.array(timeobj@assays$SCT@data),
                           cData = cell_mat, fData = feat_mat)
      save(scaRaw, 
           file = paste0("./MAST_scDE/scaRAW_obj_", type, "_", time, ".RData"))
      expressed_genes <- freq(scaRaw) > 0.05 # proportion of nonzero values across all cells
      sca <- scaRaw[expressed_genes, ]
      rm(scaRaw)
      gc()
      deg_counts[i, 3 * j - 1] <- sum(expressed_genes)
      colData(sca)$aff <- factor(colData(sca)$aff)
      colData(sca)$sex <- factor(colData(sca)$sex)
      colData(sca)$orig.ident <- factor(colData(sca)$orig.ident)
      
      zlmCond <- zlm(~aff + age + sex + orig.ident + cell.type.fraction, sca)
      save(zlmCond, file = paste0("./MAST_scDE/MAST_zlmCond_", type, "_", time, ".RData"))
      
      summaryCond <- summary(zlmCond, doLRT = 'affcontrol') 
      rm(zlmCond)
      gc()
      save(summaryCond, 
           file = paste0("./MAST_scDE/summaryCond_", type, "_", time, ".RData"))
      
      summaryDt <- summaryCond$datatable
      rm(summaryCond)
      gc()
      fcHurdle <- merge(summaryDt[contrast=='affcontrol' & component=='H',
                                  .(primerid, `Pr(>Chisq)`)], #hurdle P values
                        summaryDt[contrast=='affcontrol' & component=='logFC', 
                                  .(primerid, coef, ci.hi, ci.lo)], by='primerid')
      
      fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
      # fcHurdle <- merge(fcHurdle[fdr < 0.05], 
      #                   as.data.table(mcols(sca)), by='primerid')
      setorder(fcHurdle, fdr)
      save(fcHurdle, file = paste0("./MAST_scDE/results_", type, "_", time, ".RData"))
      deg_counts[i, 3 * j] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef > 0)
      deg_counts[i, 3 * j + 1] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef < 0)
      
      # plot
      fcHurdle$significance <- "nonsignificant"
      fcHurdle$significance[fcHurdle$fdr < 0.05 & fcHurdle$coef > 0] <- "up"
      fcHurdle$significance[fcHurdle$fdr < 0.05 & fcHurdle$coef < 0] <- "down"
      unique(fcHurdle$significance)
      fcHurdle$significance <- factor(fcHurdle$significance, levels = c("up", "nonsignificant", "down"))
      fcHurdle$neg_log_pval <- (0 - log2(fcHurdle$`Pr(>Chisq)`))
      png(paste0("./MAST_scDE/volcano_plots/", type, "_", time, 
                 "_casevsctrl_volcano_plot.png"), width = 700, height = 700)
      p <- ggplot(data = as.data.frame(fcHurdle),
                  aes(x = coef, 
                      y = neg_log_pval, 
                      color = significance)) + 
        geom_point(size = 0.2) +
        scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
        theme_minimal() +
        xlim(c(-2, 2)) +
        ggtitle(paste0("Case vs Control DE - ", type, " ", time))
      print(p)
      dev.off()
      
      pdf(paste0("./MAST_scDE/volcano_plots/", type, "_", time, 
                 "_casevsctrl_volcano_plot.pdf"))
      p <- ggplot(data = as.data.frame(fcHurdle),
                  aes(x = coef, 
                      y = neg_log_pval, 
                      color = significance)) + 
        geom_point(size = 0.2) +
        scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
        theme_minimal() +
        xlim(c(-2, 2)) +
        ggtitle(paste0("Case vs Control DE - ", type, " ", time))
      print(p)
      dev.off()
    }
  }
}
