# Chuxuan Li 05/09/2023
# Do MAST testing the difference between anesthesia and control cells in 029 
#rabbit data in oligodendrocytes

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
unique(integrated_labeled$cell.type)
oligo <- subset(integrated_labeled, cell.type == "oligodendrocyte")
sum(oligo$cond == "control")

# MAST ####
feat_mat <- data.frame(geneid = oligo@assays$SCT@data@Dimnames[[1]], 
                       row.names = oligo@assays$SCT@data@Dimnames[[1]])
cell_mat <- oligo@meta.data # has to be a data.frame, contains metadata for each cell
ncol(oligo@assays$SCT@data) == nrow(cell_mat)
nrow(oligo@assays$SCT@data) == nrow(feat_mat)
scaRaw <- FromMatrix(exprsArray = as.array(oligo@assays$SCT@data),
                     cData = cell_mat, fData = feat_mat)
save(scaRaw, 
     file = paste0("./MAST_scDE/scaRAW_obj_oligo.RData"))
expressed_genes <- freq(scaRaw) > 0.05 # proportion of nonzero values across all cells
sca <- scaRaw[expressed_genes, ]
rm(scaRaw)
gc()
zlmCond <- zlm(~cond, sca)
save(zlmCond, file = paste0("./MAST_scDE/MAST_zlmCond_oligo.RData"))

summaryCond <- summary(zlmCond, doLRT = 'condanesthesia') 
rm(zlmCond)
gc()
save(summaryCond, 
     file = paste0("./MAST_scDE/summaryCond_oligo.RData"))

summaryDt <- summaryCond$datatable
rm(summaryCond)
gc()
fcHurdle <- merge(summaryDt[contrast=='condanesthesia' & component=='H',
                            .(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='condanesthesia' & component=='logFC', 
                            .(primerid, coef, ci.hi, ci.lo)], by='primerid')

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
setorder(fcHurdle, fdr)
save(fcHurdle, file = paste0("./MAST_scDE/results_oligo.RData"))

# plot
fcHurdle$significance <- "nonsignificant"
fcHurdle$significance[fcHurdle$fdr < 0.05 & fcHurdle$coef > 0] <- "up"
fcHurdle$significance[fcHurdle$fdr < 0.05 & fcHurdle$coef < 0] <- "down"
unique(fcHurdle$significance)
fcHurdle$significance <- factor(fcHurdle$significance, levels = c("up", "nonsignificant", "down"))
fcHurdle$neg_log_pval <- (0 - log2(fcHurdle$`Pr(>Chisq)`))
png(paste0("./MAST_scDE/volcano_plots/oligo_volcano_plot.png"), width = 400, height = 400)
p <- ggplot(data = as.data.frame(fcHurdle),
            aes(x = coef, 
                y = neg_log_pval, 
                color = significance)) + 
  geom_point(size = 0.2) +
  scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
  theme_minimal() +
  xlim(c(-1, 1)) +
  ggtitle(paste0("Differential Expression between 2 conditions in rabbit\n oligodendrocytes"))
print(p)
dev.off()

pdf(paste0("./MAST_scDE/volcano_plots/oligo_volcano_plot.pdf"))
p <- ggplot(data = as.data.frame(fcHurdle),
            aes(x = coef, 
                y = neg_log_pval, 
                color = significance)) + 
  geom_point(size = 0.2) +
  scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
  theme_minimal() +
  xlim(c(-2, 2)) +
  ggtitle(paste0("Differential Expression between 2 conditions in rabbit\n oligodendrocytes"))
print(p)
dev.off()
deg_counts <- array(dim = c(1, 4), dimnames = list("oligodendrocytes", c("total genes", "passed filter", 
                                                            "upregulated (in anesthesia)", "downregulated (in anesthesia)")))

deg_counts[1, 1] <- nrow(oligo)
deg_counts[1, 2] <- sum(expressed_genes)
deg_counts[1, 3] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef > 0, na.rm = T)
deg_counts[1, 4] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef < 0, na.rm = T)
write.table(deg_counts, file = "./MAST_scDE/deg_counts_oligo.csv")
