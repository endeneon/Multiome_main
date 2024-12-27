# Chuxuan Li 04/20/2023
# re-run nmglut 0hr and npglut 0hr MAST scDE analysis to check if they are correctly done

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


# check at summaryCond level ####
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/MAST_scDE/summaryCond_nmglut_0hr.RData")
summaryCond_nm0 <- summaryCond
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/MAST_scDE/summaryCond_npglut_0hr.RData")
summaryCond_np0 <- summaryCond_0hr
rm(summaryCond_0hr)

summaryDt_nm0 <- summaryCond_nm0$datatable
summaryDt_np0 <- summaryCond_np0$datatable
fcHurdle_nm0 <- merge(summaryDt_nm0[contrast=='affcontrol' & component=='H',
                                .(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt_nm0[contrast=='affcontrol' & component=='logFC', 
                                .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
fcHurdle_nm0[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR

fcHurdle_np0 <- merge(summaryDt_np0[contrast=='affcontrol' & component=='H',
                                    .(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt_np0[contrast=='affcontrol' & component=='logFC', 
                                    .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
fcHurdle_np0[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
fcHurdle_np0$coef <- (-1) * fcHurdle_np0$coef
fcHurdle_nm0$coef <- (-1) * fcHurdle_nm0$coef
sum(fcHurdle_nm0$coef > 0 & fcHurdle_nm0$fdr < 0.05)
sum(fcHurdle_np0$coef > 0 & fcHurdle_np0$fdr < 0.05)

sum(fcHurdle_0hr$coef > 0 & fcHurdle_0hr$fdr < 0.05)

# redo nmglut 6hr, npglut 0hr ####
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/integrated_labeled_018-029_RNAseq_obj_with_scDE_covars.RData")
load("/nvmefs/scARC_Duan_018/case_control_DE_analysis/018-029_covar_table_final.RData")

types <- sort(unique(integrated_labeled$cell.type))
lines <- sort(unique(integrated_labeled$cell.line.ident))

nmglut <- subset(integrated_labeled, cell.type == "nmglut")
nmglut_6hr <- subset(nmglut, time.ident == "6hr")
npglut <- subset(integrated_labeled, cell.type == "npglut")
npglut_0hr <- subset(npglut, time.ident == "0hr")

objects <- list(nmglut_6hr, npglut_0hr)
for (i in 1:length(objects)) {
  timeobj <- objects[[i]]
  type <- unique(timeobj$cell.type)
  time <- unique(timeobj$time.ident)
  feat_mat <- data.frame(geneid = timeobj@assays$SCT@data@Dimnames[[1]], 
                         row.names = timeobj@assays$SCT@data@Dimnames[[1]])
  cell_mat <- timeobj@meta.data # has to be a data.frame, contains metadata for each cell
  stopifnot(ncol(timeobj@assays$SCT@data) == nrow(cell_mat), 
            nrow(timeobj@assays$SCT@data) == nrow(feat_mat))
  scaRaw <- FromMatrix(exprsArray = as.array(timeobj@assays$SCT@data),
                       cData = cell_mat, fData = feat_mat)
  save(scaRaw, 
       file = paste0("./MAST_scDE/scaRAW_obj_", type, "_", time, "redo.RData"))
  expressed_genes <- freq(scaRaw) > 0.05 # proportion of nonzero values across all cells
  sca <- scaRaw[expressed_genes, ]
  rm(scaRaw)
  gc()
  colData(sca)$aff <- factor(colData(sca)$aff)
  colData(sca)$sex <- factor(colData(sca)$sex)
  colData(sca)$orig.ident <- factor(colData(sca)$orig.ident)
  
  zlmCond <- zlm(~aff + age + sex + orig.ident + cell.type.fraction, sca)
  save(zlmCond, file = paste0("./MAST_scDE/MAST_zlmCond_", type, "_", time, "redo.RData"))
  
  summaryCond <- summary(zlmCond, doLRT = 'affcontrol', logFC = TRUE) 
  rm(zlmCond)
  gc()
  save(summaryCond, 
       file = paste0("./MAST_scDE/summaryCond_", type, "_", time, "redo.RData"))
  
  summaryDt <- summaryCond$datatable
  rm(summaryCond)
  gc()
  fcHurdle <- merge(summaryDt[contrast=='affcontrol' & component=='H',
                              .(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='affcontrol' & component=='logFC', 
                              .(primerid, coef, ci.hi, ci.lo, z)], by='primerid')
  
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
  # fcHurdle <- merge(fcHurdle[fdr < 0.05], 
  #                   as.data.table(mcols(sca)), by='primerid')
  setorder(fcHurdle, fdr)
  save(fcHurdle, file = paste0("./MAST_scDE/results_", type, "_", time, "redo.RData"))
} 

