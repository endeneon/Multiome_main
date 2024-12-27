 # Chuxuan Li 04/10/2023
# Use MAST to do single cell DE analysis with 018-029 data combined

# init ####
{
  library(ggplot2)
  library(ggrepel)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
  library(Seurat)
}

load("./018-029_RNA_integrated_labeled_with_harmony.RData")
load("/nvmefs/scARC_Duan_018/case_control_DE_analysis/018-029_covar_table_final.RData")

#options(mc.cores = 30)

# make object ####
unique(integrated_labeled$cell.type)
integrated_labeled$cell.type <- as.vector(integrated_labeled$cell.type)
integrated_labeled <- subset(integrated_labeled, cell.type != "unidentified" |
                               cell.type != 'glut?')
types <- sort(unique(integrated_labeled$cell.type))
lines <- sort(unique(integrated_labeled$cell.line.ident))

integrated_labeled$age <- NA
integrated_labeled$sex <- NA
integrated_labeled$aff <- NA
integrated_labeled$cell.type.fraction <- NA

for (l in lines) {
  print(l)
  age <- unique(covar_table_final$age[covar_table_final$cell_line == l])
  cat("\nage: ", age)
  integrated_labeled$age[integrated_labeled$cell.line.ident == l] <- age
  sex <- unique(covar_table_final$sex[covar_table_final$cell_line == l])
  cat("\nsex: ", sex)
  integrated_labeled$sex[integrated_labeled$cell.line.ident == l] <- sex
  aff <- unique(covar_table_final$aff[covar_table_final$cell_line == l])
  cat("\naff: ", aff, "\n")
  integrated_labeled$aff[integrated_labeled$cell.line.ident == l] <- aff
}

times <- sort(unique(integrated_labeled$time.ident))
for (y in types) {
  for (i in times) {
    for (l in lines) {
      if (y == "GABA") {
        integrated_labeled$cell.type.fraction[integrated_labeled$cell.type == y & integrated_labeled$time.ident == i & integrated_labeled$cell.line.ident == l] <-
          covar_table_final$GABA_fraction[covar_table_final$time == i & covar_table_final$cell_line == l]
      } else if (y == "nmglut") {
        integrated_labeled$cell.type.fraction[integrated_labeled$cell.type == y & integrated_labeled$time.ident == i & integrated_labeled$cell.line.ident == l] <-
          covar_table_final$nmglut_fraction[covar_table_final$time == i & covar_table_final$cell_line == l]
      } else if (y == "npglut") {
        integrated_labeled$cell.type.fraction[integrated_labeled$cell.type == y & integrated_labeled$time.ident == i & integrated_labeled$cell.line.ident == l] <-
          covar_table_final$npglut_fraction[covar_table_final$time == i & covar_table_final$cell_line == l]
      }
    }
  }
}
unique(integrated_labeled$cell.type.fraction)
save(integrated_labeled, file = "integrated_labeled_018-029_RNAseq_obj_with_scDE_covars.RData")


# a simple freq() function for Seurat objects ####
filterGenes <- function(sobj, thres) {
  # calculate the percentage of cells expressing each gene in a Seurat object
  # input: sobj, a Seurat object
  ct_mat <- GetAssayData(object = sobj, assay = "SCT", slot = "data")
  print(nrow(ct_mat))
  passfilter <- (rowSums(ct_mat) / ncol(ct_mat)) > thres
  return(passfilter)
}

# loop MAST ####
deg_counts <- array(dim = c(3, 10), 
                    dimnames = list(types,
                                    c("total genes", 
                                      "passed filter - 0hr", 
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
      fcHurdle <- merge(summaryDt[contrast == 'affcontrol' & 
                                    component == 'H',
                                  .(primerid, `Pr(>Chisq)`)], #hurdle P values
                        summaryDt[contrast == 'affcontrol' & 
                                    component == 'logFC', 
                                  .(primerid, coef, ci.hi, ci.lo)], 
                        by = 'primerid')
      
      fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
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
  

# test run ####
npglut <- subset(integrated_labeled, cell.type == "npglut")
# filtering
zero <- subset(npglut, time.ident == "0hr")
one <- subset(npglut, time.ident == "1hr")
six <- subset(npglut, time.ident == "6hr")
feat_mat <- data.frame(geneid = zero@assays$SCT@data@Dimnames[[1]], 
                       row.names = zero@assays$SCT@data@Dimnames[[1]])
cell_mat <- zero@meta.data # has to be a data.frame, contains metadata for each cell
ncol(zero@assays$SCT@data) # 128392 cells
nrow(cell_mat) # 128392
nrow(zero@assays$SCT@data) # 32823 genes
nrow(feat_mat) # 32823
scaRaw <- FromMatrix(exprsArray = as.array(zero@assays$SCT@data), # count matrix
                     cData = cell_mat, # cell level covariates
                     fData = feat_mat # feature level covariates
)

# filter: at least 0.05 in either time point
passfilter <- filterGenes(zero, 0.05)
sum(passfilter) #10650
sca <- scaRaw[passfilter, ]
thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 10, min_per_bin = 20)
par(mfrow = c(5,4))
plot(thres)

# differential expression 
colData(sca)$aff <- factor(colData(sca)$aff)
colData(sca)$sex <- factor(colData(sca)$sex)
colData(sca)$orig.ident <- factor(colData(sca)$orig.ident)

zlmCond <- zlm(~aff + age + sex + orig.ident + cell.type.fraction, sca) # time vs. , LHS is variable, RHS is predictor
save(zlmCond, file = "./MAST_scDE/MAST_zlmCond_npglut_0hr_filter_5pct.RData")

summaryCond_0hr <- summary(zlmCond, doLRT = 'affcontrol') 
save(summaryCond_0hr, file = "./MAST_scDE/summaryCond_npglut_0hr.RData")
print(summaryCond_0hr, n = 4)

# make data.table with coeff, standar errors
summaryDt <- summaryCond_0hr$datatable
fcHurdle_0hr <- merge(summaryDt[contrast=='affcontrol' & component=='H',
                                .(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='affcontrol' & component=='logFC', 
                                .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle_0hr[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
# fcHurdle_0hr <- merge(fcHurdle_0hr[fdr < 0.05], 
#                          as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdle_0hr, fdr)
save(fcHurdle_0hr, file = "./MAST_scDE/result_0hr_npglut.RData")

# count DEGs 
deg_counts[3, 2] <- sum(passfilter)
deg_counts[3, 3] <- sum(fcHurdle_0hr$fdr < 0.05 & fcHurdle_0hr$coef > 0)
deg_counts[3, 4] <- sum(fcHurdle_0hr$fdr < 0.05 & fcHurdle_0hr$coef < 0)

# plot
fcHurdle_0hr$significance <- "nonsignificant"
fcHurdle_0hr$significance[fcHurdle_0hr$fdr < 0.05 & fcHurdle_0hr$coef > 0] <- "up"
fcHurdle_0hr$significance[fcHurdle_0hr$fdr < 0.05 & fcHurdle_0hr$coef < 0] <- "down"
unique(fcHurdle_0hr$significance)
fcHurdle_0hr$significance <- factor(fcHurdle_0hr$significance, levels = c("up", "nonsignificant", "down"))
fcHurdle_0hr$neg_log_pval <- (0 - log2(fcHurdle_0hr$`Pr(>Chisq)`))
png("./MAST_scDE/volcano_plots/npglut_0hr_casevsctrl_volcano_plot.png", 
     width = 700, height = 700)
p <- ggplot(data = as.data.frame(fcHurdle_0hr),
            aes(x = coef, 
                y = neg_log_pval, 
                color = significance)) + 
  geom_point(size = 0.2) +
  scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
  theme_minimal() +
  xlim(c(-1, 1)) +
  ggtitle(paste0("Case vs Control DE - npglut_0hr"))
print(p)
dev.off()

pdf("./MAST_scDE/volcano_plots/npglut_0hr_casevsctrl_volcano_plot.pdf")
p <- ggplot(data = as.data.frame(fcHurdle_0hr),
            aes(x = coef, 
                y = neg_log_pval, 
                color = significance)) + 
  geom_point(size = 0.2) +
  scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
  theme_minimal() +
  xlim(c(-1, 1)) +
  ggtitle(paste0("Case vs Control DE - npglut_0hr"))
print(p)
dev.off()

# output DEG counts table ####
write.table(deg_counts, file = "./MAST_scDE/deg_counts.csv", 
            quote = F, sep = ",", row.names = T, col.names = T)

