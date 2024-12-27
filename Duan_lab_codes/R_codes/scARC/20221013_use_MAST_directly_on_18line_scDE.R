# Chuxuan Li 10/13/2022
# Use MAST to do single cell DE analysis without going through Seurat's findmarkers()

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
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000_with_covars_from_30Sep2022_single_cell_DE.RData")

# make object ####
unique(obj$cell.line.ident)
unique(obj$cell.type)
obj <- subset(obj, cell.type != "unknown")
obj$cell.type[obj$cell.type == "NEFM_neg_glut"] <- "nmglut"
obj$cell.type[obj$cell.type == "NEFM_pos_glut"] <- "npglut"
types <- sort(unique(obj$cell.type))

# loop MAST ####
deg_counts <- array(dim = c(3, 6), 
                    dimnames = list(types,
                                    c("total", "total after filtering", 
                                      "1v0 upregulated", "1v0 downregulated", 
                                      "6v0 upregulated", "6v0 downregulated")))
for (i in 2:length(types)) {
  type <- types[i]
  subobj <- subset(obj, cell.type == type)
  deg_counts[i, 1] <- nrow(subobj)
  feat_mat <- data.frame(geneid = subobj@assays$SCT@data@Dimnames[[1]], 
                         row.names = subobj@assays$SCT@data@Dimnames[[1]])
  cell_mat <- subobj@meta.data # has to be a data.frame, contains metadata for each cell
  stopifnot(ncol(subobj@assays$SCT@data) == nrow(cell_mat), 
            nrow(subobj@assays$SCT@data) == nrow(feat_mat))
  
  # Make single cell experiment object
  scaRaw <- FromMatrix(exprsArray = as.array(subobj@assays$SCT@data), # count matrix
                       cData = cell_mat, # cell level covariates
                       fData = feat_mat # feature level covariates
  )
  save(scaRaw, file = paste0("./MAST_direct_DE/scaRAW_obj_", type, ".RData"))
  
  # filter
  expressed_genes <- freq(scaRaw) > 0.1 # proportion of nonzero values across all cells
  sca <- scaRaw[expressed_genes, ]
  deg_counts[i, 2] <- sum(expressed_genes)
  cond <- factor(colData(sca)$time.ident)
  cond <- relevel(cond, "0hr")
  colData(sca)$condition <- cond
  
  # fit
  zlmCond <- zlm(~condition + age + sex + aff + orig.ident + cell.type.fraction, sca) # time vs. , LHS is variable, RHS is predictor
  save(zlmCond, file = paste0("./MAST_direct_DE/MAST_directly_zlmCond_", type, ".RData"))
  
  # 1hr contrast
  summaryCond_1v0 <- summary(zlmCond, doLRT = 'condition1hr') 
  save(summaryCond_1v0, file = paste0("./MAST_direct_DE/summaryCond_1v0_", type, ".RData"))
  summaryDt <- summaryCond_1v0$datatable
  fcHurdle_1v0 <- merge(summaryDt[contrast=='condition1hr' & component=='H',
                                  .(primerid, `Pr(>Chisq)`)], #hurdle P values
                        summaryDt[contrast=='condition1hr' & component=='logFC', 
                                  .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle_1v0[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
  save(fcHurdle_1v0, file = paste0("./MAST_direct_DE/result_1v0_", type, ".RData"))
  
  # 6v0 contrast
  summaryCond_6v0 <- summary(zlmCond, doLRT = 'condition6hr') 
  save(summaryCond_6v0, file = paste0("./MAST_direct_DE/summaryCond_6v0_", type, ".RData"))
  summaryDt <- summaryCond_6v0$datatable
  fcHurdle_6v0 <- merge(summaryDt[contrast=='condition6hr' & component=='H',
                                  .(primerid, `Pr(>Chisq)`)], #hurdle P values
                        summaryDt[contrast=='condition6hr' & component=='logFC', 
                                  .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle_6v0[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
  save(fcHurdle_6v0, file = paste0("./MAST_direct_DE/result_6v0_", type, ".RData"))
  
  # summary
  deg_counts[i, 3] <- sum(fcHurdle_1v0$fdr < 0.05 & fcHurdle_1v0$coef > 0)
  deg_counts[i, 4] <- sum(fcHurdle_1v0$fdr < 0.05 & fcHurdle_1v0$coef < 0)
  deg_counts[i, 5] <- sum(fcHurdle_6v0$fdr < 0.05 & fcHurdle_6v0$coef > 0)
  deg_counts[i, 6] <- sum(fcHurdle_6v0$fdr < 0.05 & fcHurdle_6v0$coef < 0)
  
  # plot
  fcHurdle_1v0$significance <- "nonsignificant"
  fcHurdle_1v0$significance[fcHurdle_1v0$fdr < 0.05 & fcHurdle_1v0$coef > 0] <- "up"
  fcHurdle_1v0$significance[fcHurdle_1v0$fdr < 0.05 & fcHurdle_1v0$coef < 0] <- "down"
  unique(fcHurdle_1v0$significance)
  fcHurdle_1v0$significance <- factor(fcHurdle_1v0$significance, levels = c("up", "nonsignificant", "down"))
  fcHurdle_1v0$neg_log_pval <- (0 - log2(fcHurdle_1v0$`Pr(>Chisq)`))
  fcHurdle_1v0$labelling <- ""
  for (j in c("FOS", "NPAS4", "IGF1", "NR4A1","FOSB", "EGR1", "BDNF")){
    fcHurdle_1v0$labelling[fcHurdle_1v0$primerid %in% j] <- j
  }
  jpeg(paste0("./MAST_direct_DE/volcano_plots/", type, "1v0_volcano_plot.jpeg"))
  p <- ggplot(data = as.data.frame(fcHurdle_1v0),
              aes(x = coef, 
                  y = neg_log_pval, 
                  color = significance,
                  label = labelling)) + 
    geom_point(size = 0.2) +
    scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
    theme_minimal() +
    geom_text_repel(box.padding = unit(0.05, 'lines'),
                    min.segment.length = 0,
                    force = 0.2,
                    max.overlaps = Inf,
                    force_pull = 0.1,
                    show.legend = F) +
    xlim(c(-2, 2)) +
    ggtitle(paste0(type, " 1v0"))
  print(p)
  dev.off()
  
  # plot 6v0
  fcHurdle_6v0$significance <- "nonsignificant"
  fcHurdle_6v0$significance[fcHurdle_6v0$fdr < 0.05 & fcHurdle_6v0$coef > 0] <- "up"
  fcHurdle_6v0$significance[fcHurdle_6v0$fdr < 0.05 & fcHurdle_6v0$coef < 0] <- "down"
  unique(fcHurdle_6v0$significance)
  fcHurdle_6v0$significance <- factor(fcHurdle_6v0$significance, levels = c("up", "nonsignificant", "down"))
  fcHurdle_6v0$neg_log_pval <- (0 - log2(fcHurdle_6v0$`Pr(>Chisq)`))
  fcHurdle_6v0$labelling <- ""
  for (j in c("FOS", "NPAS4", "IGF1", "NR4A1","FOSB", "EGR1", "BDNF")){
    fcHurdle_6v0$labelling[fcHurdle_6v0$primerid %in% j] <- j
  }
  jpeg(paste0("./MAST_direct_DE/volcano_plots/", type, "6v0_volcano_plot.jpeg"))
  p <- ggplot(data = as.data.frame(fcHurdle_6v0),
              aes(x = coef, 
                  y = neg_log_pval, 
                  color = significance,
                  label = labelling)) + 
    geom_point(size = 0.2) +
    scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
    theme_minimal() +
    geom_text_repel(box.padding = unit(0.05, 'lines'),
                    min.segment.length = 0,
                    force = 0.2,
                    max.overlaps = Inf,
                    force_pull = 0.1,
                    show.legend = F) +
    xlim(c(-2, 2)) +
    ggtitle(paste0(type, " 6v0"))
  print(p)
  dev.off()
}
write.table(deg_counts, file = "./MAST_direct_DE/deg_counts.csv", 
            quote = F, sep = ",", row.names = T, col.names = T)

# test run ####
npglut <- subset(obj, cell.type == "NEFM_pos_glut")
feat_mat <- data.frame(geneid = npglut@assays$SCT@data@Dimnames[[1]], 
                       row.names = npglut@assays$SCT@data@Dimnames[[1]])
cell_mat <- npglut@meta.data # has to be a data.frame, contains metadata for each cell
ncol(npglut@assays$SCT@data) # 41536 cells
nrow(cell_mat) # 41536
nrow(npglut@assays$SCT@data) # 29546 genes
nrow(feat_mat) # 29546
scaRaw <- FromMatrix(exprsArray = as.array(npglut@assays$SCT@data), # count matrix
                     cData = cell_mat, # cell level covariates
                     fData = feat_mat # feature level covariates
                     )
save(scaRaw, file = "./MAST_direct_DE/scaRaw_obj_npglut.RData")
aheatmap(x = assay(scaRaw[1:1000,]), 
         labRow = '', 
         annCol = as.data.frame(colData(scaRaw)[,c('time.ident', 'orig.ident', 'cell.type')]), 
         distfun = 'spearman')

# filtering
expressed_genes <- freq(scaRaw) > 0.1 # proportion of nonzero values across all cells
sca <- scaRaw[expressed_genes, ]
thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 10, min_per_bin = 20)
par(mfrow = c(5,4))
plot(thres)


# differential expression ####
cond <- factor(colData(sca)$time.ident)
cond <- relevel(cond, "0hr")
colData(sca)$condition <- cond

# For each gene in sca, fits the hurdle model in formula (linear for et>0), 
#logistic for et==0 vs et>0. Return an object of class ZlmFit containing slots 
#giving the coefficients, variance-covariance matrices, etc. After each gene, 
#optionally run the function on the fit named by 'hook'
zlmCond <- zlm(~condition + age + sex + aff + orig.ident + cell.type.fraction, sca) # time vs. , LHS is variable, RHS is predictor
save(zlmCond, file = "./MAST_direct_DE/MAST_directly_zlmCond_npglut.RData")

# 1hr ####
#only test the condition coefficient.
summaryCond_1v0 <- summary(zlmCond, doLRT = 'condition1hr') 
save(summaryCond_1v0, file = "./MAST_direct_DE/summaryCond_1v0_npglut.RData")

#print the top 4 genes by contrast using the logFC
print(summaryCond_1v0, n = 4)

# make data.table with coeff, standar errors
summaryDt <- summaryCond_1v0$datatable
fcHurdle_1v0 <- merge(summaryDt[contrast=='condition1hr' & component=='H',
                            .(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='condition1hr' & component=='logFC', 
                            .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle_1v0[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
fcHurdleSig_1v0 <- merge(fcHurdle_1v0[fdr < 0.05], 
                         as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig_1v0, fdr)
save(fcHurdle_1v0, file = "./MAST_direct_DE/result_1v0_npglut.RData")

# 6hr ####
#only test the condition coefficient.
summaryCond_6v0 <- summary(zlmCond, doLRT = 'condition6hr') 
save(summaryCond_6v0, file = "./MAST_direct_DE/summaryCond_6v0_npglut.RData")

#print the top 4 genes by contrast using the logFC
print(summaryCond_6v0, n = 4)

# make data.table with coeff, standar errors
summaryDt <- summaryCond_6v0$datatable
fcHurdle_6v0 <- merge(summaryDt[contrast=='condition6hr' & component=='H',
                                .(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='condition6hr' & component=='logFC', 
                                .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle_6v0[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
fcHurdleSig_6v0 <- merge(fcHurdle_6v0[fdr < 0.05], 
                         as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig_6v0, fdr)
save(fcHurdle_6v0, file = "./MAST_direct_DE/result_6v0_npglut.RData")



# visualization
entrez_to_plot <- fcHurdleSig[1:16, primerid]
genes_to_plot <- fcHurdleSig[1:16, geneid]
flat_dat <- as(sca[entrez_to_plot,], 'data.table')
ggbase <- ggplot(flat_dat, aes(x=condition, y=thresh,color=condition)) + 
  geom_jitter() +
  facet_wrap(~symbolid, scale='free_y') +
  ggtitle("DE Genes")
ggbase +
  geom_violin() 

mat_to_plot <- assay(sca[entrez_to_plot, ])
rownames(mat_to_plot) <- genes_to_plot
aheatmap(mat_to_plot,
         annCol = colData(sca)[ ,"condition"],
         main = "DE genes",
         col = rev(colorRampPalette(colors = brewer.pal(name = "PiYG", n = 10))(20)))
