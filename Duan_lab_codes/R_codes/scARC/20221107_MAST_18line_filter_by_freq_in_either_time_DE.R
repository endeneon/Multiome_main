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
unique(obj$cell.type)
obj$cell.type <- as.vector(obj$cell.type)
obj$cell.type[obj$cell.type == "NEFM_neg_glut"] <- "nmglut"
obj$cell.type[obj$cell.type == "NEFM_pos_glut"] <- "npglut"
obj <- subset(obj, cell.type != "unknown")
types <- sort(unique(obj$cell.type))

# a simple freq() function for Seurat objects ####
library(future.apply)
plan(multisession, workers = 2)
options(future.globals.maxSize = 3221225472)
pctCellsExpressed <- function(sobj) {
  # calculate the percentage of cells expressing each gene in a Seurat object
  # input: sobj, a Seurat object
  mat <- as.matrix(GetAssayData(object = sobj, assay = "SCT", slot = "data"))
  print(nrow(mat))
  rowPctExp <- function(row) {
    return(sum(row != 0) / length(row) * 100)
  }
  pct <- future_apply(mat, MARGIN = 1, FUN = rowPctExp)
  return(pct)
}

# loop MAST ####
deg_counts <- array(dim = c(3, 6), 
                    dimnames = list(types,
                                    c("total", "total after filtering", 
                                      "1v0 upregulated", "1v0 downregulated", 
                                      "6v0 upregulated", "6v0 downregulated")))
nonzero_times <- c("1hr", "6hr")
scaRaw_paths <- sort(list.files("./MAST_direct_DE", pattern = ".RData", full.names = T, recursive = F))
for (i in 1:length(types)) {
  type <- types[i]
  subobj <- subset(obj, cell.type == type)
  deg_counts[i, 1] <- nrow(subobj)
  
  # load single cell experiment object
  load(scaRaw_paths[i])

  for (j in 1:2) {
    zero <- subset(subobj, time.ident == "0hr")
    nonzero <- subset(subobj, time.ident == nonzero_times[j])
    # filter: at least 0.1 in either time point
    expressed_genes <- pctCellsExpressed(zero) > 5 | pctCellsExpressed(nonzero) > 5
    deg_counts[i, 2] <- sum(expressed_genes)
    sca <- scaRaw[expressed_genes, ]
    cond <- factor(colData(sca)$time.ident)
    cond <- relevel(cond, "0hr")
    colData(sca)$condition <- cond
    # fit
    zlmCond <- zlm(~condition + age + sex + aff + orig.ident + cell.type.fraction, sca) # time vs. , LHS is variable, RHS is predictor
    save(zlmCond, file = paste0("./MAST_direct_DE/filter_by_5pct_in_either_time/zlmCond_", nonzero_times[j], 
                                "_", type, ".RData"))
    
    # contrast
    dolrt <- paste0('condition', nonzero_times[j])
    summaryCond <- summary(zlmCond, doLRT = dolrt) 
    save(summaryCond, file = paste0("./MAST_direct_DE/filter_by_5pct_in_either_time/summaryCond_", 
                                    str_remove(nonzero_times[j], "hr"), "v0_", 
                                    type, ".RData"))
    summaryDt <- summaryCond$datatable
    fcHurdle <- merge(summaryDt[contrast == dolrt & component == 'H',
                                    .(primerid, `Pr(>Chisq)`)], #hurdle P values
                          summaryDt[contrast == dolrt & component == 'logFC', 
                                    .(primerid, coef, ci.hi, ci.lo)], by = 'primerid') #logFC coefficients
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
    save(fcHurdle, file = paste0("./MAST_direct_DE/filter_by_5pct_in_either_time/result_", 
                                 str_remove(nonzero_times[j], "hr"), "v0_", type, ".RData"))
    
    # summary
    deg_counts[i, 2 * j + 1] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef > 0)
    deg_counts[i, 2 * j + 2] <- sum(fcHurdle$fdr < 0.05 & fcHurdle$coef < 0)

    # plot
    fcHurdle$significance <- "nonsignificant"
    fcHurdle$significance[fcHurdle$fdr < 0.05 & fcHurdle$coef > 0] <- "up"
    fcHurdle$significance[fcHurdle$fdr < 0.05 & fcHurdle$coef < 0] <- "down"
    unique(fcHurdle$significance)
    fcHurdle$significance <- factor(fcHurdle$significance, levels = c("up", "nonsignificant", "down"))
    fcHurdle$neg_log_pval <- (0 - log2(fcHurdle$`Pr(>Chisq)`))
    fcHurdle$labelling <- ""
    for (gene in c("FOS", "NPAS4", "IGF1", "NR4A1","FOSB", "EGR1", "BDNF")){
      fcHurdle$labelling[fcHurdle$primerid %in% gene] <- gene
    }
    png(paste0("./MAST_direct_DE/filter_by_5pct_in_either_time/volcano_plots/", type, 
                str_remove(nonzero_times[j], "hr"), "v0_volcano_plot.png"))
    p <- ggplot(data = as.data.frame(fcHurdle),
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
      ggtitle(paste0(type, "\t", str_remove(nonzero_times[j], "hr"), "v0"))
    print(p)
    dev.off()
  }

}
write.table(deg_counts, file = "./MAST_direct_DE/filter_by_5pct_in_either_time//deg_counts.csv", 
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
fcHurdle <- merge(summaryDt[contrast=='condition1hr' & component=='H',
                                .(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast=='condition1hr' & component=='logFC', 
                                .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')] # FDR
fcHurdleSig_1v0 <- merge(fcHurdle[fdr < 0.05], 
                         as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig_1v0, fdr)
save(fcHurdle, file = "./MAST_direct_DE/result_1v0_npglut.RData")

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

# output files ####
results <- list.files(path = "./MAST_direct_DE/filter_by_freq_in_either_time/",
           full.names = T, pattern = "result*")
for (i in 1:length(results)) {
  load(results[i])
  write.table(fcHurdle, file = str_replace(str_replace(results[i], "RData", "csv"), "result_", "summary_csv/full_result_"),
              sep = ",", row.names = F, col.names = c("Gene", "p-value", "coef", "ci.hi", "ci.lo", "FDR"))
  write.table(fcHurdle[fcHurdle$coef > 0 & fcHurdle$fdr < 0.05, ], 
              file = str_replace(str_replace(results[i], "RData", "csv"), "result_", 
                                 "summary_csv/upregulated_significant/upregulated_"), row.names = F, 
              col.names = c("Gene", "p-value", "coef", "ci.hi", "ci.lo", "FDR"), sep = ",")
  write.table(fcHurdle[fcHurdle$coef < 0 & fcHurdle$fdr < 0.05, ], 
              file = str_replace(str_replace(results[i], "RData", "csv"), "result_", 
                                 "summary_csv/downregulated_significant/downregulated_"), row.names = F, 
              col.names = c("Gene", "p-value", "coef", "ci.hi", "ci.lo", "FDR"), sep = ",")
   
}
# visualization ####
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
