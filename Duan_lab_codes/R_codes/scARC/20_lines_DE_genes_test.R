# Siwei 18 Feb 2022

# test 20 lines DE genes

# init
library(Seurat)
library(Signac)
library(rtracklayer)
library(GenomicRanges)
library(readr)
library(BiocParallel)
# library(future)

library(edgeR)
# library(VariantAnnotation)
library(variancePartition)


# 20 lines RNA-Seq data
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/demux_20line_integrated_labeled_obj.RData")

# set parallel computing
library(future)
plan("multicore", workers = 2)

# get the cell identities associated to the read matrix
# unique(integrated_labeled$cell.line.ident)
integrated_labeled$time.ident


cell_line_ident <- 
  unique(integrated_labeled$cell.line.ident)

DefaultAssay(integrated_labeled) <- "RNA"

cell_line_ident <-
  cell_line_ident[order(cell_line_ident)]

Idents(integrated_labeled)
  
# Levels: GABA NEFM_pos_glut NEFM_neg_glut SEMA3E_pos_GABA unknown
# stimulation time
# 

stim_time <- c("0hr", "1hr", "6hr")
cell_type <- c("GABA", "NEFM_pos_glut", "NEFM_neg_glut")

i <- 1
j <- 1
k <- 1

for (j in 1:length(cell_type)) {
  for (k in 1:length(stim_time)) {
    for (i in 1:length(cell_line_ident)) {
      print(paste(j, k, i, sep = ";"))
      temp_count <-
        rowSums(subset(x = integrated_labeled,
                       idents = cell_type[j],
                       subset = (cell.line.ident == cell_line_ident[i]) &
                       (time.ident == stim_time[k])))
      if (i == 1 & j == 1 & k == 1) {
        export_df <- data.frame(temp_count,
                                stringsAsFactors = F)
        colnames(export_df) <- paste(cell_type[j],
                                     stim_time[k],
                                     cell_line_ident[i],
                                     sep = "_")
      } else {
        temp_colnames <- colnames(export_df)
        export_df <- data.frame(cbind(export_df,
                                      temp_count),
                                stringsAsFactors = F)
        colnames(export_df) <-
          c(temp_colnames, paste(cell_type[j],
                                 stim_time[k],
                                 cell_line_ident[i],
                                 sep = "_"))
        
      }
    }
  }
}

# colnames(export_df) <- cell_line_ident
rownames(export_df)

#### analyse by variancePartition
library(readxl)
library(stats)
library(ggplot2)
#### varPart data
data(varPartData)

individual_info_df <-
  read_excel("MGS_iPSC_20_lines_18Feb2022.xlsx")
individual_info_df$Aff <- as.factor(individual_info_df$Aff)
individual_info_df$sex <- as.factor(individual_info_df$sex)
individual_info_df$`Co-culture batch` <- as.factor(individual_info_df$`Co-culture batch`)
individual_info_df$`new iPSC line ID` <- as.factor(individual_info_df$`new iPSC line ID`)
individual_info_df$cell_type <- as.factor(individual_info_df$cell_type)
individual_info_df$stim_time <- as.factor(individual_info_df$stim_time)


# reorder expression read df to make it consistent with individual_info_df
export_df <- export_df[, order(colnames(export_df))]

DGE_20lines <- 
  DGEList(counts = as.matrix(export_df),
          samples = colnames(export_df),
          remove.zeros = T)

DGE_20lines <- DGE_20lines[rowSums(cpm(DGE_20lines) > 0.5) > 10, ,keep.lib.sizes = F]
DGE_20lines <- calcNormFactors(DGE_20lines)

# construct varP uncertainty matrix for voom evaluation
voom_varP_design <- 
  model.matrix( ~ cell_type + stim_time + `Co-culture batch`, 
                individual_info_df)
# use voom to estimate weights of genes
voom_DGE_20lines <- voom(DGE_20lines, design = voom_varP_design)

# run the final Var calculation with Weights on (can be turned off by useWeights=F)
par_formula <- 
  ~ (1|Aff) + 
  (1|sex) + 
  (1|`Co-culture batch`) + 
  (1|`new iPSC line ID`) +
  (1|cell_type) + 
  (1|stim_time) +
  Age

par_DGE_20lines <- 
  fitExtractVarPartModel(exprObj = voom_DGE_20lines,
                         formula = par_formula,
                         data = individual_info_df,
                         BPPARAM = MulticoreParam(workers = 8,
                                                  progressbar = T))

con_par_formula <- 
  ~ (Aff) + 
  (sex) + 
  (`Co-culture batch`) + 
  (`new iPSC line ID`) +
  (cell_type) + 
  (stim_time) +
  Age

con_DGE_20lines <-
  canCorPairs(formula = con_par_formula,
              data = individual_info_df)

plotCorrMatrix(con_DGE_20lines)


# make some plots
plotVarPart(par_DGE_20lines) +
  ggtitle("8 factors, batch effect not removed")

## check contributions
par_formula <- 
  # ~ (cell_type + 0|`Co-culture batch`) + 
  ~ (`Co-culture batch` + 0|cell_type) + 
  (1|Aff) + 
  (1|sex) + 
  (1|`Co-culture batch`) +
  (1|`new iPSC line ID`) +
  (1|stim_time) +
  Age

par_DGE_20lines <- 
  fitExtractVarPartModel(exprObj = voom_DGE_20lines,
                         formula = par_formula,
                         data = individual_info_df,
                         BPPARAM = MulticoreParam(workers = 8,
                                                  progressbar = T))


plotVarPart(sortCols(par_DGE_20lines))


# analyze after removed batch effect
library('limma')

fit_DGE_20lines <- 
  lmFit(object = voom_DGE_20lines,
        model.matrix(~ `Co-culture batch`, 
                     individual_info_df))

residual_DGE_20_lines <-
  residuals(fit_DGE_20lines,
            voom_DGE_20lines)

par_formula_batch_removed <-
  ~ (1|Aff) + 
  (1|sex) + 
  # (1|`Co-culture batch`) +
  (1|`new iPSC line ID`) +
  (1|cell_type) + 
  (1|stim_time) +
  Age


varPartResid_DGE20lines <-
  fitExtractVarPartModel(residual_DGE_20_lines,
                         formula = par_formula_batch_removed, 
                         individual_info_df, 
                         BPPARAM = MulticoreParam(workers = 16,
                                                  progressbar = T))

plotVarPart(sortCols(varPartResid_DGE20lines)) +
  ggtitle("Residual_analysis")

### remove culture batch
modelFit_DGE20lines <-
  fitVarPartModel(exprObj = voom_DGE_20lines,
                  formula = ~ (1|`Co-culture batch`),
                  data = individual_info_df, 
                  BPPARAM = MulticoreParam(workers = 16,
                                           progressbar = T))

residual_batch_removed <- 
  residuals(modelFit_DGE20lines)

par_formula_batch_removed <-
  ~ (1|Aff) + 
  # (1|sex) + 
  # (1|`Co-culture batch`) +
  # (1|`new iPSC line ID`) +
  (1|cell_type) + 
  (1|stim_time) +
  Age

varPartResid_DGE20lines_batch_removed <-
  fitExtractVarPartModel(residual_batch_removed,
                         formula = par_formula_batch_removed, 
                         individual_info_df, 
                         BPPARAM = MulticoreParam(workers = 16,
                                                  progressbar = T))
plotVarPart(sortCols(varPartResid_DGE20lines_batch_removed)) +
  ggtitle("Residual_analysis")
## try one-step model

residList_DGE20lines <- 
  fitExtractVarPartModel(exprObj = voom_DGE_20lines,
                         ~ (1|`Co-culture batch`),
                         individual_info_df,
                         fxn = residuals,
                         BPPARAM = MulticoreParam(workers = 16,
                                                  progressbar = T))
residList_combined_DGE20lines <-
  do.call(rbind, residList_DGE20lines)

plotVarPart(sortCols(residList_combined_DGE20lines)) +
  ggtitle("Residual_analysis")

# png(filename = "library_contrib.png",
#     width = 720, height = 480)
# plot()
# dev.off()
# 
# q("no")
