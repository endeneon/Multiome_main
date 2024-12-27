# Chuxuan Li 02/16/2023
# use pseudobulk samples to compute variance partition for 029 ATACseq data

# init ####
library(Seurat)
library(Signac)
library(rtracklayer)
library(GenomicRanges)
library(readr)
library(BiocParallel)
library(edgeR)
library(variancePartition)

library(readxl)
library(stats)
library(ggplot2)
library(stringr)

load("./multiomic_obj_with_new_peaks_labeled.RData")

DefaultAssay(multiomic_obj_new)

# make pseudobulk samples ####
lines <- unique(multiomic_obj_new$cell.line.ident)
Idents(multiomic_obj_new) <- "RNA.cell.type"
times <- unique(multiomic_obj_new$time.ident)
unique(multiomic_obj_new$RNA.cell.type)
multiomic_obj_new <- subset(multiomic_obj_new, RNA.cell.type != "unidentified")
types <- as.vector(unique(multiomic_obj_new$RNA.cell.type))


for (j in 1:length(types)) {
  for (k in 1:length(times)) {
    for (i in 1:length(lines)) {
      print(paste(j, k, i, sep = ";"))
      temp_count <-
        rowSums(subset(x = multiomic_obj_new,
                       idents = types[j],
                       subset = (cell.line.ident == lines[i]) &
                         (time.ident == times[k])))
      if (i == 1 & j == 1 & k == 1) {
        export_df <- data.frame(temp_count,
                                stringsAsFactors = F)
        colnames(export_df) <- paste(types[j],
                                     times[k],
                                     lines[i],
                                     sep = "_")
      } else {
        temp_colnames <- colnames(export_df)
        export_df <- data.frame(cbind(export_df,
                                      temp_count),
                                stringsAsFactors = F)
        colnames(export_df) <-
          c(temp_colnames, paste(types[j],
                                 times[k],
                                 lines[i],
                                 sep = "_"))
        
      }
    }
  }
}

rownames(export_df)

# analyse by variancePartition ####
library(readxl)
library(stats)
library(ggplot2)
library(stringr)

#data(varPartData)
load("../Duan_project_029_RNA/029_CIRM_metadata_df_repeated_for_3timepoints.RData")

# reorder expression read df to make it consistent with individual_info_df
export_df <- export_df[, order(colnames(export_df))]

# DGEList 
DGE <- 
  DGEList(counts = as.matrix(export_df),
          samples = colnames(export_df),
          remove.zeros = T)

DGE <- DGE[rowSums(cpm(DGE) > 0.5) > 10, ,keep.lib.sizes = F]
DGE <- calcNormFactors(DGE)

# construct varP uncertainty matrix for voom evaluation
voom_varP_design <- 
  model.matrix( ~ cell.type +
                  time.point +
                  batch +
                  age +
                  sex
                , meta_df)
# use voom to estimate weights of genes
voom_DGE <- voom(DGE, design = voom_varP_design)


# run the final Var calculation with Weights on (can be turned off by useWeights=F)
par_formula <- 
  ~ (1|cell.type) +
  (1|sex) +
  (1|batch) +
  #(cell.line) +
  (1|time.point) +
  age

par_DGE <- 
  fitExtractVarPartModel(exprObj = voom_DGE,
                         formula = par_formula,
                         data = meta_df,
                         BPPARAM = MulticoreParam(workers = 8,
                                                  progressbar = T))

plotVarPart(sortCols(par_DGE)) +
  ggtitle("029 ATACseq raw data variance partition")

con_par_formula <- 
  ~ (sex) +
  (batch) +
  #(cell.line) +
  (cell.type) +
  (time.point) +
  age

con_DGE <-
  canCorPairs(formula = con_par_formula,
              data = meta_df)
plotCorrMatrix(con_DGE, margins = c(12, 12), cexRow = 1, cexCol = 1) + 
  theme(text = element_text(size = 10)) 
