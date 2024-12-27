# Chuxuan Li 02/15/2023
# use pseudobulk samples to compute variance partition

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

# load RNAseq data
load("./029_RNA_integrated_labeled.RData")

# make pseudobulk samples ####
DefaultAssay(integrated_labeled) <- "RNA"
lines <- unique(integrated_labeled$cell.line.ident)
Idents(integrated_labeled)
times <- unique(integrated_labeled$time.ident)
integrated_labeled <- subset(integrated_labeled, cell.type != "unidentified")
cell_type <- as.vector(unique(integrated_labeled$cell.type))

# count number of cells in each line, type, and time
i <- 1
j <- 1
k <- 1

for (j in 1:length(cell_type)) {
  for (k in 1:length(times)) {
    for (i in 1:length(lines)) {
      print(paste(cell_type[j], times[k], lines[i], sep = ";"))
      temp_count <-
        rowSums(subset(x = integrated_labeled,
                       subset = (cell.line.ident == lines[i]) & cell.type == cell_type[j] &
                         (time.ident == times[k])))
      if (i == 1 & j == 1 & k == 1) {
        export_df <- data.frame(temp_count,
                                stringsAsFactors = F)
        colnames(export_df) <- paste(cell_type[j],
                                     times[k],
                                     lines[i],
                                     sep = "_")
      } else {
        temp_colnames <- colnames(export_df)
        export_df <- data.frame(cbind(export_df,
                                      temp_count),
                                stringsAsFactors = F)
        colnames(export_df) <-
          c(temp_colnames, paste(cell_type[j],
                                 times[k],
                                 lines[i],
                                 sep = "_"))
      }
    }
  }
}

colnames(export_df) 
rownames(export_df)
save(export_df, file = "029_raw_pseudobulk_matrix_for_varpartition.RData")

# make metadata ####
individual_info_df <-
  read_excel("../CIRM control iPSC lines_40_Duan (003).xlsx")
individual_info_df <- individual_info_df[, c(4, 6, 13, 16)]
colnames(individual_info_df) <- c("batch", "cell.line", "age", "sex")
individual_info_df <- individual_info_df[!is.na(individual_info_df$batch), ]
individual_info_df$sex <- str_replace(individual_info_df$sex, "Female", "F")
individual_info_df$sex <- str_replace(individual_info_df$sex, "Male", "M")
individual_info_df <- individual_info_df[individual_info_df$cell.line %in% lines, ]
individual_info_df <- individual_info_df[order(individual_info_df$cell.line), ]

meta_df <- as.data.frame(individual_info_df)
meta_df$aff <- "control"
df_to_bind <- meta_df
for (i in 1:(length(times) * length(cell_type) - 1)){
  print(i)
  meta_df <- rbind(meta_df, df_to_bind)
}
ctvec <- rep(cell_type, each = length(times)*length(lines))
length(ctvec)
tpvec <- rep(times, each = length(lines))
length(tpvec)
tpvec <- c(tpvec, tpvec, tpvec)
length(tpvec)

meta_df$cell.type <- ctvec
meta_df$cell.type <- as.factor(meta_df$cell.type)
meta_df$time.point <- tpvec
meta_df$time.point <- as.factor(meta_df$time.point)
meta_df$cell.line <- as.factor(meta_df$cell.line)
meta_df$aff <- as.factor(meta_df$aff)
meta_df$sex <- as.factor(meta_df$sex)
meta_df$batch <- as.factor(meta_df$batch)
meta_df$age <- as.numeric(meta_df$age)

save(individual_info_df, file = "029_CIRM_metadata_df.RData")
save(meta_df, file = "029_CIRM_metadata_df_repeated_for_3timepoints.RData")

# variance partition ####
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
  ggtitle("029 RNAseq raw data variance partition")

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


