# Chuxuan Li 11/17/2022
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

# load RNAseq data
load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/integrated_labeled.RData")

# RNAseq ####
# make pseudobulk samples
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
save(export_df, file = "025_raw_pseudobulk_matrix_for_varpartition.RData")
# analyse by variancePartition
library(readxl)
library(stats)
library(ggplot2)
library(stringr)

#data(varPartData)
individual_info_df <-
  read_excel("./MGS_iPSC lines_60samples_scRNA_ATAC-seq_bulkATAC-seq status_AK_HZ.xlsx")
individual_info_df$cell.line <- str_replace(individual_info_df$`new iPSC line ID`, "00000", "_")
individual_info_df <- individual_info_df[individual_info_df$cell.line %in% lines, ]
individual_info_df$Aff <- individual_info_df$`Aff (SZ/ctrl, AD irrelevant)`
individual_info_df$`Aff (SZ/ctrl, AD irrelevant)` <- NULL
individual_info_df$sex <- individual_info_df$`sex (1=M, 2=F)`
individual_info_df$`sex (1=M, 2=F)` <- NULL
meta_df <- as.data.frame(cbind(individual_info_df$Age, individual_info_df$`Co-culture batch`,
                 individual_info_df$cell.line, individual_info_df$Aff, individual_info_df$sex))
colnames(meta_df) <- c("Age", "batch", "cell.line", "aff", "sex")
meta_df$batch <- str_remove(meta_df$batch, "[0-9]+ & ")
meta_df$batch <- str_remove(meta_df$batch, "[0-9]+ \\(need redo\\)  & ") 
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
meta_df$Age <- as.numeric(meta_df$Age)

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
  model.matrix( ~ aff +
                  cell.type +
                  time.point +
                  batch +
                  Age +
                  sex
                , meta_df)
# use voom to estimate weights of genes
voom_DGE <- voom(DGE, design = voom_varP_design)

# run the final Var calculation with Weights on (can be turned off by useWeights=F)
par_formula <- 
  ~ (1|aff) + 
  (1|sex) +
  (1|batch) +
  #(cell.line) +
  (1|cell.type) +
  (1|time.point) +
  Age
par_DGE <- 
  fitExtractVarPartModel(exprObj = voom_DGE,
                         formula = par_formula,
                         data = meta_df,
                         BPPARAM = MulticoreParam(workers = 8,
                                                  progressbar = T))
plotVarPart(sortCols(par_DGE)) +
  ggtitle("18-line RNAseq raw data")

con_par_formula <- 
  ~ (aff) + 
  (sex) +
  (batch) +
  #(cell.line) +
  (cell.type) +
  (time.point) +
  Age

con_DGE <-
  canCorPairs(formula = con_par_formula,
              data = meta_df)
plotCorrMatrix(con_DGE) + theme(text = element_text(size = 10))

# ATACseq ####
load("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/18line_10xaggred_labeled_by_RNA_obj.RData")

# make pseudobulk samples ####
lines <- unique(ATAC_new$cell.line.ident)
Idents(ATAC_new) <- "cell.type"
times <- unique(ATAC_new$time.ident)
unique(ATAC_new$cell.type)
ATAC_new <- subset(ATAC_new, cell.type != "unknown")
cell_type <- as.vector(unique(ATAC_new$cell.type))

i <- 1
j <- 1
k <- 1

for (j in 1:length(cell_type)) {
  for (k in 1:length(times)) {
    for (i in 1:length(lines)) {
      print(paste(j, k, i, sep = ";"))
      temp_count <-
        rowSums(subset(x = ATAC_new,
                       idents = cell_type[j],
                       subset = (cell.line.ident == lines[i]) &
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

# colnames(export_df) <- cell_line_ident
rownames(export_df)

# analyse by variancePartition ####
library(readxl)
library(stats)
library(ggplot2)
library(stringr)

#data(varPartData)

individual_info_df <-
  read_excel("MGS_iPSC lines_60samples_scRNA_ATAC-seq_bulkATAC-seq status_AK_HZ.xlsx")
individual_info_df$cell.line <- str_replace(individual_info_df$`new iPSC line ID`, "00000", "_")
individual_info_df <- individual_info_df[individual_info_df$cell.line %in% lines, ]
individual_info_df$Aff <- as.factor(individual_info_df$Aff)
individual_info_df$`Co-culture batch` <- as.factor(individual_info_df$`Co-culture batch`)
individual_info_df$`new iPSC line ID` <- as.factor(individual_info_df$`new iPSC line ID`)
individual_info_df$`sex (1=M, 2=F)` <- as.factor(individual_info_df$`sex (1=M, 2=F)`)
individual_info_df$sex <- as.factor(individual_info_df$`sex (1=M, 2=F)`)

meta_df <- individual_info_df
for (i in 1:(length(times) * length(cell_type) - 1)){
  print(i)
  meta_df <- rbind(meta_df, individual_info_df)
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
  model.matrix( ~ cell.type + time.point + `Co-culture batch`, 
                meta_df)
# use voom to estimate weights of genes
voom_DGE <- voom(DGE, design = voom_varP_design)

# run the final Var calculation with Weights on (can be turned off by useWeights=F)
par_formula <- 
  ~ (1|Aff) + 
  (1|sex) + 
  (1|`Co-culture batch`) + 
  (1|`new iPSC line ID`) +
  (1|cell.type) + 
  (1|time.point) +
  Age

par_DGE <- 
  fitExtractVarPartModel(exprObj = voom_DGE,
                         formula = par_formula,
                         data = meta_df,
                         BPPARAM = MulticoreParam(workers = 8,
                                                  progressbar = T))

con_par_formula <- 
  ~ (Aff) + 
  (sex) + 
  (`Co-culture batch`) + 
  (`new iPSC line ID`) +
  (cell.type) + 
  (time.point) +
  Age

con_DGE <-
  canCorPairs(formula = con_par_formula,
              data = meta_df)

plotCorrMatrix(con_DGE)


# make some plots
plotVarPart(sortCols(par_DGE)) +
  ggtitle("18-line ATACseq")
