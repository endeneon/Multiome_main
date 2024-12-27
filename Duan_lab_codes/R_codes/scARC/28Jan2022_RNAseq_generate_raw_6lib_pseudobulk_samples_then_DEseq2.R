# Chuxuan Li 01/28/2022
# generate 5 cell line * 3 time point = 15 pseudobulk samples for DESeq2 analysis

# init ####
library(Seurat)
library(stringr)
library(readr)
library(DESeq2)
library(ggplot2)
library(ggrepel)

# load data ####
# Read files, separate the RNAseq data from the .h5 matrix
setwd("~/NVME/scARC_Duan_018/GRCh38_mapped_only")

h5list <- list.files(path = ".", 
                     pattern = "filtered_feature_bc_matrix.h5",
                     recursive = T,
                     include.dirs = T)
objlist <- vector(mode = "list", length = length(h5list))

for (i in 1:6){
  h5file <- Read10X_h5(filename = h5list[i])
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(string = h5list[i], 
                                                  pattern = "libraries_[0-9]_[0-6]"))
  print(paste0(i, " number of genes: ", nrow(obj), 
               ", number of cells: ", ncol(obj)))
  obj$lib.ident <- i
  obj$time.ident <- paste0(str_extract(string = h5list[i], 
                                pattern = "[0-6]"), "hr")
  objlist[[i]] <- 
    PercentageFeatureSet(obj,
                         pattern = c("^MT-"),
                         col.name = "percent.mt")
  
}

objlist[[1]]



# assign cell line identities ####
times <- c("0", "1", "6")
lines_g2 <- c("CD_27", "CD_54")
lines_g8 <- c("CD_08", "CD_25", "CD_26")
g_2_0_bc_lst <- vector(mode = "list", length = 2)
g_2_1_bc_lst <- vector(mode = "list", length = 2)
g_2_6_bc_lst <- vector(mode = "list", length = 2)
g_8_0_bc_lst <- vector(mode = "list", length = 3)
g_8_1_bc_lst <- vector(mode = "list", length = 3)
g_8_6_bc_lst <- vector(mode = "list", length = 3)

# import cell ident barcodes
for (t in times){
  for (i in 1:length(lines_g2)){
    l <- lines_g2[i]
    file_name <- paste0("/nvmefs/scARC_Duan_018/GRCh38_mapped_only/common_barcodes/g_2_",
                        t,
                        "_common_CD27_CD54.best.",
                        l)
    print(file_name)
    bc <- read_csv(file_name, col_names = FALSE)
    bc <- unlist(bc)
    bc <- str_sub(bc, end = -3L)
    if (t %in% "0"){
      g_2_0_bc_lst[[i]] <- bc
    } else if (t %in% "1"){
      g_2_1_bc_lst[[i]] <- bc
    } else {
      g_2_6_bc_lst[[i]] <- bc
    }
  }
}

for (t in times){
  for (i in 1:length(lines_g8)){
    l <- lines_g8[i]
    file_name <- paste0("/nvmefs/scARC_Duan_018/GRCh38_mapped_only/common_barcodes/g_8_",
                        t,
                        "_common_CD08_CD25_CD26.best.",
                        l)
    print(file_name)
    bc <- read_csv(file_name, col_names = FALSE)
    bc <- unlist(bc)
    bc <- str_sub(bc, end = -3L)
    if (t %in% "0"){
      g_8_0_bc_lst[[i]] <- bc
    } else if (t %in% "1"){
      g_8_1_bc_lst[[i]] <- bc
    } else {
      g_8_6_bc_lst[[i]] <- bc
    }
  }
}

barcodeslist <- list(g_2_0_bc_lst, g_2_1_bc_lst, g_2_6_bc_lst, 
                  g_8_0_bc_lst, g_8_1_bc_lst, g_8_6_bc_lst)
barcodeslist[[1]][[2]]

# assign cell line identities to the objects
for (i in 1:length(objlist)){
  print(i)
  obj <- objlist[[i]]
  obj$cell.line.ident <- NA
  if (i %in% 1:3){
    obj$cell.line.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                          barcodeslist[[i]][[1]]] <- "CD_27"
    obj$cell.line.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                          barcodeslist[[i]][[2]]] <- "CD_54"
  } else {
    obj$cell.line.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                          barcodeslist[[i]][[1]]] <- "CD_08"
    obj$cell.line.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                          barcodeslist[[i]][[2]]] <- "CD_25"
    obj$cell.line.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in% 
                          barcodeslist[[i]][[3]]] <- "CD_26"
  }
  
  print(sum(is.na(obj$cell.line.ident)))
  obj$cell.line.ident[is.na(obj$cell.line.ident)] <- "unmatched"
  print(unique(obj$cell.line.ident))

  objlist[[i]] <- subset(obj, cell.line.ident != "unmatched")
}


# assign known cell type idents ####
setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_pseudobulk")
load("../Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/renormalized_4celltype_obj.RData")

unique(filtered_obj$cell.type)
GABA <- subset(filtered_obj, subset = cell.type == "GABA")
nmglut <- subset(filtered_obj, subset = cell.type == "NEFM_neg_glut")
npglut <- subset(filtered_obj, subset = cell.type == "NEFM_pos_glut")
NPC <- subset(filtered_obj, cell.type == "NPC")

glut <- subset(filtered_obj, subset = cell.type %in% c("NEFM_pos_glut",
                                                       "NEFM_neg_glut"))

newobjlist <- objlist
for (i in 1:length(objlist)){
  print(i)
  obj <- objlist[[i]]
  obj$cell.type.ident <- "unmatched"
  obj$cell.type.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in%
                        str_sub(GABA@assays$RNA@counts@Dimnames[[2]], end = -5L)] <- "GABA"
  # obj$cell.type.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in%
  #                       str_sub(nmglut@assays$RNA@counts@Dimnames[[2]], end = -5L)] <- "NEFM- glut"
  # obj$cell.type.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in%
  #                       str_sub(npglut@assays$RNA@counts@Dimnames[[2]], end = -5L)] <- "NEFM+ glut"
  obj$cell.type.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in%
                        str_sub(NPC@assays$RNA@counts@Dimnames[[2]], end = -5L)] <- "NPC"
  obj$cell.type.ident[str_sub(obj@assays$RNA@counts@Dimnames[[2]], end = -3L) %in%
                        str_sub(glut@assays$RNA@counts@Dimnames[[2]], end = -5L)] <- "glut"
  print(unique(obj$cell.type.ident))
  newobjlist[[i]] <- obj
}

GABAlist <- vector(mode = "list", length = length(objlist))
nmglutlist <- vector(mode = "list", length = length(objlist))
npglutlist <- vector(mode = "list", length = length(objlist))
NPClist <- vector(mode = "list", length = length(objlist))

glutlist <- vector(mode = "list", length = length(objlist))

for (i in 1:length(objlist)){
  GABAlist[[i]] <- subset(newobjlist[[i]], subset = cell.type.ident == "GABA")
  #nmglutlist[[i]] <- subset(newobjlist[[i]], subset = cell.type.ident == "NEFM- glut")
  #npglutlist[[i]] <- subset(newobjlist[[i]], subset = cell.type.ident == "NEFM+ glut")
  NPClist[[i]] <- subset(newobjlist[[i]], subset = cell.type.ident == "NPC")
  glutlist[[i]] <- subset(newobjlist[[i]], subset = cell.type.ident == "glut")
}

genelist <- rownames(GABAlist[[1]]@assays$RNA@data)
filler <- rep(0, length(genelist))
GABA_pseudobulk <- data.frame(CD_27_0hr = filler, CD_27_1hr = filler, CD_27_6hr = filler,
                              CD_54_0hr = filler, CD_54_1hr = filler, CD_54_6hr = filler, 
                              CD_08_0hr = filler, CD_08_1hr = filler, CD_08_6hr = filler,
                              CD_25_0hr = filler, CD_25_1hr = filler, CD_25_6hr = filler,
                              CD_26_0hr = filler, CD_26_1hr = filler, CD_26_6hr = filler,
                              row.names = genelist)
genelist <- rownames(nmglutlist[[1]]@assays$RNA@data)
filler <- rep(0, length(genelist))
nmglut_pseudobulk <- data.frame(CD_27_0hr = filler, CD_27_1hr = filler, CD_27_6hr = filler,
                              CD_54_0hr = filler, CD_54_1hr = filler, CD_54_6hr = filler,
                              CD_08_0hr = filler, CD_08_1hr = filler, CD_08_6hr = filler,
                              CD_25_0hr = filler, CD_25_1hr = filler, CD_25_6hr = filler,
                              CD_26_0hr = filler, CD_26_1hr = filler, CD_26_6hr = filler, 
                              row.names = genelist)

genelist <- rownames(npglutlist[[1]]@assays$RNA@data)
filler <- rep(0, length(genelist))
npglut_pseudobulk <- data.frame(CD_27_0hr = filler, CD_27_1hr = filler, CD_27_6hr = filler,
                                CD_54_0hr = filler, CD_54_1hr = filler, CD_54_6hr = filler,
                                CD_08_0hr = filler, CD_08_1hr = filler, CD_08_6hr = filler,
                                CD_25_0hr = filler, CD_25_1hr = filler, CD_25_6hr = filler,
                                CD_26_0hr = filler, CD_26_1hr = filler, CD_26_6hr = filler, 
                                row.names = genelist)
genelist <- rownames(NPClist[[1]]@assays$RNA@data)
filler <- rep(0, length(genelist))
NPC_pseudobulk <- data.frame(CD_27_0hr = filler, CD_27_1hr = filler, CD_27_6hr = filler,
                              CD_54_0hr = filler, CD_54_1hr = filler, CD_54_6hr = filler, 
                              CD_08_0hr = filler, CD_08_1hr = filler, CD_08_6hr = filler,
                              CD_25_0hr = filler, CD_25_1hr = filler, CD_25_6hr = filler,
                              CD_26_0hr = filler, CD_26_1hr = filler, CD_26_6hr = filler,
                              row.names = genelist)
genelist <- rownames(glutlist[[1]]@assays$RNA@data)
filler <- rep(0, length(genelist))
glut_pseudobulk <- data.frame(CD_27_0hr = filler, CD_27_1hr = filler, CD_27_6hr = filler,
                             CD_54_0hr = filler, CD_54_1hr = filler, CD_54_6hr = filler, 
                             CD_08_0hr = filler, CD_08_1hr = filler, CD_08_6hr = filler,
                             CD_25_0hr = filler, CD_25_1hr = filler, CD_25_6hr = filler,
                             CD_26_0hr = filler, CD_26_1hr = filler, CD_26_6hr = filler,
                             row.names = genelist)
# make pseudobulk dataframe
for (i in 1:length(GABAlist)){
  if (i %in% 1:3){
    CD27 <- subset(GABAlist[[i]], subset = cell.line.ident == "CD_27")
    CD54 <- subset(GABAlist[[i]], subset = cell.line.ident == "CD_54")
    GABA_pseudobulk[, i] <- as.array(rowSums(CD27@assays$RNA@data))
    GABA_pseudobulk[, (i + 3)] <- as.array(rowSums(CD54@assays$RNA@data))
  } else {
    CD08 <- subset(GABAlist[[i]], subset = cell.line.ident == "CD_08")
    CD25 <- subset(GABAlist[[i]], subset = cell.line.ident == "CD_25")
    CD26 <- subset(GABAlist[[i]], subset = cell.line.ident == "CD_26")
    GABA_pseudobulk[, (i + 3)] <- as.array(rowSums(CD08@assays$RNA@data))
    GABA_pseudobulk[, (i + 6)] <- as.array(rowSums(CD25@assays$RNA@data))
    GABA_pseudobulk[, (i + 9)] <- as.array(rowSums(CD26@assays$RNA@data))
  }
}

for (i in 1:length(nmglutlist)){
  if (i %in% 1:3){
    CD27 <- subset(nmglutlist[[i]], subset = cell.line.ident == "CD_27")
    CD54 <- subset(nmglutlist[[i]], subset = cell.line.ident == "CD_54")
    nmglut_pseudobulk[, i] <- as.array(rowSums(CD27@assays$RNA@data))
    nmglut_pseudobulk[, (i + 3)] <- as.array(rowSums(CD54@assays$RNA@data))
  } else {
    CD08 <- subset(nmglutlist[[i]], subset = cell.line.ident == "CD_08")
    CD25 <- subset(nmglutlist[[i]], subset = cell.line.ident == "CD_25")
    CD26 <- subset(nmglutlist[[i]], subset = cell.line.ident == "CD_26")
    nmglut_pseudobulk[, (i + 3)] <- as.array(rowSums(CD08@assays$RNA@data))
    nmglut_pseudobulk[, (i + 6)] <- as.array(rowSums(CD25@assays$RNA@data))
    nmglut_pseudobulk[, (i + 9)] <- as.array(rowSums(CD26@assays$RNA@data))
  }
}
for (i in 1:length(npglutlist)){
  if (i %in% 1:3){
    CD27 <- subset(npglutlist[[i]], subset = cell.line.ident == "CD_27")
    CD54 <- subset(npglutlist[[i]], subset = cell.line.ident == "CD_54")
    npglut_pseudobulk[, i] <- as.array(rowSums(CD27@assays$RNA@data))
    npglut_pseudobulk[, (i + 3)] <- as.array(rowSums(CD54@assays$RNA@data))
  } else {
    CD08 <- subset(npglutlist[[i]], subset = cell.line.ident == "CD_08")
    CD25 <- subset(npglutlist[[i]], subset = cell.line.ident == "CD_25")
    CD26 <- subset(npglutlist[[i]], subset = cell.line.ident == "CD_26")
    npglut_pseudobulk[, (i + 3)] <- as.array(rowSums(CD08@assays$RNA@data))
    npglut_pseudobulk[, (i + 6)] <- as.array(rowSums(CD25@assays$RNA@data))
    npglut_pseudobulk[, (i + 9)] <- as.array(rowSums(CD26@assays$RNA@data))
  }
}

for (i in 1:length(NPClist)){
  if (i %in% 1:3){
    CD27 <- subset(NPClist[[i]], subset = cell.line.ident == "CD_27")
    CD54 <- subset(NPClist[[i]], subset = cell.line.ident == "CD_54")
    NPC_pseudobulk[, i] <- as.array(rowSums(CD27@assays$RNA@data))
    NPC_pseudobulk[, (i + 3)] <- as.array(rowSums(CD54@assays$RNA@data))
  } else {
    CD08 <- subset(NPClist[[i]], subset = cell.line.ident == "CD_08")
    CD25 <- subset(NPClist[[i]], subset = cell.line.ident == "CD_25")
    CD26 <- subset(NPClist[[i]], subset = cell.line.ident == "CD_26")
    NPC_pseudobulk[, (i + 3)] <- as.array(rowSums(CD08@assays$RNA@data))
    NPC_pseudobulk[, (i + 6)] <- as.array(rowSums(CD25@assays$RNA@data))
    NPC_pseudobulk[, (i + 9)] <- as.array(rowSums(CD26@assays$RNA@data))
  }
}

for (i in 1:length(glutlist)){
  if (i %in% 1:3){
    CD27 <- subset(glutlist[[i]], subset = cell.line.ident == "CD_27")
    CD54 <- subset(glutlist[[i]], subset = cell.line.ident == "CD_54")
    glut_pseudobulk[, i] <- as.array(rowSums(CD27@assays$RNA@data))
    glut_pseudobulk[, (i + 3)] <- as.array(rowSums(CD54@assays$RNA@data))
  } else {
    CD08 <- subset(glutlist[[i]], subset = cell.line.ident == "CD_08")
    CD25 <- subset(glutlist[[i]], subset = cell.line.ident == "CD_25")
    CD26 <- subset(glutlist[[i]], subset = cell.line.ident == "CD_26")
    glut_pseudobulk[, (i + 3)] <- as.array(rowSums(CD08@assays$RNA@data))
    glut_pseudobulk[, (i + 6)] <- as.array(rowSums(CD25@assays$RNA@data))
    glut_pseudobulk[, (i + 9)] <- as.array(rowSums(CD26@assays$RNA@data))
  }
}

save("GABA_pseudobulk", "nmglut_pseudobulk",
     "npglut_pseudobulk", "NPC_pseudobulk", 
     file = "data_frames_for_pseudobulk_raw.RData")


# DESeq2 ####
library("DESeq2")

# df_list <- list(GABA_pseudobulk, nmglut_pseudobulk,
#                 npglut_pseudobulk, NPC_pseudobulk
#                 )
df_list <- list(GABA_pseudobulk, glut_pseudobulk, NPC_pseudobulk)

# types <- c("GABA", "nmglut",
#            "npglut", "NPC")
types <- c("GABA", "glut", "NPC")

names(df_list) <- types

res_0v1_list <- vector(mode = "list", length = length(df_list))
names(res_0v1_list) <- types
res_0v6_list <- vector(mode = "list", length = length(df_list))
names(res_0v6_list) <- types

for (i in 1:length(df_list)){
  # make colData df
  print(names(df_list[[i]]))
  coldata <- data.frame(condition = str_extract(string = colnames(df_list[[i]]),
                                                     pattern = "[0-6]hr"),
                             type = str_extract(string = colnames(df_list[[i]]),
                                                pattern = "CD_[0-9][0-9]"))
  rownames(coldata) <- colnames(df_list[[i]])
  
  # check
  print(all(rownames(coldata) %in% colnames(df_list[[i]]))) # TRUE
  print(all(rownames(coldata) == colnames(df_list[[i]]))) #TRUE
  
  # make deseq2 obj
  print(i)
  dds <- DESeqDataSetFromMatrix(countData = df_list[[i]],
                                     colData = coldata,
                                     design = ~ type + condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "1hr", "0hr"))
  res <- as.data.frame(res)
  res <- res[res$baseMean > 10, ]
  res$q_value <- p.adjust(p = res$pvalue, method = "BH")
  res$gene <- rownames(res)
  res_0v1_list[[i]] <- res
  res <- results(dds, contrast = c("condition", "6hr", "0hr"))
  res <- as.data.frame(res)
  res <- res[res$baseMean > 10, ]
  res$q_value <- p.adjust(p = res$pvalue, method = "BH")
  res$gene <- rownames(res)
  res_0v6_list[[i]] <- res
}

# count up/down regulated genes and dynamic genes ####
res_0v1_up_list <- vector("list", length = length(res_0v1_list))
res_0v1_do_list <- vector("list", length = length(res_0v1_list))
res_0v6_up_list <- vector("list", length = length(res_0v6_list))
res_0v6_do_list <- vector("list", length = length(res_0v6_list))

for (i in 1:length(res_0v1_list)) {
  df <- res_0v1_list[[i]]
  print("1v0 up")
  res_0v1_up_list[[i]] <- df[df$log2FoldChange > 0 & df$q_value < 0.05, ]
  print(nrow(res_0v1_up_list[[i]]))
  print("1v0 down")
  res_0v1_do_list[[i]] <- df[df$log2FoldChange < 0 & df$q_value < 0.05, ]
  print(nrow(res_0v1_do_list[[i]]))
  df <- res_0v6_list[[i]]
  print("6v0 up")
  res_0v6_up_list[[i]] <- df[df$log2FoldChange > 0 & df$q_value < 0.05, ]
  print(nrow(res_0v6_up_list[[i]]))
  print("6v0 down")
  res_0v6_do_list[[i]] <- df[df$log2FoldChange < 0 & df$q_value < 0.05, ]
  print(nrow(res_0v6_do_list[[i]]))
}

intersect_0v1_up <- intersect(res_0v1_up_list[[1]]$gene, res_0v1_up_list[[2]]$gene)
intersect_0v1_do <- intersect(res_0v1_do_list[[1]]$gene, res_0v1_do_list[[2]]$gene)
intersect_0v6_up <- intersect(res_0v6_up_list[[1]]$gene, res_0v6_up_list[[2]]$gene)
intersect_0v6_do <- intersect(res_0v6_do_list[[1]]$gene, res_0v6_do_list[[2]]$gene)
GABAonly_0v1_up <- res_0v1_up_list[[1]]$gene[!res_0v1_up_list[[1]]$gene%in% intersect_0v1_up]
GABAonly_0v1_do <- res_0v1_do_list[[1]]$gene[!res_0v1_do_list[[1]]$gene%in% intersect_0v1_do]
GABAonly_0v6_up <- res_0v6_up_list[[1]]$gene[!res_0v6_up_list[[1]]$gene%in% intersect_0v6_up]
GABAonly_0v6_do <- res_0v6_do_list[[1]]$gene[!res_0v6_do_list[[1]]$gene%in% intersect_0v6_do]
glutonly_0v1_up <- res_0v1_up_list[[2]]$gene[!res_0v1_up_list[[2]]$gene%in% intersect_0v1_up]
glutonly_0v1_do <- res_0v1_do_list[[2]]$gene[!res_0v1_do_list[[2]]$gene%in% intersect_0v1_do]
glutonly_0v6_up <- res_0v6_up_list[[2]]$gene[!res_0v6_up_list[[2]]$gene%in% intersect_0v6_up]
glutonly_0v6_do <- res_0v6_do_list[[2]]$gene[!res_0v6_do_list[[2]]$gene%in% intersect_0v6_do]

names <- c("GABAonly_0v1_up", "GABAonly_0v1_do", "GABAonly_0v6_up", "GABAonly_0v6_do",
"glutonly_0v1_up", "glutonly_0v1_do", "glutonly_0v6_up", "glutonly_0v6_do",
"intersect_0v1_up", "intersect_0v1_do", "intersect_0v6_up", "intersect_0v6_do")
genelist <- list(GABAonly_0v1_up, GABAonly_0v1_do, GABAonly_0v6_up, GABAonly_0v6_do,
                 glutonly_0v1_up, glutonly_0v1_do, glutonly_0v6_up, glutonly_0v6_do,
                 intersect_0v1_up, intersect_0v1_do, intersect_0v6_up, intersect_0v6_do)
for (i in 1:length(genelist)) {
  lst <- genelist[[i]]
  print(length(lst))
  write.table(lst, file = paste0(names[i], "_genes.txt"), quote = F, sep = ",", row.names = F)
}
#save(res_0v1_list, res_0v6_list, file = "full_pseudobulk_DEG_lists_without_any_filtering.RData")


for (i in 1:length(res_0v1_list)){
  filename <- paste0(types[i], "_1v0_full_DEG_list_after_basemean_filtering_only.csv")
  print(filename)
  write.table(res_0v1_list[[i]], file = filename, quote = F, sep = ",", col.names = T)
  filename <- paste0(types[i], "_6v0_full_DEG_list_after_basemean_filtering_only.csv")
  print(filename)
  write.table(res_0v6_list[[i]], file = filename, quote = F, sep = ",", col.names = T)
  
  #filtered list
  filename <- paste0(types[i], "_1v0_upregulated_significant_DEG.csv")
  print(filename)
  write.table(res_0v1_list[[i]][res_0v1_list[[i]]$q_value < 0.05 & res_0v1_list[[i]]$log2FoldChange > 0, ],
              file = filename, quote = F, sep = ",", col.names = T)
  filename <- paste0(types[i], "_1v0_downregulated_significant_DEG.csv")
  print(filename)
  write.table(res_0v1_list[[i]][res_0v1_list[[i]]$q_value < 0.05 & res_0v1_list[[i]]$log2FoldChange < 0, ],
              file = filename, quote = F, sep = ",", col.names = T)
  #6v0
  filename <- paste0(types[i], "_6v0_upregulated_significant_DEG.csv")
  print(filename)
  write.table(res_0v6_list[[i]][res_0v6_list[[i]]$q_value < 0.05 & res_0v6_list[[i]]$log2FoldChange > 0, ],
              file = filename, quote = F, sep = ",", col.names = T)
  filename <- paste0(types[i], "_6v0_downregulated_significant_DEG.csv")
  print(filename)
  write.table(res_0v6_list[[i]][res_0v6_list[[i]]$q_value < 0.05 & res_0v6_list[[i]]$log2FoldChange < 0, ],
              file = filename, quote = F, sep = ",", col.names = T)}
hist(res_0v1_list[[1]]$baseMean, breaks = 10000, xlim = c(0, 10000), 
     main = "histogram of GABA 0v1 basemean distribution")

# look at specific genes
res_0v1_list[[1]][c("FOS", "NPAS4", "VGF", "BDNF"), ]
res_0v6_list[[1]][c("FOS", "NPAS4", "VGF", "BDNF"), ]

# plot bar graphs
early_list <- c("FOS", "FOSB", "NPAS4", "NR4A1")
late_list <- c("BDNF", "IGF1", "VGF")
df_to_plot <- rbind(res_0v1_list[[1]][early_list, ], res_0v6_list[[1]][early_list, ],
                    res_0v1_list[[2]][early_list, ], res_0v6_list[[2]][early_list, ],
                    res_0v1_list[[3]][early_list, ], res_0v6_list[[3]][early_list, ]#,
                    #res_0v1_list[[4]][early_list, ], res_0v6_list[[4]][early_list, ]
                    )
df_to_plot$time <- c(rep_len("0v1", length(early_list)), rep_len("0v6", length(early_list)), 
                     rep_len("0v1", length(early_list)), rep_len("0v6", length(early_list)), 
                     rep_len("0v1", length(early_list)), rep_len("0v6", length(early_list))#, 
                     #rep_len("0v1", length(early_list)), rep_len("0v6", length(early_list))
                     )                    
df_to_plot$cell.type <- c(rep_len("GABA", 2*length(early_list)), 
                          rep_len("NEFM- glut", 2*length(early_list)), 
                          rep_len("NEFM+ glut", 2*length(early_list))#, 
                          #rep_len("NPC", 2*length(early_list))
                          )
#df_to_plot$gene.name <- rownames(df_to_plot)
df_to_plot <- as.data.frame(df_to_plot)

# for early, mark if they are excitatory or inhibitory
gene_labeller <- as_labeller(c('BDNF' = "BDNF (excitatory)", 
                               'IGF1' = "IGF1 (inhibitory)", 
                               'VGF' = "VGF (shared)"))
ggplot(df_to_plot, aes(x = cell.type, 
                       y = log2FoldChange,
                       color = cell.type,
                       fill = time,
                       group = time,
                       ymax = log2FoldChange - 1/2*lfcSE, 
                       ymin = log2FoldChange + 1/2*lfcSE)) + 
  xlab("") +
  ylab("log2(FC)") +
  geom_col(position = position_dodge(0.6),
           width = 0.5,
           color = "black") +
  geom_errorbar(color = "black",
                position = position_dodge(0.6),
                width = 0.5) +
  facet_grid(cols = vars(gene)#,labeller = gene_labeller
             ) +
  scale_fill_manual(values = c("indianred", "steelblue")) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("differential expression of early response genes")

# edgeR ####

library(edgeR)

Glut_DGEList <- DGEList(counts = as.matrix(glut_pseudobulk),
                        genes = rownames(glut_pseudobulk), 
                        samples = colnames(glut_pseudobulk),
                        remove.zeros = T)

Glut_DGEList <- calcNormFactors(Glut_DGEList)
Glut_DGEList <- estimateDisp(Glut_DGEList)

hist(cpm(Glut_DGEList), breaks = 10000, xlim = c(0, 10))
sum(cpm(Glut_DGEList) < 1)

# filter DGEList
Glut_DGEList <- Glut_DGEList[rowSums(cpm(Glut_DGEList) > 1) >= 5, ,
                             keep.lib.sizes = F]
# Glut_DGEList$genes

Glut_design_matrix <- coldata_glut
Glut_design_matrix$condition <- as.factor(Glut_design_matrix$condition)
Glut_design_matrix$type <- as.factor(Glut_design_matrix$type)

Glut_design <- model.matrix(~ 0 + 
                              Glut_design_matrix$condition +
                              Glut_design_matrix$type)

Glut_fit <- glmQLFit(Glut_DGEList,
                     design = Glut_design)
Glut_1vs0 <- glmQLFTest(Glut_fit, 
                        contrast = c(-1, 1, 0, 0, 0, 0, 0))
Glut_6vs0 <- glmQLFTest(Glut_fit, 
                        contrast = c(-1, 0, 1, 0, 0, 0, 0))
Glut_6vs1 <- glmQLFTest(Glut_fit, 
                        contrast = c(0, -1, 1, 0, 0, 0, 0))

Glut_1vs0 <- Glut_1vs0$table
Glut_6vs0 <- Glut_6vs0$table
Glut_6vs1 <- Glut_6vs1$table

genes_to_check <- c("BDNF", "FOS", "GAPDH", "PNOC", "NRN1", "VGF")

Glut_1vs0[rownames(Glut_1vs0) %in% genes_to_check, ]
Glut_6vs0[rownames(Glut_6vs0) %in% genes_to_check, ]
Glut_6vs1[rownames(Glut_6vs1) %in% genes_to_check, ]



GABA_DGEList <- DGEList(counts = as.matrix(GABA_pseudobulk),
                        genes = rownames(GABA_pseudobulk), 
                        samples = colnames(GABA_pseudobulk),
                        remove.zeros = T)

GABA_DGEList <- calcNormFactors(GABA_DGEList)
GABA_DGEList <- estimateDisp(GABA_DGEList)

hist(cpm(GABA_DGEList), breaks = 10000, xlim = c(0, 10))
sum(cpm(GABA_DGEList) < 1)

# filter DGEList
GABA_DGEList <- GABA_DGEList[rowSums(cpm(GABA_DGEList) > 1) >= 5, ,
                             keep.lib.sizes = F]
# GABA_DGEList$genes

GABA_design_matrix <- coldata_GABA
GABA_design_matrix$condition <- as.factor(GABA_design_matrix$condition)
GABA_design_matrix$type <- as.factor(GABA_design_matrix$type)

GABA_design <- model.matrix(~ 0 + 
                              GABA_design_matrix$condition +
                              GABA_design_matrix$type)

GABA_fit <- glmQLFit(GABA_DGEList,
                     design = GABA_design)
GABA_1vs0 <- glmQLFTest(GABA_fit, 
                        contrast = c(-1, 1, 0, 0, 0, 0, 0))
GABA_6vs0 <- glmQLFTest(GABA_fit, 
                        contrast = c(-1, 0, 1, 0, 0, 0, 0))
GABA_6vs1 <- glmQLFTest(GABA_fit, 
                        contrast = c(0, -1, 1, 0, 0, 0, 0))

GABA_1vs0 <- GABA_1vs0$table
GABA_6vs0 <- GABA_6vs0$table
GABA_6vs1 <- GABA_6vs1$table

genes_to_check <- c("BDNF", "FOS", "GAPDH", "PNOC", "NRN1", "VGF")

GABA_1vs0[rownames(GABA_1vs0) %in% genes_to_check, ]
GABA_6vs0[rownames(GABA_6vs0) %in% genes_to_check, ]
GABA_6vs1[rownames(GABA_6vs1) %in% genes_to_check, ]


# volcano plots ####
namelist1 <- c("GABA_0v1", "nmglut_0v1", "npglut_0v1", "NPC_0v1")
namelist6 <- c("GABA_0v6", "nmglut_0v6", "npglut_0v6", "NPC_0v6")
for (i in 1:length(res_0v1_list)){
  res <- res_0v1_list[[i]]
  res$gene.symbol <- rownames(res)
  res$significance <- "nonsignificant"
  res$significance[res$q_value < 0.05 & res$log2FoldChange > 0] <- "up"
  res$significance[res$q_value < 0.05 & res$log2FoldChange < 0] <- "down"
  unique(res$significance)
  res$significance <- factor(res$significance, levels = c("up", "nonsignificant", "down"))
  res$neg_log_pval <- (0 - log2(res$pvalue))
  res$labelling <- ""
  for (j in c("FOS", "NPAS4", "IGF1", "NR4A1","FOSB", "EGR1", "BDNF")){
    res$labelling[res$gene.symbol %in% j] <- j
  }
  print(namelist1[i])
  pdf(paste0(namelist1[i], "_volcano_plot.pdf"))
  p <- ggplot(data = as.data.frame(res),
              aes(x = log2FoldChange, 
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
    xlim(c(-10, 10)) +
    ggtitle(namelist1[i])
  print(p)
  dev.off()
}

