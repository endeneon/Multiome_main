# Chuxuan Li 01/30/2023
# DE analysis comparing case and control using 018, 022, 024, 025, 029 data,
#separated by time points (so case and control will be contrasted in each of the
#3 time points)

# init ####
library(Seurat)
library(limma)
library(edgeR)
library(sva)

library(readr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(readxl)

load("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/029_RNA_integrated_labeled.RData")
obj_029 <- integrated_labeled
rm(integrated_labeled)
load("/nvmefs/scARC_Duan_018/Duan_project_025_RNA/Analysis_part2_GRCh38/integrated_labeled.RData")
obj_025 <- integrated_labeled
rm(integrated_labeled)
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000_with_covars_from_30Sep2022_single_cell_DE.RData")
obj_024 <- obj
rm(obj)
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/renormalized_4celltype_obj.RData")
obj_018 <- filtered_obj
rm(filtered_obj)
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/20line_codes_and_rdata/demux_20line_integrated_labeled_obj.RData")
obj_022 <- integrated_labeled
rm(integrated_labeled)

# clean data ####
load("redone_lines_char_vector.RData")
obj_022$cell.line.ident[obj_022$cell.line.ident %in% redone_lines] <- "redone"
unique(obj_022$cell.line.ident)
obj_022 <- subset(obj_022, cell.line.ident != "redone")

objlist <- list(obj_018, obj_022, obj_024, obj_025, obj_029)

genes.use <- rownames(obj_024@assays$RNA)
for (i in 1:length(objlist)) {
  counts <- GetAssayData(objlist[[i]], assay = "RNA")
  counts <- counts[which(rownames(counts) %in% genes.use), ]
  o <- subset(objlist[[i]], features = rownames(counts))
  objlist[[i]] <- o
}


for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  print(unique(obj$cell.type))
  DefaultAssay(obj) <- "RNA"
  types <- unique(objlist[[i]]$cell.type)
  if ("NEFM_pos_glut" %in% types) {
    obj$cell.type <- str_replace(obj$cell.type, "NEFM_pos_glut", "npglut")
    obj$cell.type <- str_replace(obj$cell.type, "NEFM_neg_glut", "nmglut")
  }
  if ("SEMA3E_pos_GABA" %in% types | "SST_pos_GABA" %in% types) {
    obj$cell.type[str_detect(obj$cell.type, "SEMA3E_pos_GABA")] <- "GABA"
    obj$cell.type[str_detect(obj$cell.type, "SST_pos_GABA")] <- "GABA"
  }
  objlist[[i]] <- obj
}

# assign time x cell line ident
for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  lines <- unique(obj$cell.line.ident)
  times <- unique(obj$time.ident)
  DefaultAssay(obj) <- "RNA"
  obj$timexline.ident <- ""
  for (l in lines) {
    print(l)
    for (t in times) {
      print(t)
      obj$timexline.ident[obj$cell.line.ident == l & obj$time.ident == t] <- 
        paste0(l, "_", t)
    }
  }
  objlist[[i]] <- obj
}

names(objlist) <- c("obj_018", "obj_022", "obj_024", "obj_025", "obj_029")
rm(obj_018, obj_022, obj_024, obj_025, obj_029)
gc()

# make count matrices ####
# total 3 matrices for 3 cell types, each matix has 5+18+20=43 lines x 3 time
# points = 129 samples/columns.
makeMat4CombatseqLine <- function(typeobj) {
  linetimes <- sort(unique(typeobj$timexline.ident))
  for (i in 1:length(linetimes)) {
    print(linetimes[i])
    obj <- subset(typeobj, timexline.ident == linetimes[i])
    if (i == 1) {
      mat <- as.array(rowSums(obj@assays$RNA@counts))
    } else {
      to_bind <- as.array(rowSums(obj@assays$RNA@counts))
      mat <- cbind(mat, to_bind)
    }
  }
  colnames(mat) <- linetimes
  return(mat)
}
GABA_list <- vector("list", length(objlist))
nmglut_list <- vector("list", length(objlist))
npglut_list <- vector("list", length(objlist))

for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  GABA_list[[i]] <- subset(obj, cell.type == "GABA")
  nmglut_list[[i]] <- subset(obj, cell.type == "nmglut")
  npglut_list[[i]] <- subset(obj, cell.type == "npglut")
}
type_obj_lsts <- list(GABA_list, nmglut_list, npglut_list)
names(type_obj_lsts) <- c("GABA", "nmglut", "npglut")

mat_lst <- vector("list", 3)
for (i in 1:length(type_obj_lsts)) {
  print(names(type_obj_lsts)[i])
  objlst <- type_obj_lsts[[i]]
  for (k in 1:length(objlst)) {
    obj <- objlst[[k]]
    DefaultAssay(obj) <- "RNA"
    if (k == 1) {
      mat <- makeMat4CombatseqLine(obj)
    } else {
      print("k > 1")
      mat <- cbind(mat, makeMat4CombatseqLine(obj))
    }
  }
  mat <- mat[, order(colnames(mat))]
  mat_lst[[i]] <- mat
}
names(mat_lst) <- names(type_obj_lsts)

mat_lst_bytime <- vector("list", length = 3 * length(mat_lst))
names(mat_lst_bytime) <- rep(names(mat_lst), each = 3)
for (i in 1:length(mat_lst)) {
  type <- names(mat_lst)[i]
  mat <- mat_lst[[i]]
  mat_lst_bytime[[i * 3 - 2]] <- mat[, str_detect(colnames(mat), "0hr")]
  names(mat_lst_bytime)[i * 3 - 2] <- paste0(names(mat_lst_bytime)[i * 3 - 2], "_0hr")
  mat_lst_bytime[[i * 3 - 1]] <- mat[, str_detect(colnames(mat), "1hr")]
  names(mat_lst_bytime)[i * 3 - 1] <- paste0(names(mat_lst_bytime)[i * 3 - 1], "_1hr")
  mat_lst_bytime[[i * 3]] <- mat[, str_detect(colnames(mat), "6hr")]
  names(mat_lst_bytime)[i * 3] <- paste0(names(mat_lst_bytime)[i * 3], "_6hr")
}
celllines <- str_remove(colnames(mat_lst_bytime[[1]]), "_[0|1|6]hr$")

for (i in 1:length(celllines)) {
  print(celllines[i])
  if (i == 1) {
    batch = unique(covar_table_final$batch[covar_table_final$cell_line == celllines[i]])
  } else {
    batch = c(batch, unique(covar_table_final$batch[covar_table_final$cell_line == celllines[i]]))
  }
}
batch <- as.factor(batch)

adj_mat_lst <- vector("list", length(mat_lst_bytime)) 
for (i in 1:length(mat_lst_bytime)) {
  adj_mat_lst[[i]] <- ComBat_seq(mat_lst_bytime[[i]], batch = batch)
}
names(adj_mat_lst) <- names(mat_lst_bytime)

save(mat_lst, file = "018-029_by_type_linextime_mat_lst.RData")
save(adj_mat_lst, file = "018-029_by_type_and_time_combat_adj_mat_lst.RData")
save(mat_lst_bytime, file = "018-029_by_type_and_time_mat_lst_bytime.RData")


# make design matrices ####
covar_table_final <- covar_table_final[order(covar_table_final$cell_line), ]
covar_table_final$batch <- factor(covar_table_final$batch)
covar_table_final$sex <- factor(covar_table_final$sex)
covar_table_final$cell_line <- factor(covar_table_final$cell_line)
covar_table_final$time <- factor(covar_table_final$time)
covar_table_final$aff <- factor(covar_table_final$aff)
covar_table_final$age_sq <- (covar_table_final$age * covar_table_final$age)
covar_table_0hr <- covar_table_final[covar_table_final$time == "0hr", ]
covar_table_1hr <- covar_table_final[covar_table_final$time == "1hr", ]
covar_table_6hr <- covar_table_final[covar_table_final$time == "6hr", ]

designs <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  type <- str_split(names(adj_mat_lst)[i], pattern = "_", n = 2, simplify = T)[,1]
  time <- str_split(names(adj_mat_lst)[i], pattern = "_", n = 2, simplify = T)[,2]
  if (time == "0hr") {
    d = covar_table_0hr
  } else if (time == "1hr") {
    d = covar_table_1hr
  } else {
    d = covar_table_6hr
  }
  if (type == "GABA") {
    f = d$GABA_fraction
  } else if (type == "nmglut") {
    f = d$nmglut_fraction
  } else {
    f = d$npglut_fraction
  }
  designs[[i]] <- model.matrix(~0 + aff + age_sq + sex + batch + f, 
                               data = d)
}
names(designs) <- names(adj_mat_lst)

# do limma voom ####
createDGE <- function(count_matrix){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = celllines)
  A <- rowSums(y$counts)
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}
filterByCpm <- function(df, cutoff, proportion = 1, index1, index2) {
  cpm <- cpm(df)
  group1 <- cpm[, index1]
  group2 <- cpm[, index2]
  passfilter <- (rowSums(group1 >= cutoff) >= ncol(group1) * proportion |
                   rowSums(group2 >= cutoff) >= ncol(group2) * proportion) # number of samples from either group > ns
  return(passfilter)
}
cnfV <- function(y, design) {
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = T)
  return(v)  
}
contrastFit <- function(fit, design) {
  contr.matrix <- makeContrasts(
    casevscontrol = affcase-affcontrol,
    levels = colnames(design))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit <- eBayes(fit, trend = T)
  print(summary(decideTests(fit)))
  return(fit)
}
designs[[1]]

DGEs <- vector("list", length(adj_mat_lst))
voomedDGEs <- vector("list", length(adj_mat_lst))
fits <- vector("list", length(adj_mat_lst))
res_tables <- vector("list", length(fits))

DEG_summary_table <- matrix(nrow = 9, ncol = 7, dimnames = list(names(adj_mat_lst),
                                         c("upregulated (FDR < 0.05)",
                                           "downregulated (FDR < 0.05)",
                                           "nonsignificant (FDR > 0.05)",
                                           "upregulated (p-value < 0.05)",
                                           "downregulated (p-value < 0.05)",
                                           "total number of genes passed filter",
                                           "total number of genes")))
for (i in 1:length(adj_mat_lst)) {
  DGEs[[i]] <- createDGE(adj_mat_lst[[i]])
  case_ind <- covar_table_0hr$aff == "case"
  ctrl_ind <- covar_table_0hr$aff == "control"
  ind.keep <- filterByCpm(df = DGEs[[i]], cutoff = 1,
                          index1 = case_ind, index2 = ctrl_ind)
  
  dge_filtered <- DGEs[[i]][ind.keep, ]
  voomedDGEs[[i]] <- cnfV(dge_filtered, designs[[i]])
  fit <- lmFit(voomedDGEs[[i]], designs[[i]])
  fits[[i]] <- contrastFit(fit, designs[[i]])
  res <- topTable(fits[[i]], p.value = Inf, number = Inf, sort.by = "P")
  res_sig <- topTable(fits[[i]], p.value = 0.05, number = Inf)
  DEG_summary_table[i, 1] <- sum(res_sig$logFC > 0)
  DEG_summary_table[i, 2] <- sum(res_sig$logFC < 0)
  DEG_summary_table[i, 3] <- sum(res$adj.P.Val > 0.05)
  DEG_summary_table[i, 4] <- sum(res$P.Value < 0.05 & res$logFC > 0)
  DEG_summary_table[i, 5] <- sum(res$P.Value < 0.05 & res$logFC < 0)
  DEG_summary_table[i, 6] <- sum(ind.keep)
  DEG_summary_table[i, 7] <- nrow(DGEs[[i]])
  res_tables[[i]] <- res
}
names(fits) <- names(adj_mat_lst)
names(res_tables) <- names(adj_mat_lst)

write.table(DEG_summary_table, 
            file = "018-029_case_v_ctrl_combat_limma_agesq_normalize_3times_together_DEG_counts.csv",
            sep = ",", quote = F)

# save DEG lists ####
save(fits, file = "018-029_case_v_control_DE_limma_fits.RData")

# saveCSV <- function(res, path, dfname) {
#   write.table(res[res$logFC > 0 & res$adj.P.Val < 0.05, ],
#               file = paste0(path, dfname, "_upregulated_genes.csv"),
#               quote = F, sep = ",", row.names = F, col.names = T)
#   write.table(res[res$logFC < 0 & res$adj.P.Val < 0.05, ],
#               file = paste0(path,  dfname, "_downregulated_genes.csv"),
#               quote = F, sep = ",", row.names = F, col.names = T)
#   write.table(res, file = paste0(path, dfname, "_all_genes.csv"),
#               quote = F, sep = ",", row.names = F, col.names = T)
# }
# path <- "./018-029_combat_limma_agesq_normalize_all_times_together_DEG_results_csv/"
# dfname <- "018-029_case_v_control_DE_"
saveCSV <- function(res, path, dfname) {
  write.table(res[res$logFC > 0 & res$P.Value < 0.05, ],
              file = paste0(path, dfname, "_upregulated_genes.csv"),
              quote = F, sep = ",", row.names = F, col.names = T)
  write.table(res[res$logFC < 0 & res$P.Value < 0.05, ],
              file = paste0(path,  dfname, "_downregulated_genes.csv"),
              quote = F, sep = ",", row.names = F, col.names = T)
  write.table(res, file = paste0(path, dfname, "_all_genes.csv"),
              quote = F, sep = ",", row.names = F, col.names = T)
}
path <- "./018-029_combat_limma_agesq_normalize_all_times_sep_DEG_results_csv/"
dfname <- "018-029_case_v_control_DE_3times_sep_"
for (i in 1:length(fits)) {
  res <- res_tables[[i]]
  n <- paste0(dfname, names(res_tables[i]))
  saveCSV(res, path, n)
}

# # plot results ####
# deg_df <- NULL
# for (i in 1:length(fits)) {
#   print(names(fits)[i])
#   df <- topTable(fits[[i]], p.value = 0.05, number = Inf)
#   if (nrow(df) != 0) {
#     df$typextime <- names(adj_mat_lst)[i]
#     if (is.null(deg_df)) {
#       deg_df <- df
#     } else (
#       deg_df <- rbind(deg_df, df)
#     )
#   } 
# }
# deg_df$neg_log_p <- (-1) * log2(deg_df$adj.P.Val)
# ggplot(deg_df, aes(x = genes, y = typextime, color = logFC, size = neg_log_p)) +
#   geom_point() +
#   scale_color_gradient2(low = "steelblue", high = "darkred") +
#   ylab("Cell type - Time point") +
#   labs(size = "-log(adjusted p-val)") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
