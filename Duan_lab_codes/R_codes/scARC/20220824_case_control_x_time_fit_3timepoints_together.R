# Chuxuan Li 08/24/2022
# DE analysis comparing case and control using 5+18+20 line data, put 3 time 
#points into one matrix to do linear model fitting to account for the effect of
#time point in the original fitted values, then extract the fitted values and 
#separate them into 0, 1, 6hrs. Finally subtract the 1/6hr from 0hr fitted values
#and fit linear model again, this time with aff + age + sex + group + fraction 
#at 0hr + fraction at 1hr, then contrast case and control.

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

load("../../../covariates_pooled_together_for_all_linextime_5+18+20.RData")
load("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/renormalized_4celltype_obj.RData")
RNA_5line <- filtered_obj
rm(filtered_obj)
load("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/labeled_nfeat3000.RData")
RNA_18line <- integrated_labeled
rm(integrated_labeled)
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/20line_codes_and_rdata/demux_20line_integrated_labeled_obj.RData")
RNA_20line <- integrated_labeled
rm(integrated_labeled)

# clean objects ####
DefaultAssay(RNA_5line) <- "RNA"
DefaultAssay(RNA_18line) <- "RNA"
DefaultAssay(RNA_20line) <- "RNA"
genes.use <- rownames(RNA_18line@assays$RNA)
counts <- GetAssayData(RNA_20line, assay = "RNA")
counts <- counts[which(rownames(counts) %in% genes.use), ]
RNA_20line <- subset(RNA_20line, features = rownames(counts))
counts <- GetAssayData(RNA_5line, assay = "RNA")
counts <- counts[which(rownames(counts) %in% genes.use), ]
RNA_5line <- subset(RNA_5line, features = rownames(counts))

unique(RNA_18line$cell.type)
RNA_18line$cell.type[RNA_18line$cell.type %in% c("SST_pos_GABA", "GABA", "SEMA3E_pos_GABA")] <- 'GABA'
RNA_18line <- subset(RNA_18line, cell.type != "unknown")
RNA_18line$cell.type <- factor(RNA_18line$cell.type, levels = c("NEFM_pos_glut", 
                                                                "GABA", "NEFM_neg_glut"))
unique(RNA_18line$cell.type)

unique(RNA_20line$cell.type)
RNA_20line$cell.type[RNA_20line$cell.type %in% c("SEMA3E_pos_GABA", "GABA")] <- 'GABA'
RNA_20line <- subset(RNA_20line, cell.type != "unknown")
RNA_20line$cell.type <- factor(RNA_20line$cell.type, levels = c("NEFM_pos_glut", 
                                                                "GABA", "NEFM_neg_glut"))
unique(RNA_20line$cell.type)

unique(RNA_5line$cell.type)
RNA_5line <- subset(RNA_5line, cell.type != "NPC")
RNA_5line$cell.type <- factor(RNA_5line$cell.type, levels = c("NEFM_pos_glut", 
                                                              "GABA", "NEFM_neg_glut"))
unique(RNA_5line$cell.type)

celltypes <- sort(unique(RNA_5line$cell.type))
times <- c("0hr", "1hr", "6hr")

orig_objs <- list(RNA_5line, RNA_18line, RNA_20line)
for (i in 1:length(orig_objs)) {
  obj <- orig_objs[[i]]
  lines <- unique(obj$cell.line.ident)
  print(lines)
  obj$timexline.ident <- ""
  for (l in lines) {
    for (t in times) {
      obj$timexline.ident[obj$cell.line.ident == l & obj$time.ident == t] <- 
        paste(l, t, sep = "_")
    }
  }
  print(unique(obj$timexline.ident))
  orig_objs[[i]] <- obj
}

type_obj_lsts <- vector("list", 3)
for (i in 1:length(celltypes)) {
  t <- celltypes[i]
  print(t)
  sublist <- vector("list", 3)
  for (j in 1:length(orig_objs)) {
    sublist[[j]] <- subset(orig_objs[[j]], cell.type == t)
  }
  type_obj_lsts[[i]] <- sublist
  names(type_obj_lsts[[i]]) <- c("5", "18", "20")
}
names(type_obj_lsts) <- celltypes

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
mat_lst <- vector("list", 3)
for (i in 1:length(type_obj_lsts)) {
  print(names(type_obj_lsts)[i])
  objlst <- type_obj_lsts[[i]]
    for (k in 1:length(objlst)) {
      obj <- objlst[[k]]
      if (k == 1) {
        mat <- makeMat4CombatseqLine(obj)
      } else {
        mat <- cbind(mat, makeMat4CombatseqLine(obj))
      }
    }
  mat <- mat[, order(colnames(mat))]
  mat_lst[[i]] <- mat
  names(mat_lst)[i] <- names(type_obj_lsts)[i]
}
celllines <- str_remove(colnames(mat_lst[[1]]), "_[0|1|6]hr$")

for (i in 1:length(celllines)) {
  print(celllines[i])
  if (i == 1) {
    batch = unique(covar_table$group[covar_table$cell_line == celllines[i]])
  } else {
    batch = c(batch, unique(covar_table$group[covar_table$cell_line == celllines[i]]))
  }
}
batch <- as.factor(batch)

adj_mat_lst <- vector("list", length(mat_lst)) 
for (i in 1:length(mat_lst)) {
  adj_mat_lst[[i]] <- ComBat_seq(mat_lst[[i]], batch = batch)
}
names(adj_mat_lst) <- names(mat_lst)
save(adj_mat_lst, file = "5+18+20line_sepby_type_colby_linxetime_adj_mat_list.RData")


# auxillary functions ####
createDGE <- function(count_matrix){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = celllines)
  A <- rowSums(y$counts)
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}
filterByCpm <- function(y, cutoff, nsample) {
  cpm <- cpm(y)
  colInd0 <- seq(1, ncol(y), by = 3)
  colInd1 <- seq(2, ncol(y), by = 3)
  colInd6 <- seq(3, ncol(y), by = 3)
  
  cpm_0hr <- cpm[, colInd0]
  cpm_1hr <- cpm[, colInd1]
  cpm_6hr <- cpm[, colInd6]
  #passfilter <- rowSums(zerohr_cpm >= cutoff) >= nsample
  passfilter <- (rowSums(cpm_0hr >= cutoff) >= nsample |
                   rowSums(cpm_1hr >= cutoff) >= nsample |
                   rowSums(cpm_6hr >= cutoff) >= nsample)
  return(passfilter)
}
cnfV <- function(y, design) {
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = T)
  return(v)  
}
contrastFit <- function(fit, design) {
  contr.matrix <- makeContrasts(
    casevscontrol = diseasecase-diseasecontrol,
    levels = colnames(design))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit <- eBayes(fit, trend = T)
  print(summary(decideTests(fit)))
  return(fit)
}

# Limma with voom####
covar_table <- covar_table[order(covar_table$cell_line), ]
designs <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  type <- str_split(names(adj_mat_lst)[i], pattern = "-", n = 2, simplify = T)[,1]
  if (type == "GABA") {
    f = covar_table$GABA_fraction
  } else if (type == "NEFM_neg_glut") {
    f = covar_table$nmglut_fraction
  } else {
    f = covar_table$npglut_fraction
  }
  designs[[i]] <- model.matrix(~0 + disease + group + age + sex + time + f, 
                               data = covar_table)
}

DGELists <- vector("list", length(adj_mat_lst))
voomedObjs <- vector("list", length(adj_mat_lst))
fits <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  DGELists[[i]] <- createDGE(adj_mat_lst[[i]]) # plug in adjusted or raw matrix
  ind.keep <- filterByCpm(DGELists[[i]], 1, 21)
  print(sum(ind.keep)) # 17496-17667-17006
  voomedObjs[[i]] <- cnfV(DGELists[[i]][ind.keep, ], designs[[i]])
  fits[[i]] <- lmFit(voomedObjs[[i]], designs[[i]])
}
names(fits) <- names(adj_mat_lst)
save(fits, file = "5+18+20_colby_linextime_fitted_values_objects.RData")

# extract fitted values ####
fits_0hr <- vector("list", length(adj_mat_lst))
fits_1hr <- vector("list", length(adj_mat_lst))
fits_6hr <- vector("list", length(adj_mat_lst))
for (i in 1:length(fits)) {y
  ftd <- fitted.MArrayLM(fits[[i]])
  fits_0hr[[i]] <- ftd[, seq(1, ncol(ftd), by = 3)]
  fits_1hr[[i]] <- ftd[, seq(2, ncol(ftd), by = 3)]
  fits_6hr[[i]] <- ftd[, seq(3, ncol(ftd), by = 3)]
}
names(fits_0hr) <- names(adj_mat_lst)
names(fits_1hr) <- names(adj_mat_lst)
names(fits_6hr) <- names(adj_mat_lst)
fitted_1v0 <- vector("list", length(adj_mat_lst))
fitted_6v0 <- vector("list", length(adj_mat_lst))
for (i in 1:length(fits_0hr)) {
  fitted_1v0[[i]] <- fits_1hr[[i]] - fits_0hr[[i]]
  fitted_6v0[[i]] <- fits_6hr[[i]] - fits_0hr[[i]]
}
names(fitted_1v0) <- names(fits_1hr)
names(fitted_6v0) <- names(fits_6hr)
save(fitted_1v0, fitted_6v0, file = "5+18+20line_fitted_1v0_6v0_matrices.RData")

# re-fit linear model ####
designs_1v0 <- vector("list", length(fitted_1v0))
designs_6v0 <- vector("list", length(fitted_6v0))

for (i in 1:length(designs_1v0)) {
  type <- names(fitted_1v0)[i]
  if (type == "GABA") {
    fraction_0hr <- covar_table$GABA_fraction[covar_table$time == "0hr"]
    fraction_1hr <- covar_table$GABA_fraction[covar_table$time == "1hr"]
  } else if (type == "NEFM_neg_glut") {
    fraction_0hr <- covar_table$nmglut_fraction[covar_table$time == "0hr"]
    fraction_1hr <- covar_table$nmglut_fraction[covar_table$time == "1hr"]
  } else {
    fraction_0hr <- covar_table$npglut_fraction[covar_table$time == "0hr"]
    fraction_1hr <- covar_table$npglut_fraction[covar_table$time == "1hr"]
}
  designs_1v0[[i]] <- model.matrix(~0 + disease + age + sex + group + 
                                     fraction_0hr + fraction_1hr, 
                                   data = covar_table[covar_table$time == "0hr", ])
}
for (i in 1:length(designs_6v0)) {
  type <- names(fitted_6v0)[i]
  if (type == "GABA") {
    fraction_0hr <- covar_table$GABA_fraction[covar_table$time == "0hr"]
    fraction_6hr <- covar_table$GABA_fraction[covar_table$time == "6hr"]
  } else if (type == "NEFM_neg_glut") {
    fraction_0hr <- covar_table$nmglut_fraction[covar_table$time == "0hr"]
    fraction_6hr <- covar_table$nmglut_fraction[covar_table$time == "6hr"]
  } else {
    fraction_0hr <- covar_table$npglut_fraction[covar_table$time == "0hr"]
    fraction_6hr <- covar_table$npglut_fraction[covar_table$time == "6hr"]
  }
  designs_6v0[[i]] <- model.matrix(~0 + disease + age + sex + group + 
                                     fraction_0hr + fraction_6hr, 
                                   data = covar_table[covar_table$time == "0hr", ])
}
names(designs_1v0) <- names(fitted_1v0)
names(designs_6v0) <- names(fitted_6v0)

contrastOnSubtractedVal <- function(df, design) {
  fit <- lmFit(df, design)
  contr.matrix <- makeContrasts(
    casevscontrol = diseasecase-diseasecontrol,
    levels = colnames(design))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit <- eBayes(fit, trend = T)
  print(summary(decideTests(fit)))
  return(fit)
}
contrasts_1v0 <- vector("list", 3)
contrasts_6v0 <- vector("list", 3)

for (i in 1:length(fitted_1v0)) {
  contrasts_1v0[[i]] <- contrastOnSubtractedVal(fitted_1v0[[i]], designs_1v0[[i]])
  contrasts_6v0[[i]] <- contrastOnSubtractedVal(fitted_6v0[[i]], designs_6v0[[i]])
}

# output results ####
setwd("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/case_v_control_DE")
cnames <- c("1v0 up", "1v0 down", "6v0 up", "6v0 down")
types <- names(fitted_1v0)
filenameElements <- c("res_1v0", "res_6v0")
for (i in 1:3){
  res_1v0 <- topTable(contrasts_1v0[[i]], sort.by = "P", number = Inf)
  res_6v0 <- topTable(contrasts_6v0[[i]], sort.by = "P", number = Inf)
  res_1v0$gene <- rownames(res_1v0)
  res_6v0$gene <- rownames(res_6v0)
  res_all <- list(res_1v0, res_6v0)
  # for (j in 1:length(res_all)) {
  #   filename <- paste0("./unfiltered_casevcontrol_DE/", types[i], "_", 
  #                      filenameElements[j], "_full_DEGs.csv")
  #   print(filename)
  #   write.table(res_all[[j]], file = filename, 
  #               quote = F, sep = ",", row.names = F, col.names = T)
  # }
  for (j in 1:length(res_all)) {
    filename <- paste0("./upregulated_significant/", types[i], "_", 
                       filenameElements[j], "_upregulated_genes_only.txt")
    print(filename)
    write.table(res_all[[j]]$gene[res_all[[j]]$logFC > 0 & res_all[[j]]$adj.P.Val < 0.05], 
                file = filename, quote = F, sep = "\t", row.names = F, col.names = F)
    
    # filename <- paste0("./upregulated_significant/", types[i], "_", 
    #                    filenameElements[j], "_upregulated_DEGs.csv")
    # print(filename)
    # write.table(res_all[[j]][res_all[[j]]$logFC > 0 & res_all[[j]]$adj.P.Val < 0.05, ], 
    #             file = filename, quote = F, sep = ",", row.names = F, col.names = T)
    
    filename <- paste0("./downregulated_significant/", types[i], "_", 
                       filenameElements[j], "_downregulated_genes_only.txt")
    print(filename)
    write.table(res_all[[j]]$gene[res_all[[j]]$logFC < 0 & res_all[[j]]$adj.P.Val < 0.05], 
                file = filename, quote = F, sep = "\t", row.names = F, col.names = F)
    
    # filename <- paste0("./downregulated_significant/", types[i], "_", 
    #                    filenameElements[j], "_downregulated_DEGs.csv")
    # print(filename)
    # write.table(res_all[[j]][res_all[[j]]$logFC < 0 & res_all[[j]]$adj.P.Val < 0.05, ], 
    #             file = filename, quote = F, sep = ",", row.names = F, col.names = T)
  }
} 

# summary table ####
deg_counts <- array(dim = c(3, 4), dimnames = list(types, c("up", "down", "nonsig", "total"))) 
for (i in 1:3){                                                               
  res <- topTable(contrasts_6v0[[i]], sort.by = "P", number = Inf, adjust.method = "BH")
  deg_counts[i, 1] <- sum(res$logFC > 0 & res$adj.P.Val < 0.05)
  deg_counts[i, 2] <- sum(res$logFC < 0 & res$adj.P.Val < 0.05)     
  deg_counts[i, 3] <- sum(res$adj.P.Val > 0.05)
  deg_counts[i, 4] <- nrow(res)
}
write.table(deg_counts, file = "./casevcontrol_6-0_deg_summary.csv", quote = F, sep = ",", 
            row.names = T, col.names = T)

# plot ####
hist(res$logFC)
max(res$logFC)
for (i in 1:length(contrasts_1v0)){
  res_1v0 <- topTable(contrasts_1v0[[i]], sort.by = "P", number = Inf)
  res_6v0 <- topTable(contrasts_6v0[[i]], sort.by = "P", number = Inf)
  res_1v0$gene <- rownames(res_1v0)
  res_6v0$gene <- rownames(res_6v0)
  res_all <- list(res_1v0, res_6v0) 
  for (k in 1:length(res_all)) {
    res <- res_all[[k]]
    res$significance <- "nonsignificant"
    res$significance[res$adj.P.Val < 0.05 & res$logFC > 0] <- "up"
    res$significance[res$adj.P.Val < 0.05 & res$logFC < 0] <- "down"
    unique(res$significance)
    res$significance <- factor(res$significance, levels = c("up", "nonsignificant", "down"))
    res$neg_log_pval <- (0 - log2(res$P.Value))
    res$labelling <- ""
    for (j in c("NAB2", "GRIN2A", "IMMP2L", "SNAP91", "PTPRK", "SP4", "RORB", "TCF4", "DIP2A")) {
      res$labelling[res$gene %in% j] <- j
    }
    pdf(paste(types[i], "casevcontrol", filenameElements[k], "volcano_plot.pdf", sep = "_"))
    p <- ggplot(data = as.data.frame(res),
              aes(x = logFC, 
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
    xlim(c(-0.05, 0.05)) +
    ggtitle(paste(types[i], "casevcontrol", filenameElements[k], sep = "-"))
    print(p)
    dev.off()
  }
}
