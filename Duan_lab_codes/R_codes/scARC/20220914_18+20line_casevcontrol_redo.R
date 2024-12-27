# Chuxuan Li 09/14/2022
# Use 18+20 line only, repeat case v control using residuals (methods from 
#Translational Psychiatry)

# init ####
library(stringr)
library(readr)
library(Seurat)
library(limma)
library(edgeR)
library(DESeq2)
library(sva)

# load data
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

celltypes <- sort(unique(RNA_18line$cell.type))
times <- c("0hr", "1hr", "6hr")
type_obj_lsts <- vector("list", 3)
for (i in 1:length(celltypes)) {
  t <- celltypes[i]
  print(t)
  sub_18 <- subset(RNA_18line, cell.type == t)
  sub_20 <- subset(RNA_20line, cell.type == t)
  subobjs <- list(sub_18, sub_20)
  zero_objs <- vector("list", 2)
  one_objs <- vector("list", 2)
  six_objs <- vector("list", 2)
  for (j in 1:2) {
    zero_objs[j] <- subset(subobjs[[j]], time.ident == "0hr") 
    one_objs[j] <- subset(subobjs[[j]], time.ident == "1hr")
    six_objs[j] <- subset(subobjs[[j]], time.ident == "6hr")
  }
  names(zero_objs) <- c("18", "20")
  names(one_objs) <- c("18", "20")
  names(six_objs) <- c("18", "20")
  type_obj_lsts[[i]] <- list(zero_objs, one_objs, six_objs)
  names(type_obj_lsts[[i]]) <- times
}
names(type_obj_lsts) <- celltypes

# make count matrices ####
#total 9, each time point has 3 cell types, each matix has 5+18+20=43 lines
#which is 43 samples, and because time point is no longer factored into the samples,
#43 is the total number of columns in each data frame.

makeMat4CombatseqLine <- function(typeobj, lines) {
  colnames <- rep_len(NA, (length(lines)))
  for (i in 1:length(lines)) {
    print(lines[i])
    colnames[i] <- lines[i]
    obj <- subset(typeobj, cell.line.ident == lines[i])
    if (i == 1) {
      mat <- as.array(rowSums(obj@assays$RNA@counts))
    } else {
      to_bind <- as.array(rowSums(obj@assays$RNA@counts))
      mat <- cbind(mat, to_bind)
    }
  }
  colnames(mat) <- colnames
  return(mat)
}

mat_lst <- vector("list", 6)

for (i in 1:length(type_obj_lsts)) {
  print(names(type_obj_lsts)[i])
  for (j in 1:length(type_obj_lsts)) {
    print(names(type_obj_lsts[[i]])[j])
    objlst <- type_obj_lsts[[i]][[j]]
    for (k in 1:length(objlst)) {
      obj <- objlst[[k]]
      ls <- sort(unique(obj$cell.line.ident))
      if (k == 1) {
        mat <- makeMat4CombatseqLine(objlst[[k]], ls)
      } else {
        mat <- cbind(mat, makeMat4CombatseqLine(objlst[[k]], ls))
      }
    }
    print(i * 3 - 3 + j)
    mat <- mat[, order(colnames(mat))]
    mat_lst[[i * 3 - 3 + j]] <- mat
    names(mat_lst)[i * 3 - 3 + j] <- paste(names(type_obj_lsts)[i], 
                                           names(type_obj_lsts[[i]])[j], sep = "-")
  }
}
lines <- colnames(mat_lst[[1]])

load("../../../covariates_pooled_together_for_all_linextime_5+18+20.RData")
for (i in 1:length(lines)) {
  print(lines[i])
  if (i == 1) {
    batch = unique(covar_table$group[covar_table$cell_line == lines[i]])
  } else {
    batch = c(batch, unique(covar_table$group[covar_table$cell_line == lines[i]]))
  }
}
batch <- as.factor(batch)

adj_mat_lst <- vector("list", length(mat_lst)) 
for (i in 1:length(mat_lst)) {
  adj_mat_lst[[i]] <- ComBat_seq(mat_lst[[i]], batch = batch)
}
names(adj_mat_lst) <- names(mat_lst)
save(adj_mat_lst, file = "5+18+20line_sepby_type_time_colby_line_adj_mat_list.RData")

# make design matrices ####
covar_table <- covar_table[order(covar_table$cell_line), ]
covar_table_0hr <- covar_table[covar_table$time == "0hr", ]
covar_table_1hr <- covar_table[covar_table$time == "1hr", ]
covar_table_6hr <- covar_table[covar_table$time == "6hr", ]

designs <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  type <- str_split(names(adj_mat_lst)[i], pattern = "-", n = 2, simplify = T)[,1]
  time <- str_split(names(adj_mat_lst)[i], pattern = "-", n = 2, simplify = T)[,2]
  if (time == "0hr") {
    d = covar_table_0hr
  } else if (time == "1hr") {
    d = covar_table_1hr
  } else {
    d = covar_table_6hr
  }
  if (type == "GABA") {
    f = d$GABA_fraction
  } else if (type == "NEFM_neg_glut") {
    f = d$nmglut_fraction
  } else {
    f = d$npglut_fraction
  }
  designs[[i]] <- model.matrix(~0 + disease + age + sex + group + f, 
                               data = d)
}

# do limma voom ####
createDGE <- function(count_matrix){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = lines)
  A <- rowSums(y$counts)
  hasant <- rowSums(is.na(y$genes)) == 0
  y <- y[hasant, , keep.lib.size = F]
  print(dim(y)) 
  return(y)
}
filterByCpm <- function(y, cutoff, nsample) {
  cpm <- cpm(y)
  passfilter <- (rowSums(cpm >= cutoff) >= nsample)
  return(passfilter)
}
cnfV <- function(y, design) {
  y <- calcNormFactors(y)
  # voom (no logcpm needed, voom does it within the function)
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
designs[[1]]

DGEs <- vector("list", length(adj_mat_lst))
voomedDGEs <- vector("list", length(adj_mat_lst))
fits <- vector("list", length(adj_mat_lst))
for (i in 1:length(adj_mat_lst)) {
  DGEs[[i]] <- createDGE(adj_mat_lst[[i]])
  voomedDGEs[[i]] <- cnfV(DGEs[[i]][filterByCpm(DGEs[[i]], 1, 43), ], designs[[i]])
  fit <- lmFit(voomedDGEs[[i]], designs[[i]])
  fits[[i]] <- contrastFit(fit, designs[[i]])
}