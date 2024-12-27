{
  # init ####
  library(limma)
  library(edgeR)
  library(sva)
  library(ggplot2)
  library(RColorBrewer)
  library(stringr)
  library(readr)
  library(Signac)
  
  # library(readr)#library(readr)
  # library(ggplot2)
  library(ggrepel)
  # library(RColorBrewer)
  library(dplyr)
  library(graphics)
  # library(stringr)
}

set.seed(2022)
# load annotated files ####

# make fit_all var #####

load("frag_count_by_cell_type_specific_peaks_22Jul2022.RData")

load("covariates_pooled_together_for_all_linextime_5+18+20.RData")

# define covariates from covar table
covar_18line <- covar_table[covar_table$cell_line %in% lines, ]
covar_18line <- covar_18line[order(covar_18line$cell_line), ]
rm(covar_table)
GABA_fraction <- covar_18line$GABA_fraction
nmglut_fraction <- covar_18line$nmglut_fraction
npglut_fraction <- covar_18line$npglut_fraction

# make combat dataframes ####
GABA <- Seurat_object_list[[2]]
cellines <- sort(unique(GABA$cell.line.ident))
times <- sort(unique(GABA$time.ident))
GABA$timexline.ident <- NA
for (l in cellines) {
  for (t in times) {
    id <- paste(l, t, sep = "-")
    print(id)
    GABA$timexline.ident[GABA$cell.line.ident == l &
                           GABA$time.ident == t] <- id
  }
}


nmglut <- Seurat_object_list[[3]]
cellines <- sort(unique(nmglut$cell.line.ident))
times <- sort(unique(nmglut$time.ident))
nmglut$timexline.ident <- NA
for (l in cellines) {
  for (t in times) {
    id <- paste(l, t, sep = "-")
    print(id)
    nmglut$timexline.ident[nmglut$cell.line.ident == l &
                             nmglut$time.ident == t] <- id
  }
}

npglut <- Seurat_object_list[[4]]
cellines <- sort(unique(npglut$cell.line.ident))
times <- sort(unique(npglut$time.ident))
npglut$timexline.ident <- NA
for (l in cellines) {
  for (t in times) {
    id <- paste(l, t, sep = "-")
    print(id)
    npglut$timexline.ident[npglut$cell.line.ident == l &
                             npglut$time.ident == t] <- id
  }
}




unique(nmglut$timexline.ident)
lts <- sort(unique(nmglut$timexline.ident))

# make raw count matrix ####
makeMat4CombatseqLine <- function(typeobj, linetimes) {
  rownames <- rep_len(NA, (length(linetimes)))
  for (i in 1:length(linetimes)) {
    print(linetimes[i])
    rownames[i] <- linetimes[i]
    obj <- subset(typeobj, timexline.ident == linetimes[i])
    count_matrix <- obj 
    if (i == 1) {
      mat <- as.array(rowSums(obj@assays$ATAC@counts))
    } else {
      to_bind <- as.array(rowSums(obj@assays$ATAC@counts))
      mat <- rbind(mat, to_bind)
    }
  }
  rownames(mat) <- rownames
  return(mat)
}

GABA_mat <- t(makeMat4CombatseqLine(GABA, lts))
nmglut_mat <- t(makeMat4CombatseqLine(nmglut, lts))
npglut_mat <- t(makeMat4CombatseqLine(npglut, lts))

# combat-seq ####
batch = rep(c("12", "11", "12", "12", "8", "8", "9", "9", "10", "8", "8", "9", 
              "10", "10", "11", "11", "11", "10"), each = 3)
batch
length(batch) #54
ncol(GABA_mat) #54

GABA_mat_adj <- ComBat_seq(GABA_mat, batch = batch, 
                           group = NULL)
nmglut_mat_adj <- ComBat_seq(nmglut_mat, batch = batch,
                             group = NULL)
npglut_mat_adj <- ComBat_seq(npglut_mat, batch = batch,
                             group = NULL)
save(GABA_mat_adj, nmglut_mat_adj, npglut_mat_adj, file = "post-combat_count_matrices.RData")


load("post-combat_count_matrices.RData")


celllines <- str_extract(lts, "^CD_[0-9][0-9]")


# auxillary functions ####
createDGE <- function(count_matrix){
  y <- DGEList(counts = count_matrix, genes = rownames(count_matrix), group = celllines)
  A <- rowSums(y$counts)
  #isexpr <- A > 10
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
  # voom (no logcpm)
  v <- voom(y, design, plot = T)
  return(v)  
}
contrastFit <- function(fit, design) {
  contr.matrix <- makeContrasts(
    onevszero = time1hr-time0hr, 
    sixvszero = time6hr-time0hr, 
    sixvsone = time6hr-time1hr, levels = colnames(design))
  fit <- contrasts.fit(fit, contrasts = contr.matrix)
  fit <- eBayes(fit, trend = T)
  print(summary(decideTests(fit)))
  return(fit)
}


# Limma with voom####
design_GABA <- model.matrix(~0 + time + group + age + sex + disease + GABA_fraction, data = covar_18line) #26353 genes
y_GABA <- createDGE(GABA_mat_adj)
ind.keep <- filterByCpm(y_GABA, 1, 9)
sum(ind.keep) #311555 -> 172759
v_GABA <- cnfV(y_GABA[ind.keep, ], design_GABA)
fit_GABA <- lmFit(v_GABA, design_GABA)
fit_GABA_contr <- contrastFit(fit_GABA, design_GABA)
plotSA(fit_GABA_contr, main = "Final model: Mean-variance trend") # same, just has a line

design_nmglut <- model.matrix(~0 + time + group + age + sex + disease + nmglut_fraction, data = covar_18line) #26353 genes
y_nmglut <- createDGE(nmglut_mat_adj)
ind.keep <- filterByCpm(y_nmglut, 1, 9)
sum(ind.keep) #325571 -> 196083
v_nmglut <- cnfV(y_nmglut[ind.keep, ], design_nmglut)
fit_nmglut <- lmFit(v_nmglut, design_nmglut)
fit_nmglut_contr <- contrastFit(fit_nmglut, design_nmglut)
# plotSA(fit_nmglut_contr, main = "Final model: Mean-variance trend") 

design_npglut <- model.matrix(~0 + time + group + age + sex + disease + npglut_fraction, data = covar_18line) #29400 genes
y_npglut <- createDGE(npglut_mat_adj)
ind.keep <- filterByCpm(y_npglut, 1, 9)
sum(ind.keep) #340779 -> 207450
v_npglut <- cnfV(y_npglut[ind.keep, ], design_npglut)
fit_npglut <- lmFit(v_npglut, design_npglut)
fit_npglut_contr <- contrastFit(fit_npglut, design_npglut)
# plotSA(fit_npglut_contr, main = "Final model: Mean-variance trend") 


# output results ####
# setwd("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/da_peaks_limma")

fit_all <- list(fit_GABA_contr, fit_nmglut_contr, fit_npglut_contr)
saveRDS(fit_all,
        file = "fit_all.RDs")
# cnames <- c("1v0 up", "1v0 down", "6v0 up", "6v0 down")
# types <- c("GABA", "nmglut", "npglut")

#####

pathlist <-
  list.files(path = "./da_peaks_limma/bed/annotated", 
             pattern = ".tsv", 
             full.names = T)
reslist <- vector("list", 
                  length(pathlist))

for (i in 1:length(pathlist)) {
  reslist[[i]] <- 
    read_delim(pathlist[i], 
               delim = "\t", 
               col_names = F, 
               skip = 1, 
               escape_double = F)
}

for (i in 1:length(reslist)) {
  reslist[[i]] <- 
    reslist[[i]][reslist[[i]]$X18 %in% c("FOS", 
                                         "NPAS4", 
                                         "BDNF"), ]
}
# reslist[[2]]$X18

# extract peaks ####
onegenelist <- vector("list", length(reslist))

for (i in 1:length(reslist)) {
  genes.exist <- 
    unique(reslist[[i]]$X18)
  if (length(genes.exist) != 0) {
    for (j in 1:length(genes.exist)) {
      tempres <- 
        reslist[[i]][reslist[[i]]$X18 == genes.exist[j], ]
      tempres <- 
        tempres[order(tempres$X1, decreasing = T), ]
      tempres <- 
        tempres[1, ]
      tempres[, 19] <- 
        paste(tempres$X4, tempres$X5 - 1, tempres$X6, sep = "-")
      tempres <- 
        tempres[, 18:19]
      if (j == 1) {
        onegenetemp <- tempres
      } else {
        onegenetemp <- 
          rbind(onegenetemp, 
                tempres)
      }
    }
    onegenelist[[i]] <- 
      as.data.frame(onegenetemp)
  }
}

names(onegenelist) <- 
  str_extract(pathlist, 
              "[A-Za-z]+_[1|6]v0hr_[a-z]+")
onegenelist <- 
  onegenelist[seq(2, length(onegenelist), 2)]
genes.use.1v0 <- 
  onegenelist[seq(1, 6, 2)]
genes.use.6v0 <- 
  onegenelist[seq(2, 6, 2)]

# plot ####
names_1v0 <- c("GABA_1v0hr", 
               "nmglut_1v0hr", 
               "npglut_1v0hr")
names_6v0 <- c("GABA_6v0hr", 
               "nmglut_6v0hr", 
               "npglut_6v0hr")
colors <- c("steelblue3", 
            "grey", 
            "indianred3")
da_peaks_by_time_list_1v0 <- 
  vector("list", length(fit_all))
da_peaks_by_time_list_6v0 <- 
  vector("list", length(fit_all))
for (i in 1:length(fit_all)) {
  da_peaks_by_time_list_1v0[[i]] <- topTable(fit_all[[i]], coef = "onevszero", p.value = Inf, 
                                             sort.by = "P", number = Inf)
  da_peaks_by_time_list_6v0[[i]] <- topTable(fit_all[[i]], coef = "sixvszero", p.value = Inf, 
                                             sort.by = "P", number = Inf)
}

# volcano plots ####
for (i in 1:3) {
  name1 <- names_1v0[i]
  print(name1)
  da_peaks_by_time_list_1v0[[i]]$peak <- rownames(da_peaks_by_time_list_1v0[[i]])
  da_peaks_by_time_list_1v0[[i]]$genes <- NULL
  df <- da_peaks_by_time_list_1v0[[i]][, c("peak", "AveExpr", "logFC", "P.Value", "adj.P.Val")]
  df$significance <- "nonsig"
  df$significance[df$adj.P.Val < 0.05 & df$logFC > 0] <- "pos"
  df$significance[df$adj.P.Val < 0.05 & df$logFC < 0] <- "neg"
  df$neg_log_p_val <- (0 - log10(df$P.Value))
  df$labelling <- ""
  for (j in 1:length(genes.use.1v0[[i]])){
    df$labelling[df$peak %in% genes.use.1v0[[i]]$X19] <- genes.use.1v0[[i]]$X18
  }
  print(unique(df$labelling))
  
  # plot
  file_name <- paste0("./da_peaks_limma/", name1, "_new.png")
  png(file_name)
  p <- ggplot(data = df,
              aes(x = logFC, 
                  y = neg_log_p_val, 
                  label = labelling,
                  color = significance)) + 
    geom_point(size = 0.1) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    ggtitle(str_replace_all(name1, "_", " ")) +
    geom_text_repel(box.padding = unit(0.5, 'lines'),
                    min.segment.length = 0,
                    force = 2,
                    force_pull = 10,
                    max.overlaps = 100000,
                    
                    color = "darkred",
                    size = 6,
                    show.legend = F) +
    xlim(-5, 5) +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          title = element_text(size = 20))
  print(p)
  dev.off()
  
  name6 <- names_6v0[i]
  print(name6)
  da_peaks_by_time_list_6v0[[i]]$peak <- rownames(da_peaks_by_time_list_6v0[[i]])
  da_peaks_by_time_list_6v0[[i]]$genes <- NULL
  df <- da_peaks_by_time_list_6v0[[i]][, c("peak", "AveExpr", "logFC", "P.Value", "adj.P.Val")]
  df$significance <- "nonsig"
  df$significance[df$adj.P.Val < 0.05 & df$logFC > 0] <- "pos"
  df$significance[df$adj.P.Val < 0.05 & df$logFC < 0] <- "neg"
  df$neg_log_p_val <- (0 - log10(df$P.Value))
  df$labelling <- ""
  for (j in 1:length(genes.use.6v0[[i]])){
    df$labelling[df$peak %in% genes.use.6v0[[i]]$X19] <- genes.use.6v0[[i]]$X18
  }
  print(unique(df$labelling))
  file_name <- paste0("./da_peaks_limma/", name6, "_new.png")
  png(file_name)
  p <- ggplot(data = df,
              aes(x = logFC, 
                  y = neg_log_p_val, 
                  label = labelling,
                  color = significance)) + 
    geom_point(size = 0.1) +
    scale_color_manual(values = colors) +
    theme_minimal() +
    ggtitle(str_replace_all(name6, "_", " ")) +
    geom_text_repel(box.padding = unit(0.5, 'lines'),
                    min.segment.length = 0,
                    force = 2,
                    force_pull = 10,
                    max.overlaps = 100000,
                    size = 6,
                    color = "darkred",
                    show.legend = F) +
    xlim(-5, 5) +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          title = element_text(size = 20))
  print(p)
  dev.off()
}
