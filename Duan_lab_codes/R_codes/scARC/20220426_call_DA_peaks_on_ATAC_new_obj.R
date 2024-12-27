# Chuxuan Li 04/26/2022
# Using the 18-line ATAC_new object that has RNAseq cell type labels, call 
#differentially accessible peaks between 0 and 1hr, 0 and 6hr for each of the 3
#cell types

# init ####
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)

library(future)
library(readr)
library(stringr)

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(graphics)

set.seed(1886)
# set threads and parallelization
plan("multisession", workers = 2)
options(expressions = 20000)
options(future.globals.maxSize = 11474836480)
load("~/NVME/scARC_Duan_018/Duan_project_024_ATAC/ATAC_new_obj_with_MACS2_new_peaks_and_RNA_labels.RData")

# list idents ####
unique(ATAC_new$timextype.ident)
idents_list <- c("GABA_0hr", "GABA_1hr", "GABA_6hr",
                 "NEFM_neg_glut_0hr", "NEFM_neg_glut_1hr", "NEFM_neg_glut_6hr",
                 "NEFM_pos_glut_0hr", "NEFM_pos_glut_1hr", "NEFM_pos_glut_6hr")

ATAC <- subset(ATAC_new, timextype.ident %in% idents_list)
DefaultAssay(ATAC) <- "peaks"

# ATAC@assays$peaks@ranges[ATAC@assays$peaks@ranges@seqnames == "chr11" &
#                            ATAC@assays$peaks@ranges@ranges@start > 27770000, ]
# sum(ATAC@assays$peaks@counts@Dimnames[[1]] == "chr11-27770213-27771025")


# auxillary functions ####
printAPeak <- function(peaklist1, peaklist6, sequence, 
                       rdataname1 = NA, rdataname6 = NA){
  #peaklist1: a list containing da peaks called for 1v0 hr
  #peaklist6: a list containing da peaks called for 6v0 hr, must have same length
  #as peaklist1
  #sequence: a string, part of the sequence you are looking for
  #rdataname1/6: a string, the name of the RData file if you want to save the RData
  for (i in 1:length(peaklist1)){
    print(paste0(i, " 1v0"))
    peaklist1[[i]]$peak <- rownames(peaklist1[[i]])
    print(peaklist1[[i]][str_detect(peaklist1[[i]]$peak, sequence), ])
    
    print(paste0(i, " 6v0"))
    peaklist6[[i]]$peak <- rownames(peaklist6[[i]])
    print(peaklist6[[i]][str_detect(peaklist6[[i]]$peak, sequence), ])
  }
  if (!is.na(rdataname1) & !is.na(rdataname6)){
    save(peaklist1, file = rdataname1)
    save(peaklist6, file = rdataname6)
  }
}
saveDApeaksResults <- function(peaklist1, peaklist6, folder){
  #peaklist1: a list containing da peaks called for 1v0 hr
  #peaklist6: a list containing da peaks called for 6v0 hr, must have same length
  #as peaklist1
  #folder: a string, the name of the output folder
  idents_list <- c("GABA_0hr", "GABA_1hr", "GABA_6hr",
                   "NEFM_neg_glut_0hr", "NEFM_neg_glut_1hr", "NEFM_neg_glut_6hr",
                   "NEFM_pos_glut_0hr", "NEFM_pos_glut_1hr", "NEFM_pos_glut_6hr")
  for (i in 1:length(peaklist1)){
    # get the correct indices
    ind0 <- 3 * i - 2
    ident0 <- idents_list[ind0]
    ind1 <- 3 * i - 1
    ident1 <- idents_list[ind1]
    ind6 <- 3 * i
    ident6 <- idents_list[ind6]
    print(paste0("0hr: ", ident0, " 1hr: ", ident1, " 6hr: ", ident6))
    
    df <- data.frame(peaks = rownames(peaklist1[[i]]),
                     p_val = peaklist1[[i]]$p_val,
                     avg_log2FC = peaklist1[[i]]$avg_log2FC,
                     pct.1 = peaklist1[[i]]$pct.1,
                     pct.2 = peaklist1[[i]]$pct.2,
                     p_val_adj = peaklist1[[i]]$p_val_adj)
    file_name <- paste0("./", folder, "/csv/",
                        ident1, "_",
                        ident0, ".csv")
    write.table(df,
                file = file_name,
                quote = F,
                sep = ",",
                row.names = F,
                col.names = T)
    file_name <- str_replace_all(file_name, "csv", "bed")
    df_bed <- data.frame(chr = str_split(string = df$peaks, pattern = "-", simplify = T)[, 1],
                         start = gsub(pattern = "-",
                                      replacement = "",
                                      x = str_extract_all(string = df$peaks,
                                                          pattern = "-[0-9]+-"
                                      )),
                         end = gsub(pattern = "-",
                                    replacement = "",
                                    x = str_extract_all(string = df$peaks,
                                                        pattern = "-[0-9]+$"
                                    )),
                         id = paste0(df$avg_log2FC, '^',
                                     df$p_val, '^',
                                     df$p_val_adj), # log2FC^pval^pval_adj, for later use
                         p_val = df$p_val,
                         strand = rep("+", length(df$p_val)))
    write.table(df_bed,
                file = file_name,
                quote = F,
                sep = "\t",
                row.names = F,
                col.names = F)
    
    # 6hr
    df <- data.frame(peaks = rownames(peaklist6[[i]]),
                     p_val = peaklist6[[i]]$p_val,
                     avg_log2FC = peaklist6[[i]]$avg_log2FC,
                     pct.1 = peaklist6[[i]]$pct.1,
                     pct.2 = peaklist6[[i]]$pct.2,
                     p_val_adj = peaklist6[[i]]$p_val_adj)
    file_name <- paste0("./", folder, "/csv/",
                        ident6, "_",
                        ident0, ".csv")
    write.table(df,
                file = file_name,
                quote = F,
                sep = ",",
                row.names = F,
                col.names = F)
    file_name <- str_replace_all(file_name, "csv", "bed")
    df_bed <- data.frame(chr = str_split(string = df$peaks, pattern = "-", simplify = T)[, 1],
                         start = gsub(pattern = "-",
                                      replacement = "",
                                      x = str_extract_all(string = df$peaks,
                                                          pattern = "-[0-9]+-"
                                      )),
                         end = gsub(pattern = "-",
                                    replacement = "",
                                    x = str_extract_all(string = df$peaks,
                                                        pattern = "-[0-9]+$"
                                    )),
                         id = paste0(df$avg_log2FC, '^',
                                     df$p_val, '^',
                                     df$p_val_adj), # log2FC^pval^pval_adj, for later use
                         p_val = df$p_val,
                         strand = rep("+", length(df$p_val)))
    write.table(df_bed,
                file = file_name,
                quote = F,
                sep = "\t",
                row.names = F,
                col.names = F)
  }
} 

# summary table
deg_counts <- array(dim = c(3, 4), dimnames = list(types, cnames)) # 3 cell types are rows, 1/6 up/down are cols
for (i in 1:3){                                                               
  res_1v0 <- topTable(fit_all[[i]], coef = "onevszero", p.value = 0.05, sort.by = "P", number = Inf)
  res_6v0 <- topTable(fit_all[[i]], coef = "sixvszero", p.value = 0.05, sort.by = "P", number = Inf)
  deg_counts[i, 1] <- sum(res_1v0$logFC > 0)
  deg_counts[i, 2] <- sum(res_1v0$logFC < 0)     
  deg_counts[i, 3] <- sum(res_6v0$logFC > 0)
  deg_counts[i, 4] <- sum(res_6v0$logFC < 0)
}
write.table(deg_counts, file = "./deg_summary_combat.csv", quote = F, sep = ",", 
            row.names = T, col.names = T)
# test with chr11 ####
# this chunk of codes is for testing only and the output was not used in the 
#actual analysis

counts <- GetAssayData(ATAC, assay = "peaks")
counts <- counts[which(str_detect(rownames(counts), "chr11")), ]
ATAC_11 <- subset(ATAC, features = rownames(counts))

ATAC_11$broad.cell.type <- "others"
ATAC_11$broad.cell.type[ATAC_11$cell.type == "GABA"] <- "GABA"
ATAC_11$broad.cell.type[ATAC_11$cell.type %in% c("NEFM_neg_glut",
                                                 "NEFM_pos_glut")] <- "glut"
ATAC_11$broad.typextime <- NA
types <- unique(ATAC_11$broad.cell.type)
times <- unique(ATAC_11$time.ident)
for (i in types){
  for (j in times){
    ATAC_11$broad.typextime[ATAC_11$broad.cell.type == i &
                              ATAC_11$time.ident == j] <- paste(i, j, sep = "_")
  }
}
unique(ATAC_11$broad.typextime)
idents_list <- c("GABA_0hr", "GABA_1hr", "GABA_6hr",
                 "glut_0hr", "glut_1hr", "glut_6hr")

Idents(ATAC_11) <- "broad.typextime"
DefaultAssay(ATAC_11)
chr11_peaks_1v0 <- vector(mode = "list", length = 2L)
chr11_peaks_6v0 <- vector(mode = "list", length = 2L)
for (i in 1:2){
  # get the correct indices
  ind0 <- 3 * i - 2
  ident0 <- idents_list[ind0]
  ind1 <- 3 * i - 1
  ident1 <- idents_list[ind1]
  ind6 <- 3 * i
  ident6 <- idents_list[ind6]
  print(paste0("0hr: ", ident0, " 1hr: ", ident1, " 6hr: ", ident6))
  
  # call FindMarkers to find differentially accessible peaks
  da_peaks_1v0 <- FindMarkers(
    object = ATAC_11,
    ident.1 = ident1,
    ident.2 = ident0,
    test.use = 'MAST',
    logfc.threshold = 0, 
    min.pct = 0.02,
    random.seed = 100)
  
  chr11_peaks_1v0[[i]] <- da_peaks_1v0
  
  da_peaks_6v0 <- FindMarkers(
    object = ATAC_11,
    ident.1 = ident6,
    ident.2 = ident0,
    test.use = 'MAST',
    logfc.threshold = 0,
    min.pct = 0.02,
    random.seed = 100)
  
  chr11_peaks_6v0[[i]] <- da_peaks_6v0
}

printAPeak(chr11_peaks_1v0, chr11_peaks_6v0, "chr11-277", 
           "chr11_peaks_1v0_0.01.RData", "chr11_peaks_6v0_0.01.RData")
saveDApeaksResults(chr11_peaks_1v0, chr11_peaks_6v0, "da_peaks_minpct0.01_chr11")
save(chr11_peaks_1v0, chr11_peaks_6v0, file = "combined_2_gluts_chr11_DApeaks.RData")


# call da peaks on all peaks ####
da_peaks_by_time_list_1v0 <- vector(mode = "list",
                                length = 3L)
da_peaks_by_time_list_6v0 <- vector(mode = "list",
                                    length = 3L)

Idents(ATAC) <- "timextype.ident"
unique(ATAC$cell.type)

for (i in 1:3){
  # get the correct indices
  ind0 <- 3 * i - 2
  ident0 <- idents_list[ind0]
  ind1 <- 3 * i - 1
  ident1 <- idents_list[ind1]
  ind6 <- 3 * i
  ident6 <- idents_list[ind6]
  print(paste0("0hr: ", ident0, " 1hr: ", ident1, " 6hr: ", ident6))
  
  # call FindMarkers to find differentially accessible peaks
  da_peaks_1v0 <- FindMarkers(
    object = ATAC,
    ident.1 = ident1,
    ident.2 = ident0,
    test.use = 'MAST',
    #features = "chr11-27770213-27771025", 
    logfc.threshold = 0, 
    min.pct = 0.01,
    random.seed = 100)

  da_peaks_by_time_list_1v0[[i]] <- da_peaks_1v0

  da_peaks_6v0 <- FindMarkers(
    object = ATAC,
    ident.1 = ident6,
    ident.2 = ident0,
    test.use = 'MAST',
    #features = "chr11-27770213-27771025",
    logfc.threshold = 0,
    min.pct = 0.01,
    random.seed = 100)

  da_peaks_by_time_list_6v0[[i]] <- da_peaks_6v0
}

printAPeak(da_peaks_by_time_list_1v0, da_peaks_by_time_list_6v0, "chr11-277")

save(da_peaks_by_time_list_1v0, da_peaks_by_time_list_6v0, 
     file = "da_peaks_by_time_1-6_dfs_0.01.RData")

# save 0.02
idents_list <- c("GABA_0hr", "GABA_1hr", "GABA_6hr",
                 "NEFM_neg_glut_0hr", "NEFM_neg_glut_1hr", "NEFM_neg_glut_6hr",
                 "NEFM_pos_glut_0hr", "NEFM_pos_glut_1hr", "NEFM_pos_glut_6hr")
for (i in 1:length(da_peaks_by_time_list_6v0)){
  # get the correct indices
  ind0 <- 3 * i - 2
  ident0 <- idents_list[ind0]
  ind1 <- 3 * i - 1
  ident1 <- idents_list[ind1]
  ind6 <- 3 * i
  ident6 <- idents_list[ind6]
  print(paste0("0hr: ", ident0, " 1hr: ", ident1, " 6hr: ", ident6))
  
  up <- da_peaks_by_time_list_6v0[[i]][da_peaks_by_time_list_6v0[[i]]$p_val_adj < 0.05 &
                                         da_peaks_by_time_list_6v0[[i]]$avg_log2FC > 0, ]
  do <- da_peaks_by_time_list_6v0[[i]][da_peaks_by_time_list_6v0[[i]]$p_val_adj < 0.05 &
                                         da_peaks_by_time_list_6v0[[i]]$avg_log2FC < 0, ]
  df <- data.frame(peaks = rownames(up),
                   p_val = up$p_val,
                   avg_log2FC = up$avg_log2FC,
                   pct.1 = up$pct.1,
                   pct.2 = up$pct.2,
                   p_val_adj = up$p_val_adj)
  file_name <- paste0("./",
                      ident6, "_",
                      ident0, "_upregulated_significant_peaks.csv")
  write.table(df,
              file = file_name,
              quote = F,
              sep = ",",
              row.names = F,
              col.names = T)
  
  df <- data.frame(peaks = rownames(do),
                   p_val = do$p_val,
                   avg_log2FC = do$avg_log2FC,
                   pct.1 = do$pct.1,
                   pct.2 = do$pct.2,
                   p_val_adj = do$p_val_adj)
  file_name <- paste0("./",
                      ident6, "_",
                      ident0, "_downregulated_significant_peaks.csv")
  write.table(df,
              file = file_name,
              quote = F,
              sep = ",",
              row.names = F,
              col.names = T)
}
dapeak_counts <- matrix(nrow = 6, ncol = 4, 
                                 dimnames = list(c("GABA 1v0", "GABA 6v0", 
                                                   "nmglut 1v0", "nmglut 6v0",
                                                   "npglut 1v0", "npglut 6v0"),
                                                 c("total", "significant", "up", "down")))
for (i in 1:length(da_peaks_by_time_list_6v0)){
  print(i)
  dapeak_counts[2*i-1, 1] <- nrow(da_peaks_by_time_list_1v0[[i]])
  dapeak_counts[2*i, 1] <- nrow(da_peaks_by_time_list_6v0[[i]])
  dapeak_counts[2*i-1, 2] <- sum(da_peaks_by_time_list_1v0[[i]]$p_val_adj < 0.05)
  dapeak_counts[2*i, 2] <- sum(da_peaks_by_time_list_6v0[[i]]$p_val_adj < 0.05)
  dapeak_counts[2*i-1, 3] <- sum(da_peaks_by_time_list_1v0[[i]]$p_val_adj < 0.05 & 
                                   da_peaks_by_time_list_1v0[[i]]$avg_log2FC > 0)
  dapeak_counts[2*i, 3] <- sum(da_peaks_by_time_list_6v0[[i]]$p_val_adj < 0.05 & 
                                 da_peaks_by_time_list_6v0[[i]]$avg_log2FC > 0)
  dapeak_counts[2*i-1, 4] <- sum(da_peaks_by_time_list_1v0[[i]]$p_val_adj < 0.05 & 
                                   da_peaks_by_time_list_1v0[[i]]$avg_log2FC < 0)
  dapeak_counts[2*i, 4] <- sum(da_peaks_by_time_list_6v0[[i]]$p_val_adj < 0.05 & 
                                 da_peaks_by_time_list_6v0[[i]]$avg_log2FC < 0)
}
write.table(dapeak_counts, file = "dapeak_counts_summary_minpct0.01.csv", 
            quote = F, sep = ",", row.names = T, col.names = T)

# save filtered peaks into bed ####
saveBed <- function(df, file_name){
  df_bed <- data.frame(chr = str_split(string = df$peaks, pattern = "-", simplify = T)[, 1],
                       start = gsub(pattern = "-",
                                    replacement = "",
                                    x = str_extract_all(string = df$peaks,
                                                        pattern = "-[0-9]+-"
                                    )),
                       end = gsub(pattern = "-",
                                  replacement = "",
                                  x = str_extract_all(string = df$peaks,
                                                      pattern = "-[0-9]+$"
                                  )),
                       id = paste0(df$avg_log2FC, '^',
                                   df$p_val, '^',
                                   df$p_val_adj), # log2FC^pval^pval_adj, for annotation
                       p_val = df$p_val,
                       strand = rep("+", length(df$p_val)))
  write.table(df_bed,
              file = file_name,
              quote = F,
              sep = "\t",
              row.names = F,
              col.names = F)
}


# divide object and calculate again, min.pct =0.01 ####
unique(Idents(ATAC))
chrs <- unique(str_extract(rownames(ATAC), "chr[0-9]+-"))
chrs <- chrs[!is.na(chrs)]
rownames(ATAC)[(str_detect(rownames(ATAC), chrs[11]))] # test

for (i in 1:length(chrs)) {
  peaks.use <- rownames(ATAC)[(str_detect(rownames(ATAC), chrs[i]))]
  obj <- ATAC[rownames(ATAC) %in% peaks.use, ]
  print(i)
  if (i == 1) {
    NEFM_neg_glut_1v0 <- FindMarkers(
      object = obj,
      ident.1 = "NEFM_neg_glut_1hr",
      ident.2 = "NEFM_neg_glut_0hr",
      test.use = 'MAST',
      logfc.threshold = 0, 
      min.pct = 0.01, 
      random.seed = 100)
    
    NEFM_neg_glut_6v0 <- FindMarkers(
      object = obj,
      ident.1 = "NEFM_neg_glut_6hr",
      ident.2 = "NEFM_neg_glut_0hr",
      test.use = 'MAST',
      logfc.threshold = 0,
      min.pct = 0.01,
      random.seed = 100)
  } else {
    NEFM_neg_glut_1v0 <- rbind(NEFM_neg_glut_1v0, FindMarkers(
      object = obj,
      ident.1 = "NEFM_neg_glut_1hr",
      ident.2 = "NEFM_neg_glut_0hr",
      test.use = 'MAST',
      logfc.threshold = 0, 
      min.pct = 0.01, 
      random.seed = 100))
    
    NEFM_neg_glut_6v0 <- rbind(NEFM_neg_glut_6v0, FindMarkers(
      object = obj,
      ident.1 = "NEFM_neg_glut_6hr",
      ident.2 = "NEFM_neg_glut_0hr",
      test.use = 'MAST',
      logfc.threshold = 0,
      min.pct = 0.01,
      random.seed = 100))
  }

}
save(NEFM_neg_glut_1v0, NEFM_neg_glut_6v0, file = "nmglut_da_peaks_minpct0.01.RData")
save(NEFM_pos_glut_1v0, NEFM_pos_glut_6v0, file = "npglut_da_peaks_minpct0.01.RData")

for (i in 1) {
  peaks.use <- rownames(ATAC)[(str_detect(rownames(ATAC), chrs[11]))]
  obj <- ATAC[rownames(ATAC) %in% peaks.use, ]
  print(i)
  
  npglut_1v0 <- FindMarkers(
    object = obj,
    ident.1 = "NEFM_pos_glut_1hr",
    ident.2 = "NEFM_pos_glut_0hr",
    test.use = 'MAST',
    logfc.threshold = 0, 
    min.pct = 0.0, 
    random.seed = 100)
  
  npglut_6v0 <- FindMarkers(
    object = obj,
    ident.1 = "NEFM_pos_glut_6hr",
    ident.2 = "NEFM_pos_glut_0hr",
    test.use = 'MAST',
    logfc.threshold = 0,
    min.pct = 0.0,
    random.seed = 100)
} 

da_peaks_by_time_list_1v0 <- list(GABA_1v0, NEFM_neg_glut_1v0, NEFM_pos_glut_1v0)
da_peaks_by_time_list_6v0 <- list(GABA_6v0, NEFM_neg_glut_6v0, NEFM_pos_glut_6v0)

lists <- list(da_peaks_by_time_list_1v0, da_peaks_by_time_list_6v0)
types <- unique(str_remove(idents_list, "_[0|1|6]hr$"))
times <- c("1", "6")
for (j in 1:length(lists)) {
  time <- times[j]
  list <- lists[[j]]
  for (i in 1:length(list)){
    type <- types[i]
    ident <- paste(type, time, sep = "_")
    print(ident)
    #list[[i]] <- list[[i]][list[[i]]$p_val_adj < 0.05, ]
    uplist <- list[[i]][list[[i]]$avg_log2FC > 0, ]
    dolist <- list[[i]][list[[i]]$avg_log2FC < 0, ]
    df <- data.frame(peaks = rownames(uplist),
                     p_val = uplist$p_val,
                     avg_log2FC = uplist$avg_log2FC,
                     pct.1 = uplist$pct.1,
                     pct.2 = uplist$pct.2,
                     p_val_adj = uplist$p_val_adj)
    file_name <- paste0("./unfiltered_bed/",
                        ident, "v0hr_full_upregulated.bed")
    saveBed(df, file_name)
    df <- data.frame(peaks = rownames(dolist),
                     p_val = dolist$p_val,
                     avg_log2FC = dolist$avg_log2FC,
                     pct.1 = dolist$pct.1,
                     pct.2 = dolist$pct.2,
                     p_val_adj = dolist$p_val_adj)
    file_name <- paste0("./bed/",
                        ident, "v0hr_full_downregulated.bed")
    saveBed(df, file_name)
  }
}
