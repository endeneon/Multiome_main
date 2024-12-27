# Chuxuan Li 06/13/2022
# Generate data for Lifan, including batch information, cell line covariate
#information, and the count matrices for each cell line

# init ####
library(readxl)
library(stringr)
library(Seurat)

MGS <- read_excel("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/MGS_iPSC lines_60samples_scRNA_ATAC-seq_bulkATAC-seq status_AK_HZ.xlsx")
load("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v2_combine_5line_18line/5line_18line_combined_labeled_nfeat5000_obj.RData")

lines <- unique(integrated_labeled$cell.line.ident)
rm(integrated_labeled)

load("~/NVME/scARC_Duan_018/Duan_project_024_RNA/Analysis_v3_combine_18line_20line/18line_20line_combined_labeled_nfeat5000_obj.RData")
used_lines <- union(lines, unique(integrated_labeled$cell.line.ident))
rm(integrated_labeled)


# clean data table ####
covar_table <- as.data.frame(cbind(MGS$`new iPSC line ID`, MGS$`RUID cell line ID`, MGS$ID,
                                   MGS$Aff, MGS$`sex (1=M, 2=F)`, MGS$Age, MGS$`Co-culture batch`))
colnames(covar_table) <- c("cell_line", "RUID", "ID", "disease_status", "sex", "age", "batch")
covar_table$cell_line <- str_replace(covar_table$cell_line, "CD00000", "CD_")
rownames(covar_table) <- covar_table$cell_line

covar_table_clean <- covar_table[covar_table$cell_line %in% used_lines, ]
covar_table_clean$batch <- str_replace_all(covar_table_clean$batch, " & 13", "")
covar_table_clean$batch <- str_replace_all(covar_table_clean$batch, "1 & ", "")
covar_table_clean$sex[covar_table_clean$sex == "1"] <- "M"
covar_table_clean$sex[covar_table_clean$sex == "2"] <- "F"


# import count matrices - RNA ####
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v4_combine_5_18_20line/5_18_20line_combined_labeled_nfeat5000_obj.RData")
unique(integrated_labeled$cell.type)
GABA <- subset(integrated_labeled, cell.type %in% c("GABA", "SST_pos_GABA"))
npglut <- subset(integrated_labeled, cell.type %in% c("NEFM_pos_glut", "NEFM_pos_glut?"))
nmglut <- subset(integrated_labeled, cell.type %in% c("NEFM_neg_glut"))
GABA_barcodes <- unlist(str_remove(GABA@assays$RNA@counts@Dimnames[[2]], "-[0-9]_[0-9]+$"))
npglut_barcodes <- unlist(str_remove(npglut@assays$RNA@counts@Dimnames[[2]], "-[0-9]_[0-9]+$"))
nmglut_barcodes <- unlist(str_remove(nmglut@assays$RNA@counts@Dimnames[[2]], "-[0-9]_[0-9]+$"))
rm(GABA, npglut, nmglut)

for (i in unique(integrated_labeled$orig.ident)) {
  for (j in unique(integrated_labeled$cell.type)) {
    print(i)
    print(j)
    print(sum(integrated_labeled$orig.ident == i & integrated_labeled$cell.type == j))
  }
}

# 18 line import data ####

setwd("/data/FASTQ/Duan_Project_024/hybrid_output/")
# note this h5 file contains both atac-seq and gex information
# Read files, separate the ATACseq data from the .h5 matrix
h5list_18line <- list.files(path = ".",
                     pattern = "filtered_feature_bc_matrix.h5",
                     recursive = T,
                     include.dirs = F, 
                     full.names = T)
h5list_18line <- h5list_18line[1:15]
h5list_18line <- sort(h5list_18line)
objlist_18line <- vector(mode = "list", length = length(h5list_18line))

for (i in 1:length(objlist_18line)){
  h5file <- Read10X_h5(filename = h5list_18line[i])
  print(str_extract(string = h5list_18line[i],
                    pattern = "[0-9]+-[0-6]"))
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_extract(string = h5list_18line[i],
                                                  pattern = "[0-9]+-[0-6]"))
  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list_18line[i],
                                                  pattern = "-[0-6]"), "-"), "hr")
  print(unique(obj$time.ident))
  # check number of genes and cells
  print(paste0(i, " number of genes: ", nrow(obj),
               ", number of cells: ", ncol(obj)))
  objlist_18line[[i]] <- obj
}

# read 18 line specific barcodes ####
setwd("/data/FASTQ/Duan_Project_024/hybrid_output")
pathlist <- list.files(path = ".", full.names = T, pattern = ".best.tsv", recursive = T)
pathlist <- pathlist[str_detect(pathlist, "demux_0.01", T)]
barcode_list <- vector(mode = "list", length = length(pathlist))
for (i in 1:length(pathlist)){
  barcode_list[[i]] <- read.delim(pathlist[i], header = F, row.names = NULL)
  colnames(barcode_list[[i]]) <- c("barcode", "line")
}

# assign the origin of the barcodes for matching later ####
for (i in 1:length(pathlist)) {
  g = str_remove(str_extract(pathlist[i], "group_[0-9]+"), "group_")
  g = str_remove(g, "^0")
  t = str_remove(str_remove_all(str_extract(pathlist[i], ".[0|1|2].best"), "\\."), "best")
  if (t == "2") {
    t = "6"
  }
  orig = paste(g, t, sep = "-")
  print(orig)
  
  barcode_list[[i]]$orig <- orig
}

for (i in 1:length(barcode_list)){
  lines <- unique(barcode_list[[i]]$line)
  for (k in 1:length(objlist_18line)){
    o <- unique(objlist_18line[[k]]$orig.ident)
    print(o)
    refo <- unique(barcode_list[[i]]$orig)
    print(refo)
    if (o == refo) break
  }
  objlist_18line[[k]]$cell.line.ident <- "unmatched"
  for (j in 1:length(lines)){
    line_spec_barcodes <- barcode_list[[i]]$barcode[barcode_list[[i]]$line == lines[j]]
    objlist_18line[[k]]$cell.line.ident[objlist_18line[[k]]@assays$RNA@counts@Dimnames[[2]] %in% 
                                        line_spec_barcodes] <- lines[j]
    print(sum(objlist_18line[[k]]$cell.line.ident == "unmatched"))
  }
}
# remove cells without barcode match ####
objlist_18line_clean <- objlist_18line
for (i in 1:length(objlist_18line)) {
  objlist_18line_clean[[i]] <- subset(objlist_18line[[i]], cell.line.ident != "unmatched")
}


# assign cell type ####
for (i in 1:length(objlist_18line_clean)) {
  objlist_18line_clean[[i]]$cell.type <- "other"
  objlist_18line_clean[[i]]$trimmed_barcodes <- 
    str_remove(objlist_18line_clean[[i]]@assays$RNA@counts@Dimnames[[2]], "-[0-9]$")
  objlist_18line_clean[[i]]$cell.type[objlist_18line_clean[[i]]$trimmed_barcodes %in%
                                        GABA_barcodes] <- "GABA"
  objlist_18line_clean[[i]]$cell.type[objlist_18line_clean[[i]]$trimmed_barcodes %in%
                                        nmglut_barcodes] <- "nmglut"
  objlist_18line_clean[[i]]$cell.type[objlist_18line_clean[[i]]$trimmed_barcodes %in%
                                        npglut_barcodes] <- "npglut"
  print(sum(objlist_18line_clean[[i]]$cell.type == "GABA"))
  print(sum(objlist_18line_clean[[i]]$cell.type == "nmglut"))
  print(sum(objlist_18line_clean[[i]]$cell.type == "npglut"))
  print(sum(objlist_18line_clean[[i]]$cell.type == "other"))
}

# separate into cell line x time points
cellline_time_sep_obj_lst_18 <- vector("list", 18*3)
k = 1
obj_names <- rep_len("", 18*3)
for (i in 1:length(objlist_18line_clean)) {
  lines = sort(unique(objlist_18line_clean[[i]]$cell.line.ident))
  t = unique(objlist_18line_clean[[i]]$time.ident)
  for (j in 1:length(lines)) {
    print(lines[j])
    subobj <- subset(objlist_18line_clean[[i]], cell.line.ident == lines[j])
    cellline_time_sep_obj_lst_18[[k]] <- subobj
    obj_names[k] <- paste(lines[j], t, sep = "-") 
    k = k + 1
  }
}
names(cellline_time_sep_obj_lst_18) <- obj_names

# 20 line ####
setwd("/data/FASTQ/Duan_Project_022_Reseq")
# note this h5 file contains both atac-seq and gex information
# Read files, separate the ATACseq data from the .h5 matrix
h5list_20line <- list.files(path = ".",
                            pattern = "filtered_feature_bc_matrix.h5",
                            recursive = T,
                            include.dirs = F, 
                            full.names = T)
h5list_20line <- h5list_20line[str_detect(h5list_20line, "hybrid", T)]
h5list_20line <- h5list_20line[2:16]
objlist_20line <- vector(mode = "list", length = length(h5list_20line))

for (i in 1:length(objlist_20line)){
  h5file <- Read10X_h5(filename = h5list_20line[i])
  print(str_replace(str_extract(string = h5list_20line[i],
                                pattern = "[0-9]+_[0-6]"), "_", "-"))
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_replace(str_extract(string = h5list_20line[i],
                                                  pattern = "[0-9]+_[0-6]"), "_", "-"))
  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list_20line[i],
                                                  pattern = "_[0|1|6]"), "_"), "hr")
  print(unique(obj$time.ident))
  # check number of genes and cells
  print(paste0(i, " number of genes: ", nrow(obj),
               ", number of cells: ", ncol(obj)))
  objlist_20line[[i]] <- obj
}


cdlist <- data.frame(lib22 = c("CD_22", "CD_23", "CD_37", "CD_40"),
                     lib36 = c("CD_36", "CD_38", "CD_43", "CD_45"),
                     lib39 = c("CD_39", "CD_42", "CD_47", "CD_48"),
                     lib44 = c("CD_31", "CD_44", "CD_46", "CD_57"),
                     lib49 = c("CD_32", "CD_49",  "CD_58", "CD_61"))

id <- c("22_0", "22_1", "22_6",
        "36_0", "36_1", "36_6", 
        "39_0", "39_1", "39_6", 
        "44_0", "44_1", "44_6",
        "49_0", "49_1","49_6")

for (i in 1:5){
  # read in barcodes for each lib
  for(j in 1:3){
    print(id[i * 3 - 3 + j])
    df <-
      read.delim(paste0("./demux_barcode_output/output_barcodes/libraries_",
                        id[i * 3 - 3 + j],
                        "_SNG_by_cell_line.tsv"),
                 quote="")
    obj <- objlist_20line[[i * 3 - 3 + j]]
    print(unique(obj$orig.ident))
    
    # assign cell line looping over 4 cell lines per lib
    obj$cell.line.ident <- "unmatched"
    for (k in 1:4){
      # assign cell line in each lib
      ident <- cdlist[k, i]
      print(ident)
      subdf <- df$BARCODE[df$SNG.1ST == ident]
      obj$cell.line.ident[obj@assays$RNA@counts@Dimnames[[2]] %in% subdf] <- ident
    }
    print(unique(obj$cell.line.ident))
    objlist_20line[[i * 3 - 3 + j]] <- obj
  }
  
}

# remove cells without barcode match ####
objlist_20line_clean <- objlist_20line
for (i in 1:length(objlist_20line)) {
  print(sum(objlist_20line[[i]]$cell.line.ident == "unmatched"))
  objlist_20line_clean[[i]] <- subset(objlist_20line[[i]], cell.line.ident != "unmatched")
}


# assign cell type ####
for (i in 1:length(objlist_20line_clean)) {
  objlist_20line_clean[[i]]$cell.type <- "other"
  objlist_20line_clean[[i]]$trimmed_barcodes <- 
    str_remove(objlist_20line_clean[[i]]@assays$RNA@counts@Dimnames[[2]], "-[0-9]+$")
  objlist_20line_clean[[i]]$cell.type[objlist_20line_clean[[i]]$trimmed_barcodes %in%
                                        GABA_barcodes] <- "GABA"
  objlist_20line_clean[[i]]$cell.type[objlist_20line_clean[[i]]$trimmed_barcodes %in%
                                        nmglut_barcodes] <- "nmglut"
  objlist_20line_clean[[i]]$cell.type[objlist_20line_clean[[i]]$trimmed_barcodes %in%
                                        npglut_barcodes] <- "npglut"
  print(sum(objlist_20line_clean[[i]]$cell.type == "GABA"))
  print(sum(objlist_20line_clean[[i]]$cell.type == "nmglut"))
  print(sum(objlist_20line_clean[[i]]$cell.type == "npglut"))
  print(sum(objlist_20line_clean[[i]]$cell.type == "other"))
}

# separate into cell line x time points ####
cellline_time_sep_obj_lst_20 <- vector("list", 20*3)
k = 1
obj_names <- rep_len("", 20*3)
for (i in 1:length(objlist_20line_clean)) {
  lines = sort(unique(objlist_20line_clean[[i]]$cell.line.ident))
  t = unique(objlist_20line_clean[[i]]$time.ident)
  print(t)
  for (j in 1:length(lines)) {
    print(lines[j])
    subobj <- subset(objlist_20line_clean[[i]], cell.line.ident == lines[j])
    cellline_time_sep_obj_lst_20[[k]] <- subobj
    obj_names[k] <- paste(lines[j], t, sep = "-") 
    k = k + 1
  }
}
names(cellline_time_sep_obj_lst_20) <- obj_names


# 5 line ####
setwd("/nvmefs/scARC_Duan_018/GRCh38_mapped_only/")
h5list_5line <- list.files(path = ".",
                            pattern = "filtered_feature_bc_matrix.h5",
                            recursive = T,
                            include.dirs = F, 
                            full.names = T)
objlist_5line <- vector(mode = "list", length = length(h5list_5line))

for (i in 1:length(objlist_5line)){
  h5file <- Read10X_h5(filename = h5list_5line[i])
  print(str_replace(str_extract(string = h5list_5line[i],
                                pattern = "[0-9]+_[0-6]"), "_", "-"))
  obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
                            project = str_replace(str_extract(string = h5list_5line[i],
                                                              pattern = "[0-9]+_[0-6]"), "_", "-"))
  # assign library identity
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list_5line[i],
                                                  pattern = "_[0|1|6]"), "_"), "hr")
  print(unique(obj$time.ident))
  # check number of genes and cells
  print(paste0(i, " number of genes: ", nrow(obj),
               ", number of cells: ", ncol(obj)))
  objlist_5line[[i]] <- obj
}

# assign cell line to 5line ####
pathlist <- list.files("/nvmefs/scARC_Duan_018/GRCh38_mapped_only/common_barcodes",
                      pattern = "*.best.", full.names = T, recursive = T)
bclist <- vector("list", length(pathlist))
bcnames <- rep_len("", length(pathlist))
for (i in 1:length(pathlist)) {
  bc <- read_csv(pathlist[i], col_names = FALSE)
  bc <- unlist(bc)
  bc <- str_remove(bc, "-[0-9]+$")
  bclist[[i]] <- bc
  bcnames[i] <- str_extract(pathlist[i], "CD_[0-9][0-9]$")
}

# assign cell line and remove unmatched in 5 line ####

for (i in 1:length(objlist_5line)){
  obj <- objlist_5line[[i]]
  obj$cell.line.ident <- "unmatched"
  for (j in 1:length(bclist)) {
    o <- unique(obj$orig.ident)
    refo <- str_replace(str_remove(str_extract(pathlist[j], "g_[2|8]_[0|1|6]"), "^g_"), "_", "-")
    print(paste(o, refo))
    print(bcnames[j])
    if (refo %in% o) {
      obj$cell.line.ident[str_remove(obj@assays$RNA@counts@Dimnames[[2]], "-[0-9]+$") %in% 
                            bclist[[j]]] <- bcnames[j]
    }
  }
  print(unique(obj$cell.line.ident))
  objlist_5line[[i]] <- obj
}

objlist_5line_clean <- objlist_5line
for (i in 1:length(objlist_5line)) {
  print(sum(objlist_5line[[i]]$cell.line.ident == "unmatched"))
  objlist_5line_clean[[i]] <- subset(objlist_5line[[i]], cell.line.ident != "unmatched")
}

# assign cell type ####
for (i in 1:length(objlist_5line_clean)) {
  objlist_5line_clean[[i]]$cell.type <- "other"
  objlist_5line_clean[[i]]$trimmed_barcodes <- 
    str_remove(objlist_5line_clean[[i]]@assays$RNA@counts@Dimnames[[2]], "-[0-9]+$")
  objlist_5line_clean[[i]]$cell.type[objlist_5line_clean[[i]]$trimmed_barcodes %in%
                                        GABA_barcodes] <- "GABA"
  objlist_5line_clean[[i]]$cell.type[objlist_5line_clean[[i]]$trimmed_barcodes %in%
                                        nmglut_barcodes] <- "nmglut"
  objlist_5line_clean[[i]]$cell.type[objlist_5line_clean[[i]]$trimmed_barcodes %in%
                                        npglut_barcodes] <- "npglut"
  print(sum(objlist_5line_clean[[i]]$cell.type == "GABA"))
  print(sum(objlist_5line_clean[[i]]$cell.type == "nmglut"))
  print(sum(objlist_5line_clean[[i]]$cell.type == "npglut"))
  print(sum(objlist_5line_clean[[i]]$cell.type == "other"))
}

# separate into cell line x time points ####
cellline_time_sep_obj_lst_5 <- vector("list", 5*3)
k = 1
obj_names <- rep_len("",5*3)
for (i in 1:length(objlist_5line_clean)) {
  lines = sort(unique(objlist_5line_clean[[i]]$cell.line.ident))
  t = unique(objlist_5line_clean[[i]]$time.ident)
  for (j in 1:length(lines)) {
    print(lines[j])
    subobj <- subset(objlist_5line_clean[[i]], cell.line.ident == lines[j])
    cellline_time_sep_obj_lst_5[[k]] <- subobj
    obj_names[k] <- paste(lines[j], t, sep = "-") 
    k = k + 1
  }
}
names(cellline_time_sep_obj_lst_20) <- obj_names


# combined ####
# remove rat genes for all objs ####
# import the 36601 genes from previous analysis
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/renormalized_4celltype_obj.RData")
prev.genes <- filtered_obj@assays$RNA@counts@Dimnames[[1]]

nlib <- length(cellline_time_sep_obj_lst_5) + 
  length(cellline_time_sep_obj_lst) + length(cellline_time_sep_obj_lst_20)
obj_by_linextime_lst <- vector(mode = "list", length = nlib)

j = 1
k = 1
for(i in 1:nlib){
  if (i <= length(cellline_time_sep_obj_lst_5)) {
    counts <- GetAssayData(cellline_time_sep_obj_lst_5[[i]], assay = "RNA")
    counts <- counts[(which(rownames(counts) %in% prev.genes)), ]
    obj_by_linextime_lst[[i]] <- subset(cellline_time_sep_obj_lst_5[[i]], features = rownames(counts))
  } else if (i <= (length(cellline_time_sep_obj_lst_5) + 
                   length(cellline_time_sep_obj_lst_18))) {
    counts <- GetAssayData(cellline_time_sep_obj_lst_18[[j]], assay = "RNA")
    counts <- counts[(which(rownames(counts) %in% prev.genes)), ]
    obj_by_linextime_lst[[i]] <- subset(cellline_time_sep_obj_lst_18[[j]], features = rownames(counts))
    j = j + 1
  } else {
    counts <- GetAssayData(cellline_time_sep_obj_lst_20[[k]], assay = "RNA")
    counts <- counts[(which(rownames(counts) %in% prev.genes)), ]
    obj_by_linextime_lst[[i]] <- subset(cellline_time_sep_obj_lst_20[[k]], features = rownames(counts))
    k = k + 1
  }
}

use.genes <- obj_by_linextime_lst[[16]]@assays$RNA@counts@Dimnames[[1]]

for(i in 1:length(obj_by_linextime_lst)){
  counts <- GetAssayData(obj_by_linextime_lst[[i]], assay = "RNA")
  counts <- counts[(which(rownames(counts) %in% use.genes)), ]
  obj_by_linextime_lst[[i]] <- subset(obj_by_linextime_lst[[i]], features = rownames(counts))
}

# cell type proportions ####
med_df <- data.frame(cell_line = rep_len("", length(obj_by_linextime_lst)),
                     time = rep_len("", length(obj_by_linextime_lst)),
                     GABA_counts = rep_len("", length(obj_by_linextime_lst)),
                     GABA_fraction = rep_len("", length(obj_by_linextime_lst)),
                     nmglut_counts = rep_len("", length(obj_by_linextime_lst)),
                     nmglut_fraction = rep_len("", length(obj_by_linextime_lst)),
                     npglut_counts = rep_len("", length(obj_by_linextime_lst)),
                     npglut_fraction = rep_len("", length(obj_by_linextime_lst)),
                     total_counts = rep_len("", length(obj_by_linextime_lst)))
for (i in 1:length(obj_by_linextime_lst)) {
  obj <- obj_by_linextime_lst[[i]]
  line <- unique(obj$cell.line.ident)
  time <- unique(obj$time.ident)
  types <- unique(obj$cell.type)
  print(line)
  print(time)
  med_df$cell_line[i] <- line
  med_df$time[i] <- time
  med_df$GABA_counts[i] <- sum(obj$cell.type == "GABA")
  med_df$nmglut_counts[i] <- sum(obj$cell.type == "nmglut")
  med_df$npglut_counts[i] <- sum(obj$cell.type == "npglut")
  med_df$total_counts[i] <- ncol(obj)
  med_df$GABA_fraction[i] <- sum(obj$cell.type == "GABA") / ncol(obj)
  med_df$nmglut_fraction[i] <- sum(obj$cell.type == "nmglut") / ncol(obj)
  med_df$npglut_fraction[i] <- sum(obj$cell.type == "npglut") / ncol(obj)
}

# combine into covariate table
distributeTime <- function(df) {
  for (i in 1:3) {
    t = df$time[i]
    if (i == 1) {
      df_out = df[df$time == t, 2:ncol(df)]
      colnames(df_out) = paste(str_split(colnames(df_out), pattern = "_", n = 2, T)[, 1],
                               t, str_split(colnames(df_out), pattern = "_", n = 2, T)[, 2],
                               sep = "_")
    } else {
      df_tobind = df[df$time == t, 2:ncol(df)]
      colnames(df_tobind) = paste(str_split(colnames(df_tobind), pattern = "_", n = 2, T)[, 1],
                               t, str_split(colnames(df_tobind), pattern = "_", n = 2, T)[, 2],
                               sep = "_")
      df_out = cbind(df_out, df_tobind)
    }
  }
  return(df_out)
}

for (i in 1:length(covar_table_clean$cell_line)) {
  line <- covar_table_clean$cell_line[i]
  print(line)
  subdf <- med_df[med_df$cell_line == line, 2:ncol(med_df)]
  if (i == 1) {
    outdf <- distributeTime(subdf)
  } else {
    df_toappend <- distributeTime(subdf)
    outdf <- rbind(outdf, df_toappend)
  }
}
covar_table_final <- cbind(covar_table_clean, outdf)

# output data ####
RNA_obj_by_linextime_lst <- obj_by_linextime_lst
write.table(covar_table_final, 
            file = "~/backup_space/Data_2_Lifan_15Jun2022_Lexi_Siwei/covariates_with_matching_cell_line_id.csv",
            quote = F, sep = ",", row.names = F)
save(RNA_obj_by_linextime_lst, file = "~/backup_space/Data_2_Lifan_15Jun2022_Lexi_Siwei/RNA_count_matrices_by_time_and_cellline.RData")
