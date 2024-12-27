# Chuxuan Li 06/16/2022
# Prepare raw ATACseq count matrices of 5, 18, 20 line data sets

# init ####
library(readxl)
library(stringr)
library(Seurat)
library(Signac)
library(future)
library(readr)

plan("multisession", workers = 2)

# import cell type object ATAC ####
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v4_combine_5_18_20line/5_18_20line_combined_labeled_nfeat5000_obj.RData")
unique(integrated_labeled$cell.type)
GABA <- subset(integrated_labeled, cell.type %in% c("GABA", "SST_pos_GABA"))
npglut <- subset(integrated_labeled, cell.type %in% c("NEFM_pos_glut", "NEFM_pos_glut?"))
nmglut <- subset(integrated_labeled, cell.type %in% c("NEFM_neg_glut"))
GABA_barcodes <- unlist(str_remove(GABA@assays$RNA@counts@Dimnames[[2]], "-[0-9]_[0-9]+$"))
npglut_barcodes <- unlist(str_remove(npglut@assays$RNA@counts@Dimnames[[2]], "-[0-9]_[0-9]+$"))
nmglut_barcodes <- unlist(str_remove(nmglut@assays$RNA@counts@Dimnames[[2]], "-[0-9]_[0-9]+$"))
rm(GABA, npglut, nmglut)


# 18 line import data ####
load("/nvmefs/scARC_Duan_018/Duan_project_024_ATAC/raw_obj_list.RData")
objlist_18line <- vector(mode = "list", length = length(pathlist))
tsvpathlist <- list.files(path = "/data/FASTQ/Duan_Project_024/hybrid_output/",
                          pattern = "atac_fragments.tsv.gz$", recursive = T,
                          full.names = T)
five = tsvpathlist[7:9]
seven = tsvpathlist[1:6]
nine = tsvpathlist[10:15]
tsvpathlist = c(five, seven, nine)

# create Seurat object from fragment files and h5 files
for (i in 1:length(objlist)){
  obj <- objlist[[i]]
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(obj),
    features = multiomic_obj@assays$peaks@ranges,
    cells = colnames(obj)
  )
  obj[["peaks"]] <- CreateChromatinAssay(
      counts = macs2_counts, 
      fragments = Fragments(obj))
  DefaultAssay(obj) <- "peaks"
  obj[["ATAC"]] <- NULL
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list_18line[i],
                                                  pattern = "-[0-6]"), "-"), "hr")
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

# Siwei 16 Jun 2022
# Set orig.ident as extracted from tsvpathlist
objlist_18line_sample_list <-
  str_extract(string = tsvpathlist,
              pattern = "[0-9]*\\-[0-9]")

for (i in 1:length(objlist_18line)) {
  objlist_18line[[i]]$orig.ident <- 
    objlist_18line_sample_list[i]
}
#####


for (i in 1:length(barcode_list)) {
  lines <- unique(barcode_list[[i]]$line)
  for (k in 1:length(objlist_18line)) {
    o <- unique(objlist_18line[[k]]$orig.ident)
    # print(o)
    refo <- unique(barcode_list[[i]]$orig)
    # print(refo)
    # objlist_18line[[k]]$cell.line.ident <- "unmatched"
    if (o == refo) break
  }
  # print(paste("k=", k))
  
  objlist_18line[[k]]$cell.line.ident <- "unmatched" #  here
  for (j in 1:length(lines)) {
    line_spec_barcodes <- barcode_list[[i]]$barcode[barcode_list[[i]]$line == lines[j]]
    objlist_18line[[k]]$cell.line.ident[objlist_18line[[k]]@assays$peaks@counts@Dimnames[[2]] %in% 
                                          line_spec_barcodes] <- lines[j]
    print(sum(objlist_18line[[k]]$cell.line.ident == "unmatched")) # should be around 1000-2000 
  }
}
save(objlist_18line, 
     file = "~/backup_space/Data_2_Lifan_15Jun2022_Lexi_Siwei/ATAC_18line_with_macs2_peaks_for_backup_only.RData")

# check unmatched cell count
for (i in 1:length(objlist_18line)) {
  print(sum(objlist_18line[[i]]$cell.line.ident == "unmatched")) 
}
## output from the above ranges from 2222 to 5762
## Should be acceptable

## check if all barcodes end with "-1", in both cases
for (i in 1:length(objlist_18line)) {
  print(length(objlist_18line[[i]]))
  # print(unique(str_split(string = colnames(objlist_18line[[i]]),
  #                        pattern = "-",
  #                        n = Inf, 
  #                        simplify = T)[, 2]))
}

for (i in 1:length(barcode_list)) {
  print(unique(str_split(string = barcode_list[[i]],
                         pattern = "-",
                         n = Inf, 
                         simplify = T)[, 2]))
}

# remove cells without barcode match ####
objlist_18line_clean <- objlist_18line
for (i in 1:length(objlist_18line)) {
  objlist_18line_clean[[i]] <- subset(objlist_18line[[i]], cell.line.ident != "unmatched")
}

# # match to RNA barcodes
# objlist_18line_RNA <- vector(mode = "list", length = length(h5list_18line))
# for (i in 1:length(objlist_18line)){
#   h5file <- Read10X_h5(filename = h5list_18line[i])
#   print(str_extract(string = h5list_18line[i],
#                     pattern = "[0-9]+-[0-6]"))
#   obj <- CreateSeuratObject(counts = h5file$`Gene Expression`,
#                             project = str_extract(string = h5list_18line[i],
#                                                   pattern = "[0-9]+-[0-6]"))
#   # assign library identity
#   obj$time.ident <- paste0(str_remove(str_extract(string = h5list_18line[i],
#                                                   pattern = "-[0-6]"), "-"), "hr")
#   print(unique(obj$time.ident))
#   # check number of genes and cells
#   print(paste0(i, " number of genes: ", nrow(obj),
#                ", number of cells: ", ncol(obj)))
#   objlist_18line_RNA[[i]] <- obj
# }
# for (i in 1:length(barcode_list)){
#   lines <- unique(barcode_list[[i]]$line)
#   for (k in 1:length(objlist_18line_RNA)){
#     o <- unique(objlist_18line_RNA[[k]]$orig.ident)
#     print(o)
#     refo <- unique(barcode_list[[i]]$orig)
#     print(refo)
#     if (o == refo) break
#   }
#   objlist_18line_RNA[[k]]$cell.line.ident <- "unmatched"
#   for (j in 1:length(lines)){
#     line_spec_barcodes <- barcode_list[[i]]$barcode[barcode_list[[i]]$line == lines[j]]
#     objlist_18line_RNA[[k]]$cell.line.ident[objlist_18line_RNA[[k]]@assays$RNA@counts@Dimnames[[2]] %in% 
#                                           line_spec_barcodes] <- lines[j]
#     print(sum(objlist_18line_RNA[[k]]$cell.line.ident == "unmatched"))
#   }
# }
# for (i in 1:length(objlist_18line_RNA)) {
#   objlist_18line_RNA[[i]] <- subset(objlist_18line_RNA[[i]], cell.line.ident != "unmatched")
# }
# 
# for (i in 1:length(objlist_18line)) {
#   objlist_18line[[i]]$RNA.match <- "unmatched"
#   objlist_18line[[i]]$RNA.match[objlist_18line[[i]]@assays$peaks@counts@Dimnames[[2]] %in% 
#                                   objlist_18line_RNA[[i]]@assays$RNA@counts@Dimnames[[2]]]<- "m"
#   print(sum(objlist_18line[[i]]$RNA.match == "unmatched"))
#   objlist_18line_clean[[i]] <- subset(objlist_18line[[i]], RNA.match != "unmatched")
# }

# assign cell type ####
for (i in 1:length(objlist_18line_clean)) {
  objlist_18line_clean[[i]]$cell.type <- "other"
  objlist_18line_clean[[i]]$trimmed_barcodes <- 
    str_remove(objlist_18line_clean[[i]]@assays$peaks@counts@Dimnames[[2]], "-[0-9]$")
  objlist_18line_clean[[i]]$cell.type[objlist_18line_clean[[i]]$trimmed_barcodes %in%
                                        GABA_barcodes] <- "GABA"
  objlist_18line_clean[[i]]$cell.type[objlist_18line_clean[[i]]$trimmed_barcodes %in%
                                        nmglut_barcodes] <- "nmglut"
  objlist_18line_clean[[i]]$cell.type[objlist_18line_clean[[i]]$trimmed_barcodes %in%
                                        npglut_barcodes] <- "npglut"
  print(sum(objlist_18line_clean[[i]]$cell.type == "GABA"))
  print(sum(objlist_18line_clean[[i]]$cell.type == "nmglut"))
  print(sum(objlist_18line_clean[[i]]$cell.type == "npglut"))
  print(sum(objlist_18line_clean[[i]]$cell.type == "other")) # check whether this number is unusually large; should be <1000
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
for (i in 1:length(cellline_time_sep_obj_lst_18)) {
  obj <- cellline_time_sep_obj_lst_18[[i]]
  obj <- subset(obj, rat.ident != "rat")
  cellline_time_sep_obj_lst_18[[i]] <- obj
}

# output files ####
setwd("~/backup_space/Data_2_Lifan_15Jun2022_Lexi_Siwei/count_matrices")
ATAC_obj_by_linextime_lst_18line <- cellline_time_sep_obj_lst_18
for (i in 1:length(ATAC_obj_by_linextime_lst_18line)) {
  obj <- ATAC_obj_by_linextime_lst_18line[[i]]
  obj$nCount_ATAC <- NULL
  obj$nFeature_ATAC <- NULL
  obj$rat.ident <- NULL
  obj$total_human_interval_reads <- NULL
  obj$total_rat_interval_reads <- NULL
  ATAC_obj_by_linextime_lst_18line[[i]] <- obj
}
save(ATAC_obj_by_linextime_lst_18line, file = "ATAC_count_matrices_by_time_and_cellline_18line.RData")


# 5 line ####
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_ATACseq_v2_call_peaks_by_celltypes/ATACseq_5line_objonly_labeled.RData")
DefaultAssay(obj_complete) <- "peaks"
obj_complete[["RNA"]] <- NULL
obj_complete[["ATAC"]] <- NULL
obj_complete[["SCT"]] <- NULL
setwd("/nvmefs/scARC_Duan_018/GRCh38_mapped_only/")
h5list_5line <- list.files(path = ".",
                           pattern = "filtered_feature_bc_matrix.h5",
                           recursive = T,
                           include.dirs = F, 
                           full.names = T)
tsvpathlist <- list.files(path = ".",
                          pattern = "atac_fragments.tsv.gz$", recursive = T,
                          full.names = T)

objlist_5line <- vector(mode = "list", length = length(h5list_5line))

# create Seurat object from fragment files and h5 files
for (i in 1:length(h5list_5line)){
  print(str_extract(string = h5list_5line[i],
                    pattern = "[0-9]+_[0-6]"))
  counts <- Read10X_h5(h5list_5line[[i]])
  counts <- counts$Peaks
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = tsvpathlist[i]
  )
  
  obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC", project = str_extract(string = h5list_5line[i],
                                          pattern = "[0-9]+_[0-6]"))
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(obj),
    features = obj_complete@assays$peaks@ranges,
    cells = colnames(obj)
  )
  obj[["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts, 
    fragments = tsvpathlist[i])
  DefaultAssay(obj) <- "peaks"
  obj[["ATAC"]] <- NULL
  obj$time.ident <- paste0(str_remove(str_extract(string = h5list_5line[i],
                                                  pattern = "_[0|1|6]"), "_"), "hr")
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
    obj$orig.ident <- str_replace_all(obj$orig.ident, "_", "-")
    o <- unique(obj$orig.ident)
    refo <- str_replace(str_remove(str_extract(pathlist[j], "g_[2|8]_[0|1|6]"), "^g_"), "_", "-")
    print(paste(o, refo))
    print(bcnames[j])
    if (refo %in% o) {
      obj$cell.line.ident[str_remove(obj@assays$peaks@counts@Dimnames[[2]], "-[0-9]+$") %in% 
                            bclist[[j]]] <- bcnames[j]
    }
  }
  print(unique(obj$cell.line.ident))
  objlist_5line[[i]] <- obj
}
save(objlist_5line, 
     file = "~/backup_space/Data_2_Lifan_15Jun2022_Lexi_Siwei/ATAC_5line_with_macs2_peaks_for_backup_only.RData")

objlist_5line_clean <- objlist_5line
for (i in 1:length(objlist_5line)) {
  print(sum(objlist_5line[[i]]$cell.line.ident == "unmatched"))
  objlist_5line_clean[[i]] <- subset(objlist_5line[[i]], cell.line.ident != "unmatched")
}

# assign cell type ####
for (i in 1:length(objlist_5line_clean)) {
  objlist_5line_clean[[i]]$cell.type <- "other"
  objlist_5line_clean[[i]]$trimmed_barcodes <- 
    str_remove(objlist_5line_clean[[i]]@assays$peaks@counts@Dimnames[[2]], "-[0-9]+$")
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
names(cellline_time_sep_obj_lst_5) <- obj_names

# output files ####
setwd("~/backup_space/Data_2_Lifan_15Jun2022_Lexi_Siwei/count_matrices")
ATAC_obj_by_linextime_lst_5line <- cellline_time_sep_obj_lst_5
for (i in 1:length(ATAC_obj_by_linextime_lst_5line)) {
  obj <- ATAC_obj_by_linextime_lst_5line[[i]]
  obj$nCount_ATAC <- NULL
  obj$nFeature_ATAC <- NULL
  ATAC_obj_by_linextime_lst_5line[[i]] <- obj
}
save(ATAC_obj_by_linextime_lst_5line, file = "ATAC_count_matrices_by_time_and_cellline_5line.RData")
