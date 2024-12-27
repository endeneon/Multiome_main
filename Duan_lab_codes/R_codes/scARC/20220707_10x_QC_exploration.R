# Chuxuan Li 07/07/2022
# determine the QC used by 10x Genomics

#load data
load("/nvmefs/scARC_Duan_018/Duan_project_024_RNA/Analysis_v1_normalize_by_lib_then_integrate/raw_data_list.RData")

# check min # of genes that have nonzero reads
sum(rowSums(objlist[[1]]@assays$RNA@counts) != 0)
sum(colSums(objlist[[1]]@assays$RNA@counts) != 0)
min(colSums(objlist[[1]]@assays$RNA@counts))
min(colSums(objlist[[2]]@assays$RNA@counts))
max(colSums(objlist[[1]]@assays$RNA@counts))
max(colSums(objlist[[2]]@assays$RNA@counts))
# check percentile of values in the array
percentile <- function(value, sample){
  s = as.vector(sort(sample, decreasing = F))
  ind = which(s == value)
  return(ind)}

for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  # print("# of cell with nonzero counts")
  # print(sum(colSums(obj@assays$RNA@counts) != 0))
  # print("# of genes with nonzero counts")
  # print(sum(rowSums(obj@assays$RNA@counts) != 0))
  # print("min")
  # print(min(colSums(obj@assays$RNA@counts)))
  # print("max")
  # print(max(colSums(obj@assays$RNA@counts)))
  
  min <- min(colSums(obj@assays$RNA@counts))
}
for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  # print("# of cell with nonzero counts")
  # print(sum(colSums(obj@assays$RNA@counts) != 0))
  # print("# of genes with nonzero counts")
  # print(sum(rowSums(obj@assays$RNA@counts) != 0))
  # print("min")
  # print(min(colSums(obj@assays$RNA@counts)))
  # print("max")
  # print(max(colSums(obj@assays$RNA@counts)))
  
  min <- min(colSums(obj@assays$RNA@counts))
}

# count % reads in mitochondrial genes ####
pctmt <- rep_len(0, length(objlist))
for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  mtgenes <- rownames(obj@assays$RNA)[str_detect(rownames(obj@assays$RNA), "^MT-*")]
  mtcounts <- obj@assays$RNA@counts[mtgenes, ]
  counts <- obj@assays$RNA@counts
  print(max(sum(mtcounts@x)/sum(counts@x)))
  pctmt[i] <- max(sum(mtcounts@x)/sum(counts@x))
}
median(pctmt)

# find RNA mapped reads min and max ####
library(readr)
metrics_files <- list.files("/data/FASTQ/Duan_Project_024/hybrid_output/", "per_barcode_metrics.csv", 
           full.names = T,
           recursive = T)
mettables <- vector("list", length(metrics_files))
filtered_tables <- vector("list", length(metrics_files))
for (i in 1:length(metrics_files)) {
  mettables[[i]] <- read_csv(metrics_files[i])
  filtered_tables[[i]] <- mettables[[i]][mettables[[i]]$is_cell == 1, ]
}

readmins <- rep_len(0, length(metrics_files))
readmaxs <- rep_len(0, length(metrics_files))

for (i in 1:length(filtered_tables)) {
  f <- filtered_tables[[i]]
  nreads <- f$gex_mapped_reads
  readmaxs[i] <- max(nreads)
  readmins[i] <- min(nreads)
}
median(readmaxs)
median(readmins)

# find RNA mapped reads min and max ####
fragmins <- rep_len(0, length(metrics_files))
fragmaxs <- rep_len(0, length(metrics_files))

for (i in 1:length(filtered_tables)) {
  f <- filtered_tables[[i]]
  nfrags <- f$atac_fragments
  fragmaxs[i] <- max(nfrags)
  fragmins[i] <- min(nfrags)
}
median(fragmaxs)
median(fragmins)
