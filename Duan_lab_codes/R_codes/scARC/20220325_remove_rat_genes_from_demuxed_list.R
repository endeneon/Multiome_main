# Chuxuan Li 03/25/2022
# Remove rat genes from the datasets

# init ####
library(Seurat)
library(stringr)

load("demux_mapped_list.RData")

# remove nonhuman cells
for (i in 1:length(cleanobj_lst)){
  cleanobj_lst[[i]] <- subset(cleanobj_lst[[i]], cell.line.ident != "unmatched")
}
save(cleanobj_lst, file = "demux_removed_unmatched_list.RData")

# generate gene list to use ####
# select human genes using gencode.v30 gtf file 
gtf.file = "/data/Databases/Genomes/CellRanger_10x/hg38_rn6_hybrid/gencode.v30.annotation.gtf"
gtf.gr = rtracklayer::import(gtf.file) # creates a GRanges object
gtf.df = as.data.frame(gtf.gr)
human.genes = unique(gtf.df[ ,"gene_name"])

# import the 36601 genes from previous analysis
load("~/NVME/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/renormalized_4celltype_obj.RData")
prev.genes <- filtered_obj@assays$RNA@counts@Dimnames[[1]]
use.genes <- intersect(human.genes, prev.genes)


# subset gene list ####
all.genes <- list_18[[2]]@assays$RNA@counts@Dimnames[[1]]
all.genes[all.genes %in% use.genes]

finalobj_lst <- vector(mode = "list", length = length(list_18))
for(i in 1:length(list_18)){
  counts <- GetAssayData(list_18[[i]], assay = "RNA")
  counts <- counts[(which(rownames(counts) %in% use.genes)), ]
  finalobj_lst[[i]] <- subset(list_18[[i]], features = rownames(counts))
}

save(finalobj_lst, file = "removed_rat_genes_list.RData")
