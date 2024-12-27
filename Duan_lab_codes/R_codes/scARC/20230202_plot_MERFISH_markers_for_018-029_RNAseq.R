# Chuxuan Li 02/02/2023
# Plot expression of marker genes from Fang et al. Science 2022 (MERFISH) of 
#018-029 RNAseq data 

# init ####
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(Seurat)
library(stringr)

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

objlist <- list(obj_018, obj_022, obj_024, obj_025, obj_029)
names(objlist) <- c("obj_018", "obj_022", "obj_024", "obj_025", "obj_029")

# add cluster x cell type ident ####
for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  obj <- subset(obj, cell.type != "unknown")
  obj <- subset(obj, cell.type != "unidentified")
  obj$subtypecluster <- ""
  clusters <- unique(obj$seurat_clusters)
  types <- sort(unique(obj$cell.type))
  for (t in types) {
    cat("type: ", t)
    cls <- unique(obj$seurat_clusters[obj$cell.type == t])
    
    for (j in 1:length(cls)) {
      cl <- cls[j]
      cat(" cluster: ", j, "\n")
      if (t == "npglut") {
        ty <- str_replace(t, "npglut", "NEFM+ glut")
      } else if (t == "nmglut") {
        ty <- str_replace(t, "nmglut", "NEFM- glut")
      } else if (t == "NEFM_pos_glut") {
        ty <- str_replace(t, "_pos_", "+ ")
      } else if (t == "NEFM_neg_glut") {
        ty <- str_replace(t, "_neg_", "- ")
      } else if (t == "SEMA3E_pos_GABA") {
        ty <- str_replace(t, "SEMA3E_pos_", "")
      } else {
        ty <- t
      }
      obj$subtypecluster[obj$seurat_clusters == cl & obj$cell.type == t] <- 
        paste0(ty, "_cluster ", j)
    }
  }
  print(unique(obj$subtypecluster))
  objlist[[i]] <- obj
}

# dotplot ####
MER_genes <- rev(c("CTGF", "PDGFRA", "OPALIN", "SELPLG", "AQP4", "SP8", "KLF5",
               "LGI2", "LAMP5", "GAD1", "TSHZ2", "CPLX3", "C1QL3", "SYT6", 
               "SMYD1", "TRABD2A", "COL21A1", "RORB", "CUX2", "SLC17A7", "SLC17A6"))
MER_genes <- as.factor(MER_genes)

for (i in 2:length(objlist)) {
  obj <- objlist[[i]]
  DefaultAssay(obj) <- "RNA"
  obj[["integrated"]] <- NULL
  # if (i != 1) {
  #   obj[["SCT"]] <- NULL
  # }
  png(filename = paste0("MERFISH_gex_dotplot_", names(objlist)[i], ".png"), 
      width = 750, height = 500)
  p <- DotPlot_new(obj, assay = "RNA", features = MER_genes, group.by = "subtypecluster") +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          axis.title = element_blank())
  print(p)
  dev.off()
}
# fix 018 separately
rowSums(obj_018@assays$RNA@counts[obj_018@assays$RNA@counts@Dimnames[[1]] %in% MER_genes, ]) /
  ncol(obj_018@assays$RNA@counts)
DefaultAssay(obj_018)
test <- CreateSeuratObject(counts = obj_018@assays$RNA@counts)
test$subtypecluster <- objlist[[1]]$subtypecluster
png(filename = paste0("MERFISH_gex_dotplot_", names(objlist)[1], ".png"), 
    width = 750, height = 500)
p <- DotPlot_new(test, assay = "RNA", features = MER_genes, group.by = "subtypecluster") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        axis.title = element_blank())
print(p)
dev.off()
