# Chuxuan Li 02/08/2023
# Use data mapped to rat + rabbit, separate rat and rabbit after removing cells
#mapped to human barcodes

# init ####
{
  library(Seurat)
  library(glmGamPoi)
  library(sctransform)
  
  library(ggplot2)
  library(RColorBrewer)
  library(viridis)
  library(graphics)
  library(ggrepel)
  
  library(dplyr)
  library(readr)
  library(stringr)
  library(readxl)
  
  library(future)
}
# set threads and parallelization
plan("multisession", workers = 1)
# plan("sequential")
# plan()
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)

setwd("/nvmefs/scARC_Duan_018/Duan_project_029_RNA/oc2_mm10_analysis")
load("./removed_human_cells_objlist.RData")

# calculate percentage of reads mapped to each genome ####
for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  all_genes <- obj@assays$RNA@counts@Dimnames[[1]]
  rabbit_genes_ind <- str_detect(all_genes, "^[A-Z][A-Z0-9]+")
  mouse_genes_ind <- !rabbit_genes_ind[str_detect(all_genes, "5S-rRNA*", T)]
  obj$nreads_mapped_to_rabbit <- colSums(obj@assays$RNA@counts[rabbit_genes_ind, ])
  obj$nreads_mapped_to_mouse <- colSums(obj@assays$RNA@counts[mouse_genes_ind, ])
  obj$pct_reads_mapped_to_rabbit <- obj$nreads_mapped_to_rabbit / colSums(obj@assays$RNA@counts)
  obj$pct_reads_mapped_to_mouse <- obj$nreads_mapped_to_mouse / colSums(obj@assays$RNA@counts)
  obj$mouse_to_rabbit_ratio <- obj$nreads_mapped_to_mouse / obj$nreads_mapped_to_rabbit 
  objlist[[i]] <- obj
}

for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  library <- unique(obj$orig.ident)
  png(paste0("./histograms/nreads_mapped_to_mouse_hist_", library, ".png"),
      width = 800, height = 400)
  p <- hist(obj$nreads_mapped_to_mouse, breaks = 200, 
       main = paste0("distribution of # of reads mapped to mouse genome - ", library), 
       xlab = "# reads mapped to mouse genome")
  print(p)
  dev.off()
  png(paste0("./histograms/nreads_mapped_to_rabbit_hist_", library, ".png"),
      width = 800, height = 400)
  p <- hist(obj$nreads_mapped_to_rabbit, breaks = 200, 
       main = paste0("distribution of # of reads mapped to rabbit genome - ", library), 
       xlab = "# reads mapped to rabbit genome")
  print(p)
  dev.off()
  
  png(paste0("./histograms/pct_reads_mapped_to_mouse_hist_", library, ".png"),
      width = 800, height = 400)
  p <- hist(obj$pct_reads_mapped_to_mouse, breaks = 200, 
       main = paste0("distribution of percent of reads mapped to mouse genome - ", library), 
       xlab = "% reads mapped to mouse genome")
  print(p)
  dev.off()
  
  png(paste0("./histograms/pct_reads_mapped_to_rabbit_hist_", library, ".png"),
      width = 800, height = 400)
  p <- hist(obj$pct_reads_mapped_to_rabbit, breaks = 200, 
       main = paste0("distribution of percent of reads mapped to rabbit genome - ", library), 
       xlab = "% reads mapped to rabbit genome")
  print(p)
  dev.off()
  
  png(paste0("./histograms/m2r_ratio_hist_", library, ".png"),
      width = 800, height = 400)
  p <- hist(obj$mouse_to_rabbit_ratio, breaks = 200, 
       main = "#reads mapped to mouse : #reads mapped to rabbit", 
       xlab = "#reads mapped to mouse : #reads mapped to rabbit")
  print(p)
  dev.off()
}


# map identity back to clustered objects ####
load("labeled_by_species_obj_list.RData")
for (i in 1:length(scaled_lst)) {
  scaled_obj <- scaled_lst[[i]]
  obj <- objlist[[i]]
  mouse_bc <- obj@assays$RNA@counts@Dimnames[[2]][obj$pct_reads_mapped_to_mouse > 0.7]
  rabbit_bc <- obj@assays$RNA@counts@Dimnames[[2]][obj$pct_reads_mapped_to_rabbit > 0.7]
  scaled_obj$species.ident.new <- "unknown"
  scaled_obj$species.ident.new[scaled_obj@assays$RNA@counts@Dimnames[[2]] %in% mouse_bc] <- "mouse"
  scaled_obj$species.ident.new[scaled_obj@assays$RNA@counts@Dimnames[[2]] %in% rabbit_bc] <- "rabbit"
  scaled_obj$species.ident.by.ratio <- "unknown"
  mouse_bc <- obj@assays$RNA@counts@Dimnames[[2]][obj$mouse_to_rabbit_ratio > 25]
  rabbit_bc <- obj@assays$RNA@counts@Dimnames[[2]][obj$mouse_to_rabbit_ratio < 25]
  scaled_obj$species.ident.by.ratio[scaled_obj@assays$RNA@counts@Dimnames[[2]] %in% mouse_bc] <- "mouse"
  scaled_obj$species.ident.by.ratio[scaled_obj@assays$RNA@counts@Dimnames[[2]] %in% rabbit_bc] <- "rabbit"
  scaled_lst[[i]] <- scaled_obj
  scaled_obj$mouse.to.rabbit.ratio <- NA
  scaled_obj$mouse.to.rabbit.ratio[scaled_obj@assays$RNA@counts@Dimnames[[2]] %in% 
                                     obj@assays$RNA@counts@Dimnames[[2]]] <- obj$mouse_to_rabbit_ratio
  scaled_obj$mouse.to.rabbit.ratio[is.infinite(scaled_obj$mouse.to.rabbit.ratio)] <- NA
  scaled_lst[[i]] <- scaled_obj
}

for (i in 1:length(scaled_lst)) {
  scaled_obj <- scaled_lst[[i]]
  # png(paste0("./dimplot_by_ident_new/dimplot_by_species_ident_sep_by_m2r_ratio_", 
  #            unique(scaled_obj$orig.ident), ".png"),
  #     width = 750, height = 750)
  # p <- DimPlot(scaled_obj, 
  #              group.by = "species.ident.by.ratio") +
  #   ggtitle("species separated by reads mapped to mouse:rabbit ratio")
  # print(p)
  # dev.off()
  # 
  # png(paste0("./dimplot_by_ident_new/dimplot_by_species_ident_sep_by_pct_mapped_", 
  #            unique(scaled_obj$orig.ident), ".png"),
  #     width = 750, height = 750)
  # p <- DimPlot(scaled_obj, 
  #              group.by = "species.ident.new") +
  #   ggtitle("species separated by % reads mapped to genome")
  # print(p)
  # dev.off()
  
  png(paste0("./dimplot_by_m2r_ratio/dimplot_by_m2r_ratio_",
             unique(scaled_obj$orig.ident), ".png"),
      width = 750, height = 750)
  p <- FeaturePlot(scaled_obj, features = "mouse.to.rabbit.ratio") +
    ggtitle("mouse:rabbit ratio by cluster")
  print(p)
  dev.off()
}

# count number of rabbit vs mouse cells ####
species_count <- matrix(nrow = 3, ncol = length(scaled_lst), 
                        dimnames = list(c("total", "rabbit", "mouse"),
                                        rep_len("", 6)))
for (i in 1:length(scaled_lst)) {
  scaled_obj <- scaled_lst[[i]]
  colnames(species_count)[i] <- as.character(unique(scaled_obj$orig.ident))
  species_count[1, i] <- length(scaled_obj$mouse.to.rabbit.ratio)
  species_count[2, i] <- sum(scaled_obj$mouse.to.rabbit.ratio > 25)
  species_count[3, i] <- sum(scaled_obj$mouse.to.rabbit.ratio < 25)
}
write.table(species_count, "mouse_rabbit_counts_new.csv", quote = F, sep = ",")

# look at total UMI per cell in rabbit vs mouse ####
rabbit <- subset(scaled_lst[[1]], species.ident.by.ratio == "rabbit")
mouse <- subset(scaled_lst[[1]], species.ident.by.ratio == "mouse")
hist(colSums(rabbit@assays$RNA@counts), breaks = 10000, xlim = c(0, 10000))
hist(colSums(mouse@assays$RNA@counts), breaks = 100)
min(colSums(rabbit@assays$RNA@counts))
max(colSums(rabbit@assays$RNA@counts))

# 04/20/2023 use pct. reads mapped to mouse/rabbit genome to assign cell identity ####
thres <- c(0.95, 0.915, 0.91, 0.92, 0.91, 0.905)
species_count <- matrix(nrow = 3, ncol = length(objlist), 
                        dimnames = list(c("total", "rabbit", "mouse"),
                                        rep_len("", 6)))
for (i in 1:length(objlist)) {
  obj <- objlist[[i]]
  obj$mouse_or_rabbit <- "rabbit"
  obj$mouse_or_rabbit[obj$pct_reads_mapped_to_mouse > thres[i]] <- "mouse"
  colnames(species_count)[i] <- as.character(unique(obj$orig.ident))
  species_count[1, i] <- length(obj$mouse_or_rabbit)
  species_count[2, i] <- sum(obj$mouse_or_rabbit == "rabbit")
  species_count[3, i] <- sum(obj$mouse_or_rabbit == "mouse")
  objlist[[i]] <- obj
}

scaled_lst <- vector(mode = "list", length(objlist))
dimreduclust <- function(obj){
  obj <- ScaleData(obj)
  obj <- RunPCA(obj,
                verbose = T,
                seed.use = 2022)
  obj <- RunUMAP(obj,
                 reduction = "pca",
                 dims = 1:30,
                 seed.use = 2022)
  obj <- FindNeighbors(obj,
                       reduction = "pca",
                       dims = 1:30)
  obj <- FindClusters(obj,
                      resolution = 0.5,
                      random.seed = 2022)
  return(obj)
  
}
for (i in 1:length(objlist)){
  print(paste0("now at ", i))
  # proceed with normalization
  obj <- PercentageFeatureSet(objlist[[i]], pattern = c("^MT-"), col.name = "percent.mt")
  obj <- FindVariableFeatures(obj, nfeatures = 3000)
  obj <- SCTransform(obj, vars.to.regress = "percent.mt", method = "glmGamPoi", 
                     return.only.var.genes = F, variable.features.n = 8000, 
                     seed.use = 2022, verbose = T)
  scaled_lst[[i]] <- dimreduclust(obj)
}

for (i in 1:length(scaled_lst)) {
  obj <- scaled_lst[[i]]
  png(paste0("./dimplot_sep_by_rabbit_mouse_20230420/dimplot_by_mouse_pct_90_",
             unique(obj$orig.ident), ".png"),
      width = 750, height = 750)
  p <- DimPlot(obj, group.by = "mouse_or_rabbit") +
    ggtitle(paste0("species separated by pct mapped to mouse > ", thres[i], " - ", 
                   as.character(unique(obj$orig.ident))))
  print(p)
  dev.off()
}
write.table(species_count, "mouse_rabbit_counts_using_mouse_pct_20230420.csv", quote = F, sep = ",")
save(scaled_lst, file = "rabbit_mouse_labeled_cells_20230420.RData")


DefaultAssay(integrated)
Idents(integrated)

integrated$condition <-
  str_split(string = integrated$group.ident,
            pattern = '-',
            simplify = T)[, 1]

integrated_rabbit <-
  integrated[, integrated$mouse_or_rabbit == "rabbit"]

DimPlot(integrated_rabbit,
        group.by = "condition")


integrated_rabbit <-
  SCTransform(integrated_rabbit, do.correct.umi = T)
integrated_rabbit <-
  RunPCA(integrated_rabbit)
integrated_rabbit <-
  RunUMAP(integrated_rabbit, 
          dims = 1:30)


integrated_rabbit
DimPlot(integrated_rabbit,
        group.by = "condition")
DimPlot(integrated_rabbit,group.by = "condition")


ncol(integrated_rabbit)
sum(integrated_rabbit$condition == "20088") # 12634
sum(integrated_rabbit$condition == "60060") # 9915

saveRDS(integrated_rabbit,
        file = "integrated_rabbit_w_conditions.RDs")
