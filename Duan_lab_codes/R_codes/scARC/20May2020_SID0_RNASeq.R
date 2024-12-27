# 20 May 2020 Siwei
# Run single-timepoint RNAseq use SID0

# init
library(Seurat)
library(future)
library(readr)
library(sctransform)
library(stringr)
library(ggplot2)


###
library(lazyeval)
library(plyr)
library(glue)
library(car)
library(dplyr)

# setup multithread
plan("multisession", workers = 8) # should not use "multicore" here
plan("multisession", workers = 1) 


# read in SID0 data
SID0.data <- Read10X("../scRNA/SID0_scRNA/outs/filtered_feature_bc_matrix/")
SID0_scRNA_095_ASE_cell_line_index <- read_table2("~/Galaxy_Drive/May15_Duan_017/demultiplex/SID0_scRNA_095_ASE/SID0_scRNA_095_ASE_cell_line_index.txt")
# remove cells without clear identities first before proceed to clustering 
# to remove inteference from doublets
SID0.data <- SID0.data[, colnames(SID0.data) %in% 
                         SID0_scRNA_095_ASE_cell_line_index$BARCODE]

# SID0_cell_line_index <- read_table2("~/Galaxy_Drive/May15_Duan_017/demultiplex/SID0_scRNA/SID0_cell_line_index.txt")

###### use new sctransform wrapper#######
SID0.Seruat <- CreateSeuratObject(counts = SID0.data, project = "SID0.RNASeq")
# store mitochondrial percentage in object meta data
# rownames(SID0.data)[str_detect(rownames(SID0.data), regex("^MT-", TRUE))]
SID0.Seruat <- PercentageFeatureSet(SID0.Seruat, pattern = "^MT-|^Mt-", # need to include rat mt
                                    col.name = "percent.mt") 
SID0.Seruat <- SCTransform(SID0.Seruat, 
                           vars.to.regress = "percent.mt", 
                           verbose = T)
## make PCA and UMAP data and store
SID0.Seruat <- RunPCA(SID0.Seruat, verbose = T)
SID0.Seruat <- RunUMAP(SID0.Seruat, dims = 1:30, verbose = T)

SID0.Seruat <- FindNeighbors(SID0.Seruat, 
                             # k.param = 15, # set k here
                             dims = 1:30, 
                             verbose = T)
SID0.Seruat <- FindClusters(SID0.Seruat, 
                            resolution = 0.5, # adjust cluster number here, 0.5 ~ 14 clusters
                            verbose = T)

DimPlot(SID0.Seruat, label = T) +
  NoLegend()

########### find some marker genes ########
all.markers.SID0 <- FindAllMarkers(SID0.Seruat, 
                                   logfc.threshold = 0.5, 
                                   min.pct = 0.25,
                                   verbose = T)



########### check some marker genes #######
## at 14 clusters (resolution = 0.5)
VlnPlot(SID0.Seruat, features = c("GAD2", "GAD1", "HOXA5", # 0+7, GABA
                                  "SLC17A7", "GRIN1", "GLS", # 1+5+10+12, glut
                                  "VIM", "NES", "SOX2", # 2+11+13+3+6, NPC and immature GABA
                                  "S100b", "Slc1a3", "Gfap", # 4+8, rat astrocyte
                                  "MT-CO3", "MT-CO1", "MT-CYB"), # 9 dead cells
        pt.size = 0, ncol = 4)

VlnPlot(SID0.Seruat, features = c("SATB2", "CBLN2", "RASGRF2", "RORB", "FOXP2", "TSHZ2", "HTR2C", 
                                  "PCP4", "NEFM", "NR4A2", "GAD1", 
                                  "GAD2", "CALB2", "CNR1", "STXBP6", "SST", "MOBP", "S100B", "SLC1A3", "GFAP",
                                  "NES", "SOX2", "VIM"), # 9 dead cells
        pt.size = 0, ncol = 5)

VlnPlot(SID0.Seruat, features = c("PNOC"), # 9 dead cells
        pt.size = 0, ncol = 1)

FeaturePlot(SID0.Seruat, features = c("SST", "PVALB"), 
            ncol = 2)

FeaturePlot(SID0.Seruat, features = c("SLC17A7", "SLC17A6", 
                                      "BCL11B", "SOX2"), 
            ncol = 2)

FeaturePlot(SID0.Seruat, features = c("CUX1", "POU3F2",
                                      "SATB2", "BCL11B", 
                                      "TBR1", "DLG4", 
                                      "SLC32A1", "DLX2"), 
            ncol = 3)

VlnPlot(SID0.Seruat, features = c("CUX1", "POU3F2",
                                  "SATB2", "BCL11B", 
                                  "TBR1", "DLG4", 
                                  "SLC32A1", "DLX2"), # 9 dead cells
        pt.size = 0, ncol = 3)


top100_markers <- FindAllMarkers(SID0.Seruat, only.pos = T, 
                                 logfc.threshold = 0.25) %>% 
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_logFC)
write.table(top100_markers, file = "top100_markers_14_clusters.txt",
            quote = F, sep = "\t", row.names = T, col.names = T)



SID0_umap$UMAP_1 <- FetchData(SID0.Seruat, vars = "UMAP_1")
SID0_umap$UMAP_2 <- FetchData(SID0.Seruat, vars = "UMAP_2")
SID0_umap$clusters <- FetchData(SID0.Seruat, vars = "ident")
SID0_umap$barcodes <- rownames(SID0.Seruat@reductions$umap)


SID0_cell_individuals <- base::merge(x = SID0_umap, 
                                     y = SID0_scRNA_095_ASE_cell_line_index, 
                                     by.x = "barcodes", 
                                     by.y = "BARCODE")
SID0_cell_individuals$UMAP_1 <- SID0_cell_individuals$UMAP_1$UMAP_1
SID0_cell_individuals$UMAP_2 <- SID0_cell_individuals$UMAP_2$UMAP_2
SID0_cell_individuals$clusters <- SID0_cell_individuals$clusters$ident
SID0_cell_individuals[SID0_cell_individuals$BEST %in% "SNG-iPS_line_15", 9] <- "SNG-CN_line_15"
SID0_cell_individuals[SID0_cell_individuals$clusters %in% c(4,8), 9] <- "Rat Astrocytes"


ggplot(SID0_cell_individuals, aes(x = UMAP_1, y = UMAP_2, 
                                  colour = BEST)) +
  geom_point(size = 0.5) +
  scale_colour_manual(name = "Cell type", 
                      values = c("black", "red", "green4", "darkblue")) +
  ggtitle("Unstimulated") +
  theme_classic()

return_vlnplot <- VlnPlot(SID0.Seruat, features = c("GAD2", "GAD1", "HOXA5", # 0+7, GABA
                                                    "SLC17A7", "GRIN1", "GLS", # 1+5+10+12, glut
                                                    "VIM", "NES", "SOX2", # 2+11+13+3+6, NPC and immature GABA
                                                    "S100b", "Slc1a3", "Gfap", # 4+8, rat astrocyte
                                                    "MT-CO3", "MT-CO1", "MT-CYB"), # 9 dead cells
                          pt.size = 0, ncol = 4)



###################### make large violin plot #####
violin_plot_scaled <- SID0.Seruat@assays$SCT@scale.data[
  rownames(SID0.Seruat@assays$SCT@scale.data) %in% 
    c("SLC17A1", "SATB2", "CBLN2", "RASGRF2", "RORB", "FOXP2", "TSHZ2", "HTR2C", 
      "PCP4", "NEFM", "NR4A2", "GAD1", 
      "GAD2", "CALB2", "CNR1", "STXBP6", "SST", "MOBP", "S100B", "SLC1A3", "GFAP",
      "NES", "SOX2", "VIM"), ]

violin_plot_scaled <- FetchData(SID0.Seruat, vars = "SCT")

SID0.Seruat <- BuildClusterTree(SID0.Seruat)
PlotClusterTree(SID0.Seruat)

####### test customised Vlnplot_custom ########
Vlnplot_custom(SID0.Seruat, features = c("GAD2", "GAD1", "HOXA5", # 0+7, GABA
                                  "SLC17A7", "GRIN1", "GLS", # 1+5+10+12, glut
                                  "VIM", "NES", "SOX2", # 2+11+13+3+6, NPC and immature GABA
                                  "S100b", "Slc1a3", "Gfap", # 4+8, rat astrocyte
                                  "MT-CO3", "MT-CO1", "MT-CYB"), # 9 dead cells
        pt.size = 0, ncol = 1)







######### customise Vlnplot_custom ########

Vlnplot_custom <- function (object, features, cols = NULL, pt.size = 1, idents = NULL, 
                            sort = FALSE, assay = NULL, group.by = NULL, split.by = NULL, 
                            adjust = 1, y.max = NULL, same.y.lims = FALSE, log = FALSE, 
                            ncol = NULL, slot = "data", split.plot = FALSE, combine = TRUE) 
{
  if (!is.null(x = split.by) & getOption(x = "Seurat.warn.vlnplot.split", 
                                         default = TRUE)) {
    message("The default behaviour of split.by has changed.\n", 
            "Separate violin plots are now plotted side-by-side.\n", 
            "To restore the old behaviour of a single split violin,\n", 
            "set split.plot = TRUE.\n      \nThis message will be shown once per session.")
    options(Seurat.warn.vlnplot.split = FALSE)
  }
  return(ExIPlot_custom(object = object, type = ifelse(test = split.plot, 
                                                yes = "splitViolin", no = "violin"), features = features, 
                 idents = idents, ncol = ncol, sort = sort, assay = assay, 
                 y.max = y.max, same.y.lims = same.y.lims, adjust = adjust, 
                 pt.size = pt.size, cols = cols, group.by = group.by, 
                 split.by = split.by, log = log, slot = slot, combine = combine))
}

##### customised ExIPlot_custom #####
ExIPlot_custom <- function (object, features, type = "violin", idents = NULL, 
                            ncol = NULL, sort = FALSE, assay = NULL, y.max = NULL, same.y.lims = FALSE, 
                            adjust = 1, cols = NULL, pt.size = 0, group.by = NULL, split.by = NULL, 
                            log = FALSE, slot = "data", combine = TRUE) 
{
  # assay <- assay || DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  # ncol <- ncol || ifelse(test = length(x = features) > 9, 
  #                          yes = 4, no = min(length(x = features), 3))
  ncol <- min(length(x = features))
  data <- FetchData(object = object, vars = features, slot = slot)
  features <- colnames(x = data)
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  }
  else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% 
                                                 idents])
  }
  data <- data[cells, , drop = FALSE]
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  }
  else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  if (is.null(x = split.by)) {
    split <- NULL
  }
  else {
    split <- object[[split.by, drop = TRUE]][cells]
    if (!is.factor(x = split)) {
      split <- factor(x = split)
    }
    if (is.null(x = cols)) {
      cols <- hue_pal()(length(x = levels(x = idents)))
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    else if (length(x = cols) == 1 && cols == "interaction") {
      split <- interaction(idents, split)
      cols <- hue_pal()(length(x = levels(x = idents)))
    }
    else {
      cols <- Col2Hex(cols)
    }
    if (length(x = cols) < length(x = levels(x = split))) {
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
    names(x = cols) <- sort(x = levels(x = split))
    if ((length(x = cols) > 2) & (type == "splitViolin")) {
      warning("Split violin is only supported for <3 groups, using multi-violin.")
      type <- "violin"
    }
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  plots <- lapply(X = features, FUN = function(x) {
    return(SingleExIPlot_custom(type = type, data = data[, x, drop = FALSE], 
                                idents = idents, split = split, sort = sort, y.max = y.max, 
                                adjust = adjust, cols = cols, pt.size = pt.size, 
                                log = log))
  })
  label.fxn <- switch(EXPR = type, violin = ylab, splitViolin = ylab, 
                      ridge = xlab, stop("Unknown ExIPlot type ", type, call. = FALSE))
  for (i in 1:length(x = plots)) {
    key <- paste0(unlist(x = strsplit(x = features[i], split = "_"))[1], 
                  "_")
    obj <- names(x = which(x = Key(object = object) == key))
    if (length(x = obj) == 1) {
      if (inherits(x = object[[obj]], what = "DimReduc")) {
        plots[[i]] <- plots[[i]] + label.fxn(label = "Embeddings Value")
      }
      else if (inherits(x = object[[obj]], what = "Assay")) {
        next
      }
      else {
        warning("Unknown object type ", class(x = object), 
                immediate. = TRUE, call. = FALSE)
        plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
      }
    }
    else if (!features[i] %in% rownames(x = object)) {
      plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
    }
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = ncol)
    if (length(x = features) > 1) {
      plots <- plots & NoLegend()
    }
  }
  return(plots)
}

######## customised SingleExIPlot_custom #######
SingleExiPlot_custom <- function (data, idents, split = NULL, type = "violin", sort = FALSE, 
                                  y.max = NULL, adjust = 1, pt.size = 0, cols = NULL, seed.use = 42, 
                                  log = FALSE) 
{
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!is.data.frame(x = data) || ncol(x = data) != 1) {
    stop("'SingleExIPlot requires a data frame with 1 column")
  }
  feature <- colnames(x = data)
  data$ident <- idents
  if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
    data$ident <- factor(x = data$ident, levels = names(x = rev(x = sort(x = tapply(X = data[, 
                                                                                             feature], INDEX = data$ident, FUN = mean), decreasing = grepl(pattern = paste0("^", 
                                                                                                                                                                            tolower(x = sort)), x = "decreasing")))))
  }
  if (log) {
    noise <- rnorm(n = length(x = data[, feature]))/200
    data[, feature] <- data[, feature] + 1
  }
  else {
    noise <- rnorm(n = length(x = data[, feature]))/1e+05
  }
  if (all(data[, feature] == data[, feature][1])) {
    warning(paste0("All cells have the same value of ", 
                   feature, "."))
  }
  else {
    data[, feature] <- data[, feature] + noise
  }
  axis.label <- "Expression Level"
  y.max <- y.max %||% max(data[, feature][is.finite(x = data[, 
                                                             feature])])
  if (type == "violin" && !is.null(x = split)) {
    data$split <- split
    vln.geom <- geom_violin
    fill <- "split"
  }
  else if (type == "splitViolin" && !is.null(x = split)) {
    data$split <- split
    vln.geom <- geom_split_violin
    fill <- "split"
    type <- "violin"
  }
  else {
    vln.geom <- geom_violin
    fill <- "ident"
  }
  switch(EXPR = type, violin = {
    x <- "ident"
    y <- paste0("`", feature, "`")
    # xlab <- "Identity"
    # ylab <- axis.label
    geom <- list(vln.geom(scale = "width", adjust = adjust, 
                          trim = TRUE), theme(axis.text.x = element_blank()))
    if (is.null(x = split)) {
      jitter <- geom_jitter(height = 0, size = pt.size)
    } else {
      jitter <- geom_jitter(position = position_jitterdodge(jitter.width = 0.4, 
                                                            dodge.width = 0.9), size = pt.size)
    }
    log.scale <- scale_y_log10()
    axis.scale <- ylim
  }, ridge = {
    x <- paste0("`", feature, "`")
    y <- "ident"
    xlab <- axis.label
    ylab <- "Identity"
    geom <- list(geom_density_ridges(scale = 4), theme_ridges(), 
                 scale_y_discrete(expand = c(0.01, 0)), scale_x_continuous(expand = c(0, 
                                                                                      0)))
    jitter <- geom_jitter(width = 0, size = pt.size)
    log.scale <- scale_x_log10()
    axis.scale <- function(...) {
      invisible(x = NULL)
    }
  }, stop("Unknown plot type: ", type))
  plot <- ggplot(data = data, mapping = aes_string(x = x, 
                                                   y = y, fill = fill)[c(2, 3, 1)]) + 
    labs(x = xlab, y = ylab, 
    # title = feature, 
    fill = NULL) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
  plot <- do.call(what = "+", args = list(plot, geom))
  plot <- plot + if (log) {
    log.scale
  }
  else {
    axis.scale(min(data[, feature]), y.max)
  }
  if (pt.size > 0) {
    plot <- plot + jitter
  }
  if (!is.null(x = cols)) {
    if (!is.null(x = split)) {
      idents <- unique(x = as.vector(x = data$ident))
      splits <- unique(x = as.vector(x = data$split))
      labels <- if (length(x = splits) == 2) {
        splits
      }
      else {
        unlist(x = lapply(X = idents, FUN = function(pattern, 
                                                     x) {
          x.mod <- gsub(pattern = paste0(pattern, "."), 
                        replacement = paste0(pattern, ": "), x = x, 
                        fixed = TRUE)
          x.keep <- grep(pattern = ": ", x = x.mod, 
                         fixed = TRUE)
          x.return <- x.mod[x.keep]
          names(x = x.return) <- x[x.keep]
          return(x.return)
        }, x = unique(x = as.vector(x = data$split))))
      }
      if (is.null(x = names(x = labels))) {
        names(x = labels) <- labels
      }
    }
    else {
      labels <- levels(x = droplevels(data$ident))
    }
    plot <- plot + scale_fill_manual(values = cols, labels = labels)
  }
  return(plot)
}



########### old setup ###########

sum(colnames(SID0.data) %in% SID0_scRNA_dp01_cell_line_index$BARCODE)
sum(colnames(SID0.data) %in% SID0_scRNA_095_ASE_imputed_cell_line_index$BARCODE)
sum(colnames(SID0.data) %in% SID0_scRNA_095_ASE_cell_line_index$BARCODE)
sum(colnames(SID0.data) %in% SID0_scRNA_099_ASE_cell_line_index$BARCODE)

table(SID0_scRNA_095_ASE_cell_line_index$BEST[
       SID0_scRNA_095_ASE_cell_line_index$BARCODE %in% 
         colnames(SID0.data)])
table(SID0.Seruat@active.ident)

# ## normalise data
# SID0.normalised <- NormalizeData(SID0.data, 
#                                  normalization.method = "LogNormalize",   
#                                  scale.factor = 10000)
# ## apply a linear transformation to scale the data
# SID0.normalised <- ScaleData(SID0.normalised, 
#                              features = rownames(SID0.normalised))
# filter out cells without clear identity
# SID0.data.identity <- SID0.data[, colnames(SID0.data) %in% SID0_cell_line_index$BARCODE]
