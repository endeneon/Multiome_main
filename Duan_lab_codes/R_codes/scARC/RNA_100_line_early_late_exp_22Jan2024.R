# Siwei 22 Jan 2024
# Use Lexi's 100 line RNA-seq data 
# Confirm that 100 line RNA-seq GEX is consistent with 18 line GEX



# init ####
library(Seurat)
library(Signac)
# library(sctransform)
# library(glmGamPoi)
# library(EnsDb.Hsapiens.v86)
# library(GenomeInfoDb)
# library(GenomicFeatures)
# library(AnnotationDbi)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

# library(patchwork)
# library(readr)
library(ggplot2)
library(RColorBrewer)
library(gplots)

# library(dplyr)
# library(viridis)
# library(graphics)
library(stringr)
library(future)

library(grDevices)
library(viridis)
library(colorspace)

plan("multisession", workers = 4)
set.seed(42)
options(future.globals.maxSize = 229496729600)

## load data #####
load("~/NVME/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/018-030_RNA_integrated_labeled_with_harmony.RData")

## calc cellular composition
DefaultAssay(integrated_labeled) <- "SCT"

# cellular composition ####
cell_lines <- 
  sort(unique(integrated_labeled$cell.line.ident))
time_type_sum <- 
  vector(mode = "list", length = 3L)
numlines <- 
  length(cell_lines)
cell_types <- 
  as.vector(unique(integrated_labeled$cell.type.forplot))
cell_types <-
  cell_types[c(1, 2, 4)]
numtypes <- 
  length(cell_types)
cell_times <- 
  unique(integrated_labeled$time.ident)

for (k in 1:length(cell_times)) {
  print(cell_times[k])
  df <- data.frame(cell.type = rep_len(NA, numlines * numtypes),
                   cell.line = rep_len(NA, numlines * numtypes),
                   counts = rep_len(0, numlines*numtypes))
  for (i in 1:length(cell_lines)) {
    subobj <- subset(integrated_labeled, 
                     subset = (cell.line.ident == cell_lines[i] & 
                       time.ident == cell_times[k]))
    for (j in 1:numtypes) {
      print(numtypes * i - numtypes + j)
      df$cell.type[numtypes * i - numtypes + j] <- cell_types[j]
      df$cell.line[numtypes * i - numtypes + j] <- cell_lines[i]
      df$counts[numtypes * i - numtypes + j] <- 
        sum(subobj$cell.type.forplot == cell_lines[j])
    }
  }
  time_type_sum[[k]] <- df
}

# plot early-responsive genes
Idents(integrated_labeled) <- "cell.type.forplot"

# DefaultAssay(integrated_labeled) <- "integrated"
# ! Use "SCT" slot !
DefaultAssay(integrated_labeled) <- "SCT"


DotPlot(integrated_labeled,
        features = c("FOS", "JUNB", "NR4A1", "NR4A3",
                     "BTG2", "ATF3", "DUSP1", "EGR1"),
        cols = c(brewer.pal(n = 4, name = "Set1")[1], 
                 brewer.pal(n = 4, name = "Set1")[1], 
                 brewer.pal(n = 4, name = "Set1")[1]),
        idents = cell_types,
        # group.by = "cell.type.forplot",
        split.by = "time.ident") +
  scale_y_discrete(limits = rev) +
  RotatedAxis() +
  ggtitle("Early response, 100 lines")

DotPlot(integrated_labeled,
        features = c("VGF", "BDNF", "PCSK1", "DUSP4",
                     "ATP1B1", "SLC7A5", "NPTX1", "SCG2", 
                     "CREM"),
        cols = c(brewer.pal(n = 4, name = "Set1")[1], 
                 brewer.pal(n = 4, name = "Set1")[1], 
                 brewer.pal(n = 4, name = "Set1")[1]),
        idents = cell_types,
        # group.by = "cell.type.forplot",
        split.by = "time.ident") +
  scale_y_discrete(limits = rev) +
  RotatedAxis() +
  ggtitle("Late response, 100 lines")

DotPlot(integrated_labeled,
        features = c("SLC17A6", "SLC17A7",
                     "GAD1", "GAD2"),
        cols = c(brewer.pal(n = 4, name = "Set1")[1], 
                 brewer.pal(n = 4, name = "Set1")[1], 
                 brewer.pal(n = 4, name = "Set1")[1]),
        idents = cell_types,
        # group.by = "cell.type.forplot",
        split.by = "time.ident") +
  scale_y_discrete(limits = rev) +
  RotatedAxis() +
  ggtitle("Cell Type-specific Genes, 100 lines")

## Get Duan_024 groups
Duan_024_libs <-
  unique(integrated_labeled$lib.ident)
Duan_024_libs <-
  Duan_024_libs[str_detect(string = Duan_024_libs,
                           pattern = "^5-|^33-|^35-|^51-|^63-")]

Duan_024_barcodes <-
  vector(mode = "list",
         length = length(Duan_024_libs))
names(Duan_024_barcodes) <-
  Duan_024_libs

for (i in 1:length(Duan_024_barcodes)) {
  Duan_024_barcodes[[i]] <-
    vector(mode = "list",
           length = length(cell_types))
  # names(Duan_024_barcodes)[[i]] <-
  #   cell_types
}

for (i in 1:length(Duan_024_barcodes)) {
  # Duan_024_barcodes[[i]] <-
  #   vector(mode = "list",
  #          length = length(cell_types))
  names(Duan_024_barcodes[[i]]) <-
    cell_types
}
# names(Duan_024_barcodes[1])

for (i in 1:length(Duan_024_barcodes)) {
  for (j in 1:length(cell_types)) {
    print(paste(i, j,
                sep = ","))
     temp_barcode <-
      colnames(integrated_labeled)[(integrated_labeled$lib.ident == names(Duan_024_barcodes[i])) &
                                     (integrated_labeled$cell.type.forplot == names(Duan_024_barcodes[[i]])[j])]
     temp_barcode <-
       str_split(string = temp_barcode,
                 pattern = "_",
                 simplify = T)[, 1]
    Duan_024_barcodes[[i]][[j]] <- temp_barcode
  }
}

cell_saving_types <-
  c("npglut",
    "GABA",
    "nmglut")

if (!dir.exists("Duan_024_barcodes")) {
  dir.create("Duan_024_barcodes")
}

for (i in 1:length(Duan_024_barcodes)) {
  for (j in 1:length(cell_types)) {
    print(paste(i, j,
                sep = ","))
    
    table_2_write <-
      Duan_024_barcodes[[i]][[j]] 
    write.table(table_2_write,
                file = paste0("Duan_024_barcodes",
                              "/",
                              "Duan_024_",
                              names(Duan_024_barcodes[i]),
                              "_",
                              names(Duan_024_barcodes[[i]])[j],
                              ".txt"),
                row.names = F, col.names = F,
                sep = "\t", quote = F)
  }
}

# Remake DotPlot, use 3 cell types only #####
# 30 Jan 2024
unique(integrated_labeled$cell.type)
unique(integrated_labeled$time.ident)

## use the same name as used in Gene Activity Assay to facilitate code ####
subset_multiomic_obj <-
  integrated_labeled[, integrated_labeled$cell.type %in% c("npglut",
                                                           "nmglut",
                                                           "GABA")]

subset_multiomic_obj$time.ident <-
  factor(subset_multiomic_obj$time.ident,
         levels = c("0hr", "1hr", "6hr"))

subset_multiomic_obj$timextype.ident <- NA

# cell_type <-
#   c('NEFM+ glut',
#     'GABA',
#     'NEFM- glut')
cell_type <-
  c("GABA",
    "nmglut",
    "npglut")
cell_time <-
  c("0hr", "1hr", "6hr")

for (i in 1:length(cell_type)) {
  for (j in 1:length(cell_time)) {
    print(paste("i =", i,
                "j = ", j))
    subset_multiomic_obj$timextype.ident[(subset_multiomic_obj$cell.type == cell_type[i]) &
                                           (subset_multiomic_obj$time.ident == cell_time[j])] <-
      str_c(cell_type[i],
            "_",
            cell_time[j])
  }
}
unique(subset_multiomic_obj$timextype.ident)

Idents(subset_multiomic_obj) <- "timextype.ident"
Idents(subset_multiomic_obj) <- 
  factor(subset_multiomic_obj$timextype.ident,
         levels = sort(unique(subset_multiomic_obj$timextype.ident)))

DotPlot(subset_multiomic_obj,
        features = c("FOS", "JUNB", "NR4A1", "NR4A3",
                     "BTG2", "ATF3", "DUSP1", "EGR1"),
        cols = c("lightgrey",
                 brewer.pal(n = 4, name = "Set1")[1])) +
  guides(fill = "Gene Expression") +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Early response, gene expression")


DotPlot(subset_multiomic_obj,
        features = c("VGF", "BDNF", "PCSK1", "DUSP4",
                     "ATP1B1", "SLC7A5", "NPTX1", "SCG2", 
                     "CREM"),
        cols = c("lightgrey",
                 brewer.pal(n = 4, name = "Set1")[1])) +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Late response, gene expression")

DotPlot(subset_multiomic_obj,
        features = c("SLC17A6", "SLC17A7",
                     "GAD1", "GAD2"),
        cols = c("lightgrey",
                 brewer.pal(n = 4, name = "Set1")[1])) +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Cell type-specific genes, gene expression")

# save(subset_multiomic_obj,
#      file = "multiomic_obj_4_plotting_gene_expression_30Jan2024.RData")

plot_values <-
  DotPlot(subset_multiomic_obj,
          features = c("SLC17A6", "SLC17A7",
                       "GAD1", "GAD2"),
          cols = c("lightgrey",
                   brewer.pal(n = 4, name = "Set1")[1])) +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Cell type-specific genes, gene expression")
plot_values$data

plot_values <-
  DotPlot(subset_multiomic_obj,
          features = c("FOS", "JUNB", "NR4A1", "NR4A3",
                       "BTG2", "ATF3", "DUSP1", "EGR1"),
          cols = c("lightgrey",
                   brewer.pal(n = 4, name = "Set1")[1])) +
  guides(fill = "Gene Expression") +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Early response, gene expression")

df_collect <-
  plot_values$data

plot_values <-
  DotPlot(subset_multiomic_obj,
          features = c("VGF", "BDNF", "PCSK1", "DUSP4",
                       "ATP1B1", "SLC7A5", "NPTX1", "SCG2", 
                       "CREM"),
          cols = c("lightgrey",
                   brewer.pal(n = 4, name = "Set1")[1])) +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Late response, gene expression")

df_collect <-
  as.data.frame(rbind(df_collect,
                      plot_values$data))

# 
# df_2_plot <-
#   plot_values$data

ggplot(df_collect,
       aes(x = id,
           y = features.plot,
           fill = avg.exp)) +
  geom_bin2d() +
  scale_fill_gradientn(colours = diverge_hcl(n = 12),
                       guide = "colourbar") +
  # scale_fill_gradient2(low = "#80FFFF",
  #                      mid = "#FFFFFF",
  #                      high = "#FF80FF",
  #                      midpoint = 0,
  #                      guide = "colourbar") +
  theme_classic() +
  labs(x = "Cell type x Time",
       y = "Gene") +
  # guides(fill = guide_legend("Expression value")) +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))

df_2_heatmap2 <-
  as.data.frame(cbind(df_collect$avg.exp[df_collect$id == "GABA_0hr"],
                      df_collect$avg.exp[df_collect$id == "GABA_1hr"],
                      df_collect$avg.exp[df_collect$id == "GABA_6hr"]))
rownames(df_2_heatmap2) <-
  df_collect$features.plot[df_collect$id == "GABA_0hr"]
colnames(df_2_heatmap2) <-
  c("GABA_0hr",
    "GABA_1hr",
    "GABA_6hr")


df_2_heatmap2 <-
  as.data.frame(cbind(df_collect$avg.exp[df_collect$id == "nmglut_0hr"],
                      df_collect$avg.exp[df_collect$id == "nmglut_1hr"],
                      df_collect$avg.exp[df_collect$id == "nmglut_6hr"]))
rownames(df_2_heatmap2) <-
  df_collect$features.plot[df_collect$id == "nmglut_0hr"]
colnames(df_2_heatmap2) <-
  c("nmglut_0hr",
    "nmglut_1hr",
    "nmglut_6hr")

df_2_heatmap2 <-
  as.data.frame(cbind(df_collect$avg.exp[df_collect$id == "npglut_0hr"],
                      df_collect$avg.exp[df_collect$id == "npglut_1hr"],
                      df_collect$avg.exp[df_collect$id == "npglut_6hr"]))
rownames(df_2_heatmap2) <-
  df_collect$features.plot[df_collect$id == "npglut_0hr"]
colnames(df_2_heatmap2) <-
  c("npglut_0hr",
    "npglut_1hr",
    "npglut_6hr")

dev.off()
heatmap.2(as.matrix(df_2_heatmap2),
          Rowv = F, Colv = F,
          dendrogram = "none",
          scale = "row",
          col = colorRampPalette(colors = colorspace::diverge_hcl(n = 12))(100),
          trace = "none",
          density.info = "none",
          cexCol = 1,
          margins = c(5, 5))


## plot chl genes exp
load("~/NVME/scARC_Duan_018/case_control_DE_analysis/018-029_covar_table_final.RData")
DefaultAssay(subset_multiomic_obj) <- "SCT"

meta_100_lines <-
  readr::read_csv("100line_metadata_df.csv")
meta_100_lines_0hr <-
  meta_100_lines[meta_100_lines$time_point == "0hr", ]
covar_line_aff <-
  data.frame(cell.line = meta_100_lines_0hr$cell_line,
             aff = meta_100_lines_0hr$aff)

npglut_all <-
  subset_multiomic_obj[, subset_multiomic_obj$cell.type == "npglut"]
unique(npglut_all$cell.line.ident)
npglut_all$aff <-
  covar_line_aff$aff[match(x = npglut_all$cell.line.ident,
                           table = covar_line_aff$cell.line)]
npglut_all$aff_time.ident <-
  str_c(npglut_all$aff,
        npglut_all$time.ident,
        sep = "_")

npglut_all$aff_time.ident <-
  factor(npglut_all$aff_time.ident,
         levels = c("control_0hr",
                    "case_0hr",
                    "control_1hr",
                    "case_1hr",
                    "control_6hr",
                    "case_6hr"))
Idents(npglut_all) <- "aff_time.ident"

# unique(subset_multiomic_obj$cell.line.ident)

DotPlot(npglut_all,
        features = c("ACAT2",
                     "HMGCR",
                     "SQLE",
                     "MSMO1",
                     "FDFT1",
                     "HMGCS1",
                     "INSIG1",
                     # "CYP51A1",
                     "DHCR24",
                     "SCD",
                     # "NSDHL",
                     "SC5D"),
        cols = c("lightgrey",
                 brewer.pal(n = 4, name = "Set1")[1])) +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Chl genes, npglut\ngene expression")

p_results <-
  DotPlot(npglut_all,
          features = c("ACAT2",
                       "HMGCR",
                       "SQLE",
                       "MSMO1",
                       "FDFT1",
                       "HMGCS1",
                       "INSIG1",
                       # "CYP51A1",
                       "DHCR24",
                       "SCD",
                       # "NSDHL",
                       "SC5D"),
          cols = c("lightgrey",
                   brewer.pal(n = 4, name = "Set1")[1])) +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Chl genes, npglut\ngene expression")

View(p_results$data)
write.table(p_results$data,
            file = "chl_gex_output_table.txt",
            quote = F, sep = "\t",
            row.names = F, col.names = T)
