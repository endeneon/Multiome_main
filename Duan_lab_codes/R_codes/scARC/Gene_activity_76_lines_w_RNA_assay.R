# Siwei 22 Jan 2024
# Use 18 lines, add GeneActivity assay

# init ####
library(Seurat)
library(Signac)

# library(harmony)

library(future)

library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)

library(stringr)
library(ggplot2)
library(scales)
library(RColorBrewer)
# library(plyranges)
library(gplots)
library(grDevices)
library(viridis)
library(colorspace)

plan("multisession", workers = 4)
set.seed(42)
options(future.globals.maxSize = 229496729600)

##### load data
# load("~/NVME/scARC_Duan_018/018-029_combined_analysis/integrated_labeled_018-029_RNAseq_obj_with_scDE_covars.RData")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/df_gene_activity_counts_76_lines_26Jan2024.RData")
# use this object #
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/018-029_multiomic_obj_clustered_added_metadata.RData")
# load("~/NVME/scARC_Duan_018/018-029_combined_analysis/018-029_RNA_integrated_labeled_with_harmony.RData")
# load("~/NVME/scARC_Duan_018/018-029_combined_analysis/multiomic_obj_with_new_peaks_labeled.RData")

# DefaultAssay(integrated_labeled) <- "integrated"
# DimPlot(integrated_labeled)
# ncol(counts) # 477107
# ncol(integrated_labeled) # 533333
# 
# sum(colnames(integrated_labeled) %in% 
#       colnames(counts))
# head(colnames(integrated_labeled))
# # > head(colnames(integrated_labeled))
# # [1] "AAACAGCCAAGGTATA-1_1" "AAACAGCCAGGATGGC-1_1" "AAACAGCCAGTAAAGC-1_1" "AAACAGCCATCCAGGT-1_1"
# # [5] "AAACAGCCATGGCCTG-1_1" "AAACATGCAACCTGGT-1_1"
# head(colnames(counts))
# # > head(colnames(counts))
# # [1] "AAACAGCCAAACAACA-62" "AAACAGCCAAACCTAT-59" "AAACAGCCAAACCTTG-24" "AAACAGCCAAACTAAG-40"
# # [5] "AAACAGCCAAACTCAT-15" "AAACAGCCAAACTGTT-31"
# unique(str_split(string = colnames(integrated_labeled),
#                  pattern = "_",
#                  simplify = T)[, 2])
# # integrated_labeled_backup <-
# #   integrated_labeled
# 
# new_colnames_integrated_labeled <-
#   str_c(str_split(string = colnames(integrated_labeled),
#                   pattern = "1_",
#                   simplify = T)[, 1],
#         str_split(string = colnames(integrated_labeled),
#                   pattern = "1_",
#                   simplify = T)[, 2])
# length(unique(new_colnames_integrated_labeled))
# sum(new_colnames_integrated_labeled %in% 
#       colnames(counts))
# unique(integrated_labeled$orig.ident)
# 
# ## trim excessive barcodes
# counts_subsetted <-
#   counts[, colnames(counts) %in% colnames(multiomic_obj_new)]
# # Now need to figure out
# 
# 
# ## construct gact
# integrated_labeled[['gact']] <-
#   CreateAssay5Object()

DefaultAssay(multiomic_obj) <- "ATAC"
DimPlot(object = multiomic_obj, 
        label = F, 
        cols = brewer.pal(n = 3,
                          name = "Set1"),
        group.by = "time.ident",
        pt.size = 1,
        alpha = 0.5) + 
  ggtitle("by time") +
  theme(text = element_text(size = 10), 
        axis.text = element_text(size = 10))

DefaultAssay(integrated_labeled_backup)
DimPlot(object = integrated_labeled_backup, 
        label = F, 
        cols = brewer.pal(n = 3,
                          name = "Set1"),
        group.by = "time.ident",
        pt.size = 1,
        alpha = 0.5) + 
  ggtitle("by time") +
  theme(text = element_text(size = 10), 
        axis.text = element_text(size = 10))

ncol(multiomic_obj)
ncol(counts)

colnames(multiomic_obj)
colnames(counts)

counts_subsetted <-
  counts[, match(x = colnames(multiomic_obj),
                 table = colnames(counts),
                 nomatch = 0)]
ncol(counts_subsetted)

rownames(counts_subsetted)
Signac::Annotation(multiomic_obj)

multiomic_obj_backup <-
  multiomic_obj

multiomic_obj[['gact']] <-
  CreateAssayObject(counts = counts_subsetted)
DefaultAssay(multiomic_obj) <- "gact"
multiomic_obj <-
  NormalizeData(object = multiomic_obj,
                assay = "gact",
                normalization.method = 'LogNormalize',
                scale.factor = median(multiomic_obj$nCount_gact))
rownames(multiomic_obj)
save(multiomic_obj,
     file = "multiomic_obj_w_gene_activity_normalised_29Jan2024.RData")

multiomic_obj_backup <-
  multiomic_obj

FeaturePlot(object = multiomic_obj,
            features = "GAD1",
            max.cutoff = 'q95')

Features(multiomic_obj) <-
  gene_activity_names

# rownames(multiomic_obj) <- gene_activity_names
# multiomic_obj <-
#   UpdateSeuratObject(multiomic_obj)
# Features(multiomic_obj)
# rownames(multiomic_obj) <- gene_activity_names
# Features(multiomic_obj)

multiomic_obj[['gact5']] <-
  as(object = multiomic_obj[['gact']],
     Class = "Assay5")
DefaultAssay(multiomic_obj) <- "gact5"
rownames(multiomic_obj) <- gene_activity_names
Features(multiomic_obj)

FeaturePlot(object = multiomic_obj,
            features = "GAD1",
            # cols = c("lightgrey", "darkblue"),
            max.cutoff = "q90")
FeaturePlot(object = multiomic_obj,
            features = "GAD2",
            # cols = c("lightgrey", "darkblue"),
            max.cutoff = "q90")

FeaturePlot(object = multiomic_obj,
            features = "SLC17A6",
            # cols = c("lightgrey", "darkblue"),
            max.cutoff = "q90")
FeaturePlot(object = multiomic_obj,
            features = "SLC17A7",
            # cols = c("lightgrey", "darkblue"),
            max.cutoff = "q95")

FeaturePlot(object = multiomic_obj,
            features = "FOS",
            # cols = c("lightgrey", "darkblue"),
            max.cutoff = "q95")
FeaturePlot(object = multiomic_obj,
            features = "NEFM",
            # cols = c("lightgrey", "darkblue"),
            max.cutoff = "q90")


DefaultAssay(multiomic_obj) <- "gact5"
unique(multiomic_obj$RNA.cell.type)
unique(multiomic_obj$time.ident)
rownames(multiomic_obj) <- gene_activity_names
multiomic_obj$time.ident <-
  factor(multiomic_obj$time.ident,
         levels = c("0hr", "1hr", "6hr"))

# "npglut"       "GABA"         "unidentified" "nmglut" 
## make dot plots
DefaultAssay(multiomic_obj) <- "gact5"
Idents(multiomic_obj) <- "RNA.cell.type"
sum(multiomic_obj$RNA.cell.type %in% "npglut")
Features(multiomic_obj)
unique(multiomic_obj$RNA.cell.type)
Idents(multiomic_obj) <- "RNA.cell.type"

subset_multiomic_obj <-
  subset(multiomic_obj,
         cells = (multiomic_obj$RNA.cell.type %in% c("npglut", "nmglut", "GABA")))
rownames(multiomic_obj)
subset_multiomic_obj <-
  multiomic_obj[, multiomic_obj$RNA.cell.type %in% c("npglut", "nmglut", "GABA")]

DefaultAssay(subset_multiomic_obj) <- "gact5"
Idents(subset_multiomic_obj) <- "RNA.cell.type"

# multiomic_obj
subset_multiomic_obj$time.ident <-
  factor(subset_multiomic_obj$time.ident,
         levels = c("0hr", "1hr", "6hr"))

cell_time <-
  c("0hr", "1hr", "6hr")
cell_type <-
  c('NEFM+ glut',
    'GABA',
    'NEFM- glut')
unique(subset_multiomic_obj$timextype.ident)

Idents(subset_multiomic_obj) <- 
  factor(subset_multiomic_obj$timextype.ident,
         levels = sort(unique(subset_multiomic_obj$timextype.ident)))


DotPlot(subset_multiomic_obj,
        features = c("FOS", "JUNB", "NR4A1", "NR4A3",
                     "BTG2", "ATF3", "DUSP1", "EGR1"),
        cols = c("lightgrey",
                 brewer.pal(n = 4, name = "Set1")[1])) +
  guides(fill = "Gene Activity") +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Early response, gene activity")


DotPlot(subset_multiomic_obj,
        features = c("VGF", "BDNF", "PCSK1", "DUSP4",
                     "ATP1B1", "SLC7A5", "NPTX1", "SCG2", 
                     "CREM"),
        cols = c("lightgrey",
                 brewer.pal(n = 4, name = "Set1")[1])) +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Late response, gene activity")

DotPlot(subset_multiomic_obj,
        features = c("SLC17A6", "SLC17A7",
                     "GAD1", "GAD2"),
        cols = c("lightgrey",
                 brewer.pal(n = 4, name = "Set1")[1])) +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Cell Type-specific Genes, gene activity")

save(subset_multiomic_obj,
     file = "multiomic_obj_4_plotting_gene_activity_29Jan2024.RData")


## plot heatmap
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
  ggtitle("Early response, gene activity")

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
  ggtitle("Late response, gene activity")

df_collect <-
  as.data.frame(rbind(df_collect,
                      plot_values$data))


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

dev.off()
pdf(file = "Gene_act_heatmap_GABA.pdf",
    width = 2,
    height = 5.5,
    onefile = T)
heatmap.2(as.matrix(df_2_heatmap2),
          Rowv = F, Colv = F,
          dendrogram = "none",
          scale = "row",
          col = colorRampPalette(colors = viridis(n = 12))(100),
          trace = "none",
          density.info = "none",
          cexCol = 1,
          margins = c(5, 5))
dev.off()

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

dev.off()
pdf(file = "Gene_act_heatmap_nmglut.pdf",
    width = 2,
    height = 5.5,
    onefile = T)
heatmap.2(as.matrix(df_2_heatmap2),
          Rowv = F, Colv = F,
          dendrogram = "none",
          scale = "row",
          col = colorRampPalette(colors = viridis(n = 12))(100),
          trace = "none",
          density.info = "none",
          cexCol = 1,
          margins = c(5, 5))
dev.off()


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
pdf(file = "Gene_act_heatmap_npglut.pdf",
    width = 2,
    height = 5.5,
    onefile = T)
heatmap.2(as.matrix(df_2_heatmap2),
          Rowv = F, Colv = F,
          dendrogram = "none",
          scale = "row",
          col = colorRampPalette(colors = viridis(n = 12))(100),
          trace = "none",
          density.info = "none",
          cexCol = 1,
          margins = c(5, 5))
dev.off()


### subset multiomics_obj for Alexi's GeneActivity analysis
DefaultAssay(multiomic_obj) <- "gact5"

npglut_geneActivity <-
  vector(mode = "list",
         length = 3L)
names(npglut_geneActivity) <-
  c('hr_0', 'hr_1', 'hr_6')

reduced_obj <-
  multiomic_obj
reduced_obj@assays$ATAC <- NULL
reduced_obj@assays$gact <- NULL

npglut_geneActivity[['hr_0']] <-
  reduced_obj[, reduced_obj$timextype.ident %in% "npglut_0hr"]
npglut_geneActivity[['hr_1']] <-
  reduced_obj[, reduced_obj$timextype.ident %in% "npglut_1hr"]
npglut_geneActivity[['hr_6']] <-
  reduced_obj[, reduced_obj$timextype.ident %in% "npglut_6hr"]
unique(reduced_obj$timextype.ident)

save(npglut_geneActivity,
     file = "npglut_geneActivity_list_all_times.RData")
