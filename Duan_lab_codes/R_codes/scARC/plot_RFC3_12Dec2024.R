# Siwei 21 Dec 2024
# make violin plot of RFC3 in nmglut

# init ####
{
  library(Seurat)
  library(Signac)
  
  library(ggplot2)
  library(ggrepel)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(stringr)
  # library(NMF)
  # library(rsvd)
  library(RColorBrewer)
  library(MAST)

}

{
  plan("multisession", workers = 4)
  options(expressions = 20000)
  options(future.globals.maxSize = 207374182400)
  
}



integrated_100_lines_aff <-
  readRDS("integrated_018_030_21_clusters_w_aff.RData")

# downsample to 20K cells
# downsampled_20K_integrated <-
#   integrated_100_lines_aff[, 
#                            sample(colnames(integrated_100_lines_aff),
#                                   size = 20000,
#                                   replace = F), ]
# 
# Idents(downsampled_20K_integrated) <- "cell.type"
# ncol(integrated_100_lines_aff)


Idents(integrated_100_lines_aff) <- "cell.type"
DefaultAssay(integrated_100_lines_aff) <- "SCT"
DefaultAssay(integrated_100_lines_aff) <- "integrated"

VlnPlot(integrated_100_lines_aff,
        features = "RFC3",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "nmglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0,
        # y.max = 5,
        # log = T,
        same.y.lims = T) +
  ggtitle("RFC3_nmglut")


VlnPlot(integrated_100_lines_aff,
        features = "RFC3",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "npglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0,
        # y.max = 5,
        log = T,
        same.y.lims = T) +
  ggtitle("RFC3_npglut")

VlnPlot(integrated_100_lines_aff,
        features = "RFC3",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "GABA",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0,
        # y.max = 5,
        log = T,
        same.y.lims = T) +
  ggtitle("RFC3_GABA")

# NAALADL2
VlnPlot(integrated_100_lines_aff,
        features = "NAALADL2",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "nmglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0,
        # y.max = 5,
        # log = T,
        same.y.lims = T) +
  ggtitle("NAALADL2_nmglut")


VlnPlot(integrated_100_lines_aff,
        features = "NAALADL2",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "npglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0,
        # y.max = 5,
        log = T,
        same.y.lims = T) +
  ggtitle("NAALADL2_npglut")

VlnPlot(integrated_100_lines_aff,
        features = "NAALADL2",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "GABA",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0,
        # y.max = 5,
        log = T,
        same.y.lims = T) +
  ggtitle("NAALADL2_GABA")

# HMGCS1
DefaultAssay(integrated_100_lines_aff) <- "SCT"

VlnPlot(integrated_100_lines_aff,
        features = "HMGCS1",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "nmglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0,
        # y.max = 5,
        log = F,
        same.y.lims = T) +
  ggtitle("HMGCS1_nmglut")


VlnPlot(integrated_100_lines_aff,
        features = "HMGCS1",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "npglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0,
        # y.max = 5,
        log = F,
        same.y.lims = T) +
  ggtitle("HMGCS1_npglut")

VlnPlot(integrated_100_lines_aff,
        features = "HMGCS1",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "GABA",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0,
        # y.max = 5,
        log = T,
        same.y.lims = T) +
  ggtitle("HMGCS1_GABA")


# GRIK2
DefaultAssay(integrated_100_lines_aff) <- "integrated"

VlnPlot(integrated_100_lines_aff,
        features = "GRIK2",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "nmglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0.01,
        alpha = 0.01,
        # y.max = 5,
        log = F,
        # adjust = c(scale = "width"),
        same.y.lims = T) +
  ggtitle("GRIK2_nmglut")


VlnPlot(integrated_100_lines_aff,
        features = "HMGCS1",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "npglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0,
        # y.max = 5,
        log = F,
        same.y.lims = T) +
  ggtitle("HMGCS1_npglut")

VlnPlot(integrated_100_lines_aff,
        features = "HMGCS1",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "GABA",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0,
        # y.max = 5,
        log = T,
        same.y.lims = T) +
  ggtitle("HMGCS1_GABA")

# PTPRK
DefaultAssay(integrated_100_lines_aff) <- "integrated"

VlnPlot(integrated_100_lines_aff,
        features = "PTPRK",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "nmglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0.0,
        alpha = 0.01,
        # y.max = 5,
        log = F,
        # adjust = c(scale = "width"),
        same.y.lims = T) +
  ggtitle("PTPRK_nmglut")

# INSIG1
DefaultAssay(integrated_100_lines_aff) <- "SCT"
DefaultAssay(integrated_100_lines_aff) <- "integrated"
VlnPlot(integrated_100_lines_aff,
        features = "INSIG1",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "nmglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0.0, 
        alpha = 0.01,
        # y.max = 5,
        log = T,
        # adjust = c(scale = "width"),
        same.y.lims = T) +
  ggtitle("INSIG1_nmglut")

# OPCML
# DefaultAssay(integrated_100_lines_aff) <- "SCT"
DefaultAssay(integrated_100_lines_aff) <- "integrated"
VlnPlot(integrated_100_lines_aff,
        features = "OPCML",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "npglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0.0, 
        alpha = 0.01,
        # y.max = 5,
        log = F,
        # adjust = c(scale = "width"),
        same.y.lims = T) +
  ggtitle("OPCML_npglut")




summary(data.table(integrated_100_lines_aff$orig.ident))

Idents(integrated_100_lines_aff) <- "aff"
UMAPPlot(integrated_100_lines_aff,
         pt.size = 0.5,
         alpha = 0.5)

length(unique(integrated_100_lines_aff$orig.ident))

table(integrated_100_lines_aff$orig.ident)

# load the "harmony" object


load("018-030_RNA_integrated_labeled_with_harmony.RData")
table(integrated_labeled$cell.line.ident)



# CD_02, CD_03, CD_04, CD_06, CD_07, CD_09, CD_10, CD_11, CD_12, CD_13, CD_14, CD_15, CD_16, CD_17, CD_18, CD_19, CD_21, CD_22, CD_23, CD_31, CD_32, CD_33, CD_34, CD_35, CD_36, CD_37, CD_38, CD_39, CD_40, CD_42, CD_43, CD_44, CD_45, CD_46, CD_47, CD_48, CD_49, CD_50, CD_51, CD_52, CD_53, CD_55, CD_56, CD_57, CD_58, CD_60, CD_61, CD_62, CD_63, CD_64, CD_65, CD_66, CW20058, CW20063, CW20079, CW20112

sample_all <-
  "CD_02, CD_03, CD_04, CD_06, CD_07, CD_09, CD_10, CD_11, CD_12, CD_13, CD_14, CD_15, CD_16, CD_17, CD_18, CD_19, CD_21, CD_22, CD_23, CD_31, CD_32, CD_33, CD_34, CD_35, CD_36, CD_37, CD_38, CD_39, CD_40, CD_42, CD_43, CD_44, CD_45, CD_46, CD_47, CD_48, CD_49, CD_50, CD_51, CD_52, CD_53, CD_55, CD_56, CD_57, CD_58, CD_60, CD_61, CD_62, CD_63, CD_64, CD_65, CD_66, CW20058, CW20063, CW20079, CW20112"

sample_all <-
  str_split_1(string = sample_all,
              pattern = ', ')

integrated_56_matched_samples <-
  integrated_labeled[, integrated_labeled$cell.line.ident %in% sample_all]
saveRDS(integrated_56_matched_samples,
        file = "integrated_56_samples_case_ctrl_21Dec2024.RDs")

integrated_56_matched_samples <-
  readRDS("integrated_56_samples_case_ctrl_21Dec2024.RDs")

aff_ident <-
  read.csv("100line_metadata_df.csv")
case_samples <-
  aff_ident$cell_line[aff_ident$aff == "case"]

integrated_56_matched_samples$aff <-
  "ctrl"
integrated_56_matched_samples$aff[integrated_56_matched_samples$cell.line.ident %in% case_samples] <-
  "case"


Idents(integrated_56_matched_samples) <- "aff"
UMAPPlot(integrated_56_matched_samples,
         cols = c("darkred",
                  "darkblue"),
         pt.size = 0.5,
         shuffle = T,
         alpha = 0.5)
UMAPPlot(integrated_56_matched_samples,
         cols = c("darkred",
                  "darkblue"),
         split.by = "aff",
         pt.size = 0.5,
         shuffle = T,
         alpha = 0.5)

# saveRDS(integrated_56_matched_samples,
#         file = "integrated_56_samples_case_ctrl_21Dec2024.RDs")

integrated_56_matched_samples$aff <-
  factor(integrated_56_matched_samples$aff,
         levels = c("ctrl",
                    "case"))

DefaultAssay(integrated_56_matched_samples) <- "integrated"
Idents(integrated_56_matched_samples) <- "cell.type"
VlnPlot(integrated_56_matched_samples,
        features = "RFC3",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "GABA",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0.0,
        alpha = 0.0,
        # y.max = 5,
        log = F,
        same.y.lims = T) +
  geom_boxplot(width = 0.25,
               outliers = F,
               position = position_dodge(0.9)) +
  ggtitle("RFC3_GABA")

VlnPlot(integrated_56_matched_samples,
        features = "PTPRK",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "nmglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0.0,
        alpha = 0.0,
        # y.max = 5,
        log = F,
        same.y.lims = T) +
  geom_boxplot(width = 0.25,
               outliers = F,
               position = position_dodge(0.9)) +
  ggtitle("PTPRK_nmglut")

DefaultAssay(integrated_56_matched_samples) <- "SCT"
VlnPlot(integrated_56_matched_samples,
        features = "HMGCS1",
        cols = brewer.pal(n = 3,
                          name = "Dark2"),
        idents = "nmglut",
        group.by = "time.ident",
        split.by = "aff",
        pt.size = 0.0,
        alpha = 0.0,
        # y.max = 5,
        log = T,
        same.y.lims = T) +
  geom_boxplot(width = 0.25,
               outliers = F,
               position = position_dodge(0.9),
               alpha = 0.1) +
  ggtitle("HMGCS1_nmglut")
