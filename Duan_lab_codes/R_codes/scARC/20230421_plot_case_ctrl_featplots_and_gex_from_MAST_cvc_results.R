# Chuxuan Li 04/21/2023
# plot feature plot of case and control overlayed on clustering
#and gene expresion of some genes from MAST case vs control results

# init ####
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(stringr)
library(MAST)
library(GGally)
library(limma)
library(data.table)
library(NMF)
library(rsvd)
library(reshape2)
library(knitr)
library(GSEABase)

load("~/NVME/scARC_Duan_018/018-029_combined_analysis/integrated_labeled_018-029_RNAseq_obj_with_scDE_covars.RData")

# plot dimplot to check case and control status overlapping with UMAP ####
DimPlot(integrated_labeled, group.by = "aff", cols = c("red", "transparent")) +
  ggtitle("018-029 integrated RNAseq object\ndistribution of SZ case cells") +
  theme(text = element_text(size = 10)) +
  NoLegend()
DimPlot(integrated_labeled, group.by = "aff", cols = c("transparent", "blue"))+
  ggtitle("018-029 integrated RNAseq object\ndistribution of control cells") +
  theme(text = element_text(size = 10)) +
  NoLegend()

# feature plot of certain genes ####
unique(integrated_labeled$cell.type)
GABA_case_cells <- colnames(integrated_labeled)[integrated_labeled$cell.type == "GABA" &
                                   integrated_labeled$aff == "case"]
nmglut_case_cells <- colnames(integrated_labeled)[integrated_labeled$cell.type == "nmglut" &
                                     integrated_labeled$aff == "case"]
npglut_case_cells <- colnames(integrated_labeled)[integrated_labeled$cell.type == "npglut" &
                                     integrated_labeled$aff == "case"]
GABA_control_cells <- colnames(integrated_labeled)[integrated_labeled$cell.type == "GABA" &
                                                     integrated_labeled$aff == "control"]
nmglut_control_cells <- colnames(integrated_labeled)[integrated_labeled$cell.type == "nmglut" &
                                                       integrated_labeled$aff == "control"]
npglut_control_cells <- colnames(integrated_labeled)[integrated_labeled$cell.type == "npglut" &
                                                       integrated_labeled$aff == "control"]
# CSMD1, GABA
FeaturePlot(integrated_labeled, features = "CSMD1", cells = GABA_case_cells,
                  split.by = "time.ident", keep.scale = "all") + 
    theme(legend.position = "right") +
    scale_y_continuous(sec.axis = dup_axis(name = ""))

FeaturePlot(integrated_labeled, features = "CSMD1", cells = GABA_control_cells,
            split.by = "time.ident", keep.scale = "all") + 
  theme(legend.position = "right") +
  scale_y_continuous(sec.axis = dup_axis(name = ""))

# CACNA1C, nmglut
FeaturePlot(integrated_labeled, features = "CACNA1C", cells = nmglut_case_cells,
             split.by = "time.ident", keep.scale = "all") + 
  theme(legend.position = "right") +
  scale_y_continuous(sec.axis = dup_axis(name = ""))

FeaturePlot(integrated_labeled, features = "CACNA1C", cells = nmglut_control_cells,
            split.by = "time.ident", keep.scale = "all") + 
  theme(legend.position = "right") +
  scale_y_continuous(sec.axis = dup_axis(name = ""))

# OPCML, npglut
FeaturePlot(integrated_labeled, features = "OPCML", cells = npglut_case_cells,
            split.by = "time.ident", keep.scale = "all") + 
  theme(legend.position = "right") +
  scale_y_continuous(sec.axis = dup_axis(name = ""))

FeaturePlot(integrated_labeled, features = "OPCML", cells = npglut_control_cells,
            split.by = "time.ident", keep.scale = "all") + 
  theme(legend.position = "right") +
  scale_y_continuous(sec.axis = dup_axis(name = ""))

# ACAT2, npglut
FeaturePlot(integrated_labeled, features = "ACAT2", cells = npglut_case_cells,
            split.by = "time.ident", keep.scale = "all") + 
  theme(legend.position = "right") +
  scale_y_continuous(sec.axis = dup_axis(name = ""))

FeaturePlot(integrated_labeled, features = "ACAT2", cells = npglut_control_cells,
            split.by = "time.ident", keep.scale = "all") + 
  theme(legend.position = "right") +
  scale_y_continuous(sec.axis = dup_axis(name = ""))

# vlnplot of certain genes (from MAST) ####
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/MAST_scDE/scaRAW_obj_GABA_0hr.RData")
df <- as(scaRaw["CSMD1", ], "data.table")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/MAST_scDE/scaRAW_obj_GABA_1hr.RData")
df_to_append <- as(scaRaw["CSMD1", ], "data.table")
df <- rbind(df, df_to_append)
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/MAST_scDE/scaRAW_obj_GABA_6hr.RData")
df_to_append <- as(scaRaw["CSMD1", ], "data.table")
df <- rbind(df, df_to_append)
g = ggplot(df, aes(x = aff, y = value, color = aff)) + geom_jitter() + 
  facet_wrap(~time.ident) + 
  ggtitle("CSMD1 GABA") +
  geom_violin()


df <- zlmCond@coefD["CSMD1", ]
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/MAST_scDE/MAST_zlmCond_GABA_1hr.RData")
df_to_append <- as(scaRaw["CSMD1", ], "data.table")
df <- rbind(df, df_to_append)
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/MAST_scDE/scaRAW_obj_GABA_6hr.RData")
df_to_append <- as(scaRaw["CSMD1", ], "data.table")
df <- rbind(df, df_to_append)
ggplot(df, aes(x = aff, y = value, color = aff)) + geom_jitter() + 
  facet_wrap(~time.ident) + 
  ggtitle("CSMD1 GABA 0hr") +
  geom_violin()

# scatter plot 
df[, lmPred := lm(value~aff)$fitted, key = primerid]
g + aes(x=aff) + geom_line(aes(y = lmPred), lty=1)

# dotplot of certain genes (SCT normalized counts) ####
GABA <- subset(integrated_labeled, cell.type == "GABA")
GABA_case <- subset(GABA, aff == "case")
GABA_ctrl <- subset(GABA, aff == "control")

nmglut <- subset(integrated_labeled, cell.type == "nmglut")
nmglut_case <- subset(nmglut, aff == "case")
nmglut_ctrl <- subset(nmglut, aff == "control")

npglut <- subset(integrated_labeled, cell.type == "npglut")
npglut_case <- subset(npglut, aff == "case")
npglut_ctrl <- subset(npglut, aff == "control")

# CSMD1, GABA
p_case <- DotPlot(GABA_case, features = "CSMD1", group.by = "time.ident", assay = "SCT") + coord_flip()
p_ctrl <- DotPlot(GABA_ctrl, features = "CSMD1", group.by = "time.ident", assay = "SCT") + coord_flip()
df_case <- p_case$data[order(p_case$data$id), ]
df_case$aff <- "case"
df_ctrl <- p_ctrl$data[order(p_ctrl$data$id), ]
df_ctrl$aff <- "control"
df <- rbind(df_case, df_ctrl)
df$aff <- factor(df$aff, levels = c("control", "case"))
ggplot(df, aes(x = id, y = aff)) +
  geom_point(aes(color = avg.exp.scaled, size = pct.exp)) +
  theme_light() +
  scale_color_gradientn(colours = brewer.pal(5, "Reds")) +
  labs(title = "Expression of CSMD1 in GABA", x = "time point", 
       color = "average\nexpression", size = "percent cells\nexpressing gene") + 
  theme(axis.text = element_text(size = 10)) 

# CACNA1C, nmglut
p_case <- DotPlot(nmglut_case, features = "CACNA1C", group.by = "time.ident") + coord_flip()
p_ctrl <- DotPlot(nmglut_ctrl, features = "CACNA1C", group.by = "time.ident") + coord_flip()
df_case <- p_case$data[order(p_case$data$id), ]
df_case$aff <- "case"
df_ctrl <- p_ctrl$data[order(p_ctrl$data$id), ]
df_ctrl$aff <- "control"
df <- rbind(df_case, df_ctrl)
df$aff <- factor(df$aff, levels = c("control", "case"))
ggplot(df, aes(x = id, y = aff)) +
  geom_point(aes(color = avg.exp.scaled, size = pct.exp)) +
  theme_light() +
  scale_color_gradientn(colours = brewer.pal(5, "Reds")) +
  labs(title = "Expression of CACNA1C in nmglut", x = "time point", 
       color = "average\nexpression", size = "percent cells\nexpressing gene") + 
  theme(axis.text = element_text(size = 10)) 

# OPCML, npglut
p_case <- DotPlot(npglut_case, features = "OPCML", group.by = "time.ident") + coord_flip()
p_ctrl <- DotPlot(npglut_ctrl, features = "OPCML", group.by = "time.ident") + coord_flip()
df_case <- p_case$data[order(p_case$data$id), ]
df_case$aff <- "case"
df_ctrl <- p_ctrl$data[order(p_ctrl$data$id), ]
df_ctrl$aff <- "control"
df <- rbind(df_case, df_ctrl)
df$aff <- factor(df$aff, levels = c("control", "case"))
ggplot(df, aes(x = id, y = aff)) +
  geom_point(aes(color = avg.exp, size = pct.exp)) +
  theme_light() +
  scale_color_gradientn(colours = brewer.pal(5, "Reds")) +
  labs(title = "Expression of OPCML in npglut", x = "time point", 
       color = "average\nexpression", size = "percent cells\nexpressing gene") + 
  theme(axis.text = element_text(size = 10)) 

# ACAT2, npglut
p_case <- DotPlot(npglut_case, features = "ACAT2", group.by = "time.ident") + coord_flip()
p_ctrl <- DotPlot(npglut_ctrl, features = "ACAT2", group.by = "time.ident") + coord_flip()
df_case <- p_case$data[order(p_case$data$id), ]
df_case$aff <- "case"
df_ctrl <- p_ctrl$data[order(p_ctrl$data$id), ]
df_ctrl$aff <- "control"
df <- rbind(df_case, df_ctrl)
df$aff <- factor(df$aff, levels = c("control", "case"))
ggplot(df, aes(x = id, y = aff)) +
  geom_point(aes(color = avg.exp, size = pct.exp)) +
  theme_light() +
  scale_color_gradientn(colours = brewer.pal(5, "Reds")) +
  labs(title = "Expression of ACAT2 in npglut", x = "time point", 
       color = "average\nexpression", size = "percent cells\nexpressing gene") + 
  theme(axis.text = element_text(size = 10)) 
