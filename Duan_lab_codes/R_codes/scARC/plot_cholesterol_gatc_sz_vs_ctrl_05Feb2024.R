# Siwei 05 Feb 2024
# Plot cholesterol gene activity scores, SZ vs control,
# npglut, 0, 1, 6 hr

# init ####
library(Seurat)
library(Signac)

# library(harmony)

library(future)

library(Matrix)
library(tidyverse)

library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)

library(stringr)
library(ggplot2)
# library(scales)
library(RColorBrewer)
# library(plyranges)
# library(gplots)
# library(grDevices)
# library(viridis)
# library(colorspace)

plan("multisession", workers = 4)
set.seed(42)
options(future.globals.maxSize = 229496729600)
options(Seurat.object.assay.version = "v5")

# temp function

# load data ####
# load("~/NVME/scARC_Duan_018/018-029_combined_analysis/npglut_geneActivity_list_all_times.RData")
load("~/NVME/scARC_Duan_018/018-029_combined_analysis/multiomic_obj_4_plotting_gene_activity_29Jan2024.RData")
load("~/NVME/scARC_Duan_018/case_control_DE_analysis/018-029_covar_table_final.RData")

# make affxtime ident
# for (i in 1:length(npglut_geneActivity)) {
#   DefaultAssay(npglut_geneActivity[[i]]) <-
#     "gact5"
# }
# npglut_all <-
#   merge(npglut_geneActivity[[1]],
#         y = c(npglut_geneActivity[[2]],
#               npglut_geneActivity[[3]]))
DefaultAssay(subset_multiomic_obj) <- "gact5"
npglut_all <-
  subset_multiomic_obj[, subset_multiomic_obj$RNA.cell.type == "npglut"]
# make affxtime ident

covar_table_final_0hr <-
  covar_table_final[covar_table_final$time == "0hr", ]
covar_line_aff <-
  data.frame(cell.line = covar_table_final_0hr$cell_line,
             aff = covar_table_final_0hr$aff)
length(match(x = covar_line_aff$cell.line,
             table = npglut_all$cell.line.ident))
head(match(x = npglut_all$cell.line.ident,
             table = covar_line_aff$cell.line))
covar_line_aff$aff[match(x = npglut_all$cell.line.ident,
           table = covar_line_aff$cell.line)]
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
DefaultAssay(npglut_all) <- "gact5"

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
  # scale_fill_gradient(legend = "Gene Activity") +
  guides(fill = guide_legend(title = "Gene Activity")) +
  labs(x = "Gene",
       y = "") +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Chl genes, gene activity")

p_dotplot <-
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
  # scale_fill_gradient(legend = "Gene Activity") +
  guides(fill = guide_legend(title = "Gene Activity")) +
  labs(x = "Gene",
       y = "") +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("Chl genes, npglut\nGene activity")
# p_dotplot$labels <- "Gene Activity"
p_dotplot


## plot BDNF in npglut+nmglut, MvsF, 016 hrs
glut_all <-
  subset_multiomic_obj[, subset_multiomic_obj$RNA.cell.type %in% c("nmglut",
                                                                   "npglut")]

covar_line_sex <-
  data.frame(cell.line = covar_table_final_0hr$cell_line,
             sex = covar_table_final_0hr$sex)

glut_all$sex <-
  covar_line_sex$sex[match(x = glut_all$cell.line.ident,
                           table = covar_line_sex$cell.line)]
glut_all$sex_time.ident <-
  str_c(glut_all$sex,
        glut_all$time.ident,
        sep = "_")

glut_all$aff <-
  covar_line_aff$aff[match(x = glut_all$cell.line.ident,
                           table = covar_line_aff$cell.line)]


glut_all$cell.type_sex_time.ident <-
  str_c(glut_all$RNA.cell.type,
        glut_all$sex,
        glut_all$time.ident,
        sep = "_")

unique(glut_all$cell.type_sex_time.ident)
glut_all$cell.type_sex_time.ident <-
  factor(glut_all$cell.type_sex_time.ident,
         levels = c("nmglut_F_0hr",
                    "nmglut_M_0hr",
                    "nmglut_F_1hr",
                    "nmglut_M_1hr",
                    "nmglut_F_6hr",
                    "nmglut_M_6hr",
                    "npglut_F_0hr",
                    "npglut_M_0hr",
                    "npglut_F_1hr",
                    "npglut_M_1hr",
                    "npglut_F_6hr",
                    "npglut_M_6hr"))
Idents(glut_all) <- "cell.type_sex_time.ident"
DefaultAssay(glut_all) <- "gact5"

DotPlot(glut_all,
        features = c("BDNF"),
        cols = c("lightgrey",
                 brewer.pal(n = 4, name = "Set1")[1])) +
  # scale_fill_gradient(legend = "Gene Activity") +
  guides(fill = guide_legend(title = "Gene Activity")) +
  labs(x = "Gene",
       y = "") +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("stim_response genes, gene activity")

p_data <-
  DotPlot(glut_all,
          features = c("BDNF"),
          cols = c("lightgrey",
                   brewer.pal(n = 4, name = "Set1")[1])) +
  # scale_fill_gradient(legend = "Gene Activity") +
  guides(fill = guide_legend(title = "Gene Activity")) +
  labs(x = "Gene",
       y = "") +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("stim_response genes, gene activity")

df_2_plot <-
  p_data$data
df_2_plot <-
  df_2_plot[, c(2, 4, 5)]

df_2_plot$cell.type <-
  str_split(df_2_plot$id,
            pattern = "_",
            simplify = T)[, 1]
df_2_plot$sex <-
  str_split(df_2_plot$id,
            pattern = "_",
            simplify = T)[, 2]
df_2_plot$time <-
  str_split(df_2_plot$id,
            pattern = "_",
            simplify = T)[, 3]
df_2_plot$cell.typexsex <-
  str_c(df_2_plot$cell.type,
        df_2_plot$sex,
        sep = "_")

df_2_plot <-
  reshape2::melt(data = df_2_plot,
                 id = c("pct.exp",
                        "avg.exp.scaled"))
df_2_plot <-
  reshape2::melt(data = df_2_plot,
                 id = c("id",
                        "cell_type",
                        "sex",
                        "time"))

ggplot(df_2_plot,
       aes(x = cell.typexsex,
           y = time,
           size = pct.exp,
           fill = avg.exp.scaled)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low = "lightgrey",
                      high = brewer.pal(n = 4, name = "Set1")[1]) +
  labs(x = "",
       y = "Time") +
  # guides(fill = guide_legend(title = "avg.gact.scaled")) +
  scale_y_discrete(limits = rev) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  ggtitle("BDNF, gene activity")




glut_all$timextype.ident <-
  factor(glut_all$timextype.ident,
         levels = c("nmglut_0hr",
                    "nmglut_1hr",
                    "nmglut_6hr",
                    "npglut_0hr",
                    "npglut_1hr",
                    "npglut_6hr"))
Idents(glut_all) <- "timextype.ident"
DotPlot(glut_all,
        features = c("FOS", 
                     "JUNB",
                     "NR4A1",
                     "NR4A3",
                     "VGF",
                     "BDNF",
                     "PCSK1",
                     "DUSP4"),
        cols = c("lightgrey",
                 brewer.pal(n = 4, name = "Set1")[1])) +
  # scale_fill_gradient(legend = "Gene Activity") +
  guides(fill = guide_legend(title = "Gene Activity")) +
  labs(x = "Gene",
       y = "") +
  scale_x_discrete(limits = rev) +
  RotatedAxis() +
  coord_flip() +
  ggtitle("stim_response genes, gene activity")
