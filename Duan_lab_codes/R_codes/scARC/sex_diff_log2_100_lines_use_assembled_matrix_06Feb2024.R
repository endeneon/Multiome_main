# Siwei 06 Feb 2024
# Test if there is any log2 diff changes caused by sex factor

# Need to use sva/Combat to get normalised counts
# Adapt to 100 lines,
# 100 lines = 76 lines from Lexi's 76 line assay + remaining lines from 
# the 100 line assay

# init ####
library(Seurat)
library(Signac)

library(limma)
library(edgeR)
library(ggplot2)
library(stringr)

library(sva)
# library(harmony)

library(reshape2)

library(readr)
library(readxl)

plan("multisession", workers = 4)
set.seed(42)
options(future.globals.maxSize = 229496729600)
options(Seurat.object.assay.version = "v5")

# Aux functions ####
createDGE <- 
  function(count_matrix) {
    y <- 
      DGEList(counts = count_matrix, 
              genes = rownames(count_matrix), 
              group = cell_lines)
    A <- rowSums(y$counts)
    hasant <- rowSums(is.na(y$genes)) == 0
    y <- y[hasant, , keep.lib.size = F]
    print(dim(y))
    return(y)
  }

filterByCpm <- 
  function(df1, df2, cutoff, proportion = 1) {
    cpm1 <- cpm(df1)
    cpm2 <- cpm(df2)
    passfilter <- (rowSums(cpm1 >= cutoff) >= ncol(cpm1) * proportion |
                     # rowSums(cpm2 >= cutoff) >= ncol(cpm2) * proportion |
                     (rowSums(cpm1 > 0) & (rowSums(cpm2 >= cutoff) >= ncol(cpm2) * proportion))) # number of samples from either group > ns
    return(passfilter)
  }

cnfV <- 
  function(y, design) {
    y <- calcNormFactors(y)
    v <- voom(y, design, plot = F)
    return(v)  
  }
# deal three cell types one by one, 
# extract them from adj_mat_lst

# load data ####
load("npglut_030_matrix.RData")
load("npglut_count_matrix_76_lines.RData")

covar_table_full <- 
  read_csv("100line_metadata_df.csv")
covar_table_final <-
  covar_table_full
covar_table_final <-
  covar_table_final[order(covar_table_final$cell_line), ]

ncol(npglut_030)
ncol(npglut_count_76)

colnames(npglut_030)
colnames(npglut_count_76)

npglut_only_030_samples <-
  npglut_030[, !(colnames(npglut_030) %in% colnames(npglut_count_76))]
colnames(npglut_only_030_samples)

# calculate a shared gene list
nrow(npglut_030) # 36601
nrow(npglut_count_76) # 34283
length(rownames(npglut_030) %in% rownames(npglut_count_76))

shared_gene_list <-
  rownames(npglut_030)[rownames(npglut_030) %in% rownames(npglut_count_76)]
shared_gene_list <-
  shared_gene_list[shared_gene_list %in% rownames(npglut_030)]

# cbind the two halves of count matrices
new_npglut_100_count_pt1 <-
  npglut_count_76[rownames(npglut_count_76) %in% shared_gene_list, ]
new_npglut_100_count_pt1 <-
  npglut_count_76[order(rownames(npglut_count_76)), ]

new_npglut_100_count_pt2 <-
  npglut_only_030_samples[rownames(npglut_only_030_samples) %in% shared_gene_list, ]
nrow(new_npglut_100_count_pt2)
new_npglut_100_count_pt2 <-
  new_npglut_100_count_pt2[order(rownames(new_npglut_100_count_pt2)), ]
rownames(new_npglut_100_count_pt2)

merged_new_npglut_100 <-
  as.data.frame(cbind(new_npglut_100_count_pt1,
                      new_npglut_100_count_pt2))
rownames(merged_new_npglut_100)


####
nrow(merged_new_npglut_100[rowSums(merged_new_npglut_100) > 0, ])

{
  master_df_raw <-
    merged_new_npglut_100
  master_df_raw <-
    master_df_raw[, order(colnames(master_df_raw))]
  master_df_raw <-
    master_df_raw[rowSums(master_df_raw) > 0, ]
  master_df_raw <-
    master_df_raw + 1
  master_list_raw <-
    vector(mode = "list",
           length = 3L)
  
  cell_time <-
    c("0hr",
      "1hr",
      "6hr")
  
  for (i in 1:3L) {
    master_list_raw[[i]] <-
      as.data.frame(master_df_raw[, str_detect(string = colnames(master_df_raw),
                                               pattern = cell_time[i])])
  }
  names(master_list_raw) <-
    cell_time
  
  ## make design matrices ####
  covar_table_final <- 
    covar_table_final[order(covar_table_final$cell_line), ]
  covar_table_final$age_sq <- 
    covar_table_final$age * covar_table_final$age
  
  df_design <-
    vector(mode = "list",
           length = length(master_list_raw))
  names(df_design) <- 
    c("covar_table_0hr", 
      "covar_table_1hr", 
      "covar_table_6hr")
  for (i in 1:3L) {
    df_design[[i]] <-
      covar_table_final[covar_table_final$time_point %in% cell_time[i], ]
  }
  
}

{
  cell_lines <- 
    unique(str_remove(colnames(master_df_raw), 
                      pattern = "_[0|1|6]hr$"))
  colnames(master_df_raw)
  
  
  passfilter_1vs0 <-
    filterByCpm(df1 = master_list_raw[[1]],
                df2 = master_list_raw[[2]],
                cutoff = 0.25,
                proportion = 1)
  passfilter_6vs0 <-
    filterByCpm(df1 = master_list_raw[[1]],
                df2 = master_list_raw[[3]],
                cutoff = 0.25,
                proportion = 1)
  passfilter_all <-
    passfilter_1vs0 & passfilter_6vs0
  sum(passfilter_all)
  passfilter_all[names(passfilter_all) == "FOSB"]
  
  matrix_design <-
    vector(mode = "list",
           length = length(master_list_raw))
  for (i in 1:length(master_list_raw)) {
    matrix_design[[i]] <-
      model.matrix(~ 0 + 
                     sex,
                   data = df_design[[i]])
  }
  
  voom_list_post_sva <-
    vector(mode = "list",
           length = length(master_list_raw))
  names(voom_list_post_sva) <-
    cell_time
  
  for (i in 1:length(master_list_raw)) {
    print(i)
    matrix_2_voom <-
      ComBat_seq(counts = as.matrix(master_list_raw[[i]][passfilter_all, ]),
                 batch = factor(df_design[[i]]$seq_batch))
    # print(head(matrix_2_voom))
    print("==1==")
    matrix_2_voom <-
      createDGE(matrix_2_voom)
    print("==2==")
    voom_list_post_sva[[i]] <-
      cnfV(y = matrix_2_voom,
           design = model.matrix(~ 0 + 
                                   GABA_porportion + # ! here change
                                   aff +
                                   sex +
                                   age_sq,
                                 data = df_design[[i]]))
    
  }
}

{
  
  log2FC_list <-
    vector(mode = "list",
           length = 2L)
  names(log2FC_list) <-
    c("log2FC_1vs0",
      "log2FC_6vs0")
  
  # genes_to_plot <-
  #   c("FOS", "JUNB", "NR4A1", "NR4A3",
  #     "BTG2", "ATF3", "DUSP1", "EGR1",
  #     "VGF", "BDNF", "PCSK1", "DUSP4",
  #     "ATP1B1", "SLC7A5", "NPTX1", "SCG2", 
  #     "CREM","SLC17A6", "SLC17A7",
  #     "GAD1", "GAD2")
  # genes_to_plot <-
  #   c("BDNF", "JUNB")
  # genes_to_plot <-
  #   c("BDNF", "JUNB",
  #     "FOSB", "FOSL2",
  #     "HOMER1", "NR4A2",
  #     "NPAS4", "MEF2C")
  genes_to_plot <-
    c("BDNF")
  
  interim_matrix <-
    voom_list_post_sva[['1hr']]$E - voom_list_post_sva[['0hr']]$E
  interim_matrix <-
    interim_matrix[rownames(interim_matrix) %in% genes_to_plot, ]
  print(nrow(interim_matrix))
  
  cell_line_type_specific <-
    colnames(interim_matrix)
  cell_line_type_specific <-
    str_remove_all(string = cell_line_type_specific,
                   pattern = "_1hr")
  
  # df_writeout <-
  #   as.data.frame(t(interim_matrix))
  # df_writeout$cell_line <-
  #   cell_line_type_specific
  # df_writeout$sex <-
  #   df_design$covar_table_0hr$sex
  # write.table(df_writeout,
  #             file = "npglut_1vs0hr_BDNF_JUNB_log2CPM.txt",
  #             quote = F, sep = "\t",
  #             row.names = F, col.names = T)
  
  
  df_2_plot <-
    melt(interim_matrix)
  
  
  
  df_2_plot$cell_line <-
    rep(cell_line_type_specific,
        each = 1)
  df_2_plot$sex <-
    rep(df_design$covar_table_1hr$sex,
        each = 1)
}

ggplot(df_2_plot,
       aes(x = factor(sex),
           y = value,
           fill = factor(sex))) +
  geom_violin(scale = "width",
              trim = F) +
  geom_boxplot(aes(ymin = min(value, na.rm = T),
                   lower = quantile(value, 0.25, na.rm = T),
                   middle = median(value, na.rm = T),
                   upper = quantile(value, 0.75, na.rm = T),
                   ymax = max(value, na.rm = T)),
               # stat = "identity",
               width = 0.3, 
               # notch = T,
               na.rm = T) +
  geom_dotplot(binaxis = "y",
               dotsize = 0.2,
               stackdir = "centerwhole") +
  geom_point(position = "dodge",
             stat = "summary",
             fun = "mean",
             shape = 4,
             colour = "white") +
  facet_wrap(~ Var1, 
             nrow = 2) +
  labs(x = "sex",
       y = "log2FC") +
  theme_classic() +
  ggtitle(paste0("npglut 1vs0 hr",
                 ", 76 lines + batch 030"))

ggplot(df_2_plot,
       aes(x = factor(sex),
           y = value,
           fill = factor(sex))) +
  geom_violin(scale = "width",
              trim = F) +
  geom_boxplot(aes(ymin = min(value, na.rm = T),
                   lower = quantile(value, 0.25, na.rm = T),
                   middle = median(value, na.rm = T),
                   upper = quantile(value, 0.75, na.rm = T),
                   ymax = max(value, na.rm = T)),
               # stat = "identity",
               width = 0.3, 
               notch = T,
               na.rm = T) +
  geom_dotplot(binaxis = "y",
               dotsize = 0.2,
               stackdir = "centerwhole") +
  geom_point(position = "dodge",
             stat = "summary",
             fun = "mean",
             shape = 4,
             colour = "white") +
  # facet_wrap(~ Var1, 
  #            nrow = 2) +
  labs(x = "sex",
       y = "log2FC") +
  guides(fill = "none") +
  theme_classic() +
  theme(axis.text = element_text(size = 12)) +
  ggtitle(paste0("BDNF, npglut 1vs0 hr",
                 ";\nP=0.012; F=44, M=56"))
  # ggtitle(paste0("npglut 1vs0 hr",
  #                ", P = 0.012"))

for (i in 1:length(genes_to_plot)) {
  try({
    print(paste(genes_to_plot[i],
                t.test(x = df_2_plot$value[(df_2_plot$Var1 == genes_to_plot[i]) &
                                             df_2_plot$sex == "M"],
                       y = df_2_plot$value[(df_2_plot$Var1 == genes_to_plot[i]) &
                                             df_2_plot$sex == "F"],
                       alternative = "t")$p.value))
  })
}


interim_matrix <-
  voom_list_post_sva[['6hr']]$E - voom_list_post_sva[['0hr']]$E
interim_matrix <-
  interim_matrix[rownames(interim_matrix) %in% genes_to_plot, ]

# df_writeout <-
#   as.data.frame(t(interim_matrix))
# df_writeout$cell_line <-
#   cell_line_type_specific
# df_writeout$sex <-
#   df_design$covar_table_0hr$sex
# df_writeout <-
#   as.data.frame(t(df_writeout))
# colnames(df_writeout) <-
#   cell_line_type_specific
# # str_split(string = colnames(df_writeout),
# #           pattern = "_",
# #           simplify = T)[, 1]
# write.table(df_writeout,
#             file = "npglut_6vs0hr_all_log2FC.txt",
#             quote = F, sep = "\t",
#             row.names = F, col.names = T)


df_2_plot <-
  melt(interim_matrix)
df_2_plot$cell_line <-
  rep(cell_lines,
      each = 8)
df_2_plot$sex <-
  rep(df_design$covar_table_0hr$sex,
      each = 8)
df_2_plot$sex <-
  rep(df_design$covar_table_0hr$sex,
      each = 1)

ggplot(df_2_plot,
       aes(x = factor(sex),
           y = value,
           fill = factor(sex))) +
  geom_violin(scale = "width",
              trim = F) +
  geom_boxplot(aes(ymin = min(value, na.rm = T),
                   lower = quantile(value, 0.25, na.rm = T),
                   middle = median(value, na.rm = T),
                   upper = quantile(value, 0.75, na.rm = T),
                   ymax = max(value, na.rm = T)), 
               width = 0.3,
               # width = position_dodge(width = 0.5),
               # stat = "identity",
               na.rm = T) +
  geom_dotplot(binaxis = "y",
               dotsize = 0.2,
               stackdir = "centerwhole") +
  geom_point(position = "dodge",
             stat = "summary",
             fun = "mean",
             shape = 4,
             colour = "white") +
  # position_dodge(width = 1) +
  facet_wrap(~ Var1, 
             nrow = 2) +
  labs(x = "sex",
       y = "log2FC") +
  theme_classic() +
  ggtitle(paste0("npglut 6vs0 hr",
                 ", 76 lines + batch 030"))

ggplot(df_2_plot,
       aes(x = factor(sex),
           y = value,
           fill = factor(sex))) +
  geom_violin(scale = "width",
              trim = F) +
  geom_boxplot(aes(ymin = min(value, na.rm = T),
                   lower = quantile(value, 0.25, na.rm = T),
                   middle = median(value, na.rm = T),
                   upper = quantile(value, 0.75, na.rm = T),
                   ymax = max(value, na.rm = T)), 
               width = 0.3,
               notch = T,
               # width = position_dodge(width = 0.5),
               # stat = "identity",
               na.rm = T) +
  geom_dotplot(binaxis = "y",
               dotsize = 0.2,
               stackdir = "centerwhole") +
  geom_point(position = "dodge",
             stat = "summary",
             fun = "mean",
             shape = 4,
             colour = "white") +
  # position_dodge(width = 1) +
  # facet_wrap(~ Var1, 
  #            nrow = 2) +
  guides(fill = "none") +
  labs(x = "sex",
       y = "log2FC") +
  theme_classic() +
  theme(axis.text = element_text(size = 12)) +
  ggtitle(paste0("BDNF, npglut 6vs0 hr",
                 ";\nP=0.040; F=44, M=56"))
sum(df_2_plot$sex == "M")

for (i in 1:length(genes_to_plot)) {
  try({
    print(paste(genes_to_plot[i],
                t.test(x = df_2_plot$value[(df_2_plot$Var1 == genes_to_plot[i]) &
                                             df_2_plot$sex == "M"],
                       y = df_2_plot$value[(df_2_plot$Var1 == genes_to_plot[i]) &
                                             df_2_plot$sex == "F"],
                       alternative = "t")$p.value))
  })
}

