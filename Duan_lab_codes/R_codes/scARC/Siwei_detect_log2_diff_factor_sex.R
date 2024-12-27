# Siwei 30 Jan 2024
# Test if there is any log2 diff changes caused by sex factor

# Need to use sva/Combat to get normalised counts
# Adapt to 100 lines

# init ####
library(Seurat)
library(Signac)

library(limma)
library(edgeR)
library(ggplot2)
library(stringr)

library(sva)
library(harmony)

library(reshape2)

library(readr)
library(readxl)

plan("multisession", workers = 4)
set.seed(42)
options(future.globals.maxSize = 229496729600)
options(Seurat.object.assay.version = "v5")
# load data
# load("./018-029_covar_table_final.RData")
# load("/nvmefs/scARC_Duan_018/case_control_DE_analysis/018-029_by_type_linextime_combat_adj_mat_lst.RData")
load("/nvmefs/scARC_Duan_018/Integrated_all_RNAseq_data_analysis/multiomic_obj_4_plotting_gene_expression_30Jan2024.RData")
# Step 1: make the count matrices list as in Lexi's code #####
## Use "RNA" assay to get raw counts (will use ComBat-Seq later)
DefaultAssay(subset_multiomic_obj) <- "RNA"

unique(subset_multiomic_obj$timextype.ident)
subset_multiomic_obj$timextype.ident <-
  factor(subset_multiomic_obj$timextype.ident,
         levels = sort(unique(subset_multiomic_obj$timextype.ident)))

unique(subset_multiomic_obj$cell.type)
subset_multiomic_obj$cell.type <-
  factor(subset_multiomic_obj$cell.type,
         levels = sort(unique(subset_multiomic_obj$cell.type)))

unique(subset_multiomic_obj$time.ident)
subset_multiomic_obj$time.ident <-
  factor(subset_multiomic_obj$time.ident,
         levels = sort(unique(subset_multiomic_obj$time.ident)))

unique(subset_multiomic_obj$cell.line.ident)
subset_multiomic_obj$cell.line.ident <-
  factor(subset_multiomic_obj$cell.line.ident,
         levels = sort(unique(subset_multiomic_obj$cell.line.ident)))
# make a function to return count matrices
# df_ident_2_batch <-
#   data.frame(ident)
# names(unique(subset_multiomic_obj$timextype.ident))
# subset_multiomic_obj$seq.batch.ident[unique(subset_multiomic_obj$timextype.ident)]

makeMat4CombatseqLine <-
  function(typeobj) {
    cell_line_count <- 
      sort(unique(typeobj$cell.line.ident))
    cell_time_count <-
      sort(unique(typeobj$time.ident))
    for (a in 1:length(cell_line_count)) {
      print(paste0("cell_time_count=",
                   a, 
                   ", ",
                   cell_line_count[a]))
      obj <-
        typeobj[,
                 typeobj$cell.line.ident == cell_line_count[a]]
      # obj <- subset(typeobj, 
      #               cell.line.ident == cell_line_count[a])
      if (a == 1) {
        for (b in 1:length(cell_time_count)) {
          print(paste0("cell_time_count=", 
                       b,
                       ", ",
                       cell_time_count[b]))
          single_time_obj <-
            obj[,
                obj$time.ident == cell_time_count[b]]
            # subset(obj,
            #        time.ident = cell_time_count[b])
          if (b == 1) {
            time_mat <- as.array(rowSums(single_time_obj@assays$RNA@counts))
            # names(time_mat) <
          } else {
            to_bind <- as.array(rowSums(single_time_obj@assays$RNA@counts))
            time_mat <- cbind(time_mat, 
                              to_bind)
          }
        }
        colnames(time_mat) <-
          str_c(cell_line_count[a],
                cell_time_count,
                sep = "_")
        mat <- time_mat
      } else {
        for (b in 1:length(cell_time_count)) {
          print(paste0("cell_time_count=", 
                       b,
                       ", ",
                       cell_time_count[b]))
          single_time_obj <-
            obj[,
                obj$time.ident == cell_time_count[b]]
          # single_time_obj <-
          #   subset(obj,
          #          time.ident = cell_time_count[b])
          if (b == 1) {
            time_mat <- as.array(rowSums(single_time_obj@assays$RNA@counts))
          } else {
            to_bind <- as.array(rowSums(single_time_obj@assays$RNA@counts))
            time_mat <- cbind(time_mat, 
                              to_bind)
          }
        }
        colnames(time_mat) <-
          str_c(cell_line_count[a],
                cell_time_count,
                sep = "_")
        mat <-
          cbind(mat,
                time_mat)
      }
    }
    return(mat)
  }
# colnames(mat) <- linetimes

#}

#####
type_obj_lsts <- 
  vector(mode = "list",
         length = 3L)
  
names(type_obj_lsts) <- 
  # c("GABA_list", "nmglut_list", npglut_list)
  sort(unique(subset_multiomic_obj$cell.type))


mat_lst <- vector("list", 3)
names(mat_lst) <-
  sort(unique(subset_multiomic_obj$cell.type))

for (i in 1:length(type_obj_lsts)) {
  print(names(type_obj_lsts)[i])
  # objlst <- type_obj_lsts[[i]]
  # for (k in 1:length(objlst)) {
  obj_input <- 
    subset_multiomic_obj[,
                         subset_multiomic_obj$cell.type == names(type_obj_lsts)[i]]
    # subset(subset_multiomic_obj,
    #             cell.type = names(type_obj_lsts)[i])
  DefaultAssay(obj_input) <- "RNA"
  mat <- makeMat4CombatseqLine(obj_input)
  # if (i == 1) {
  #   mat <- makeMat4CombatseqLine(obj_input)
  # } else {
  #   print("k > 1")
  #   mat <- cbind(mat, 
  #                makeMat4CombatseqLine(obj_input))
  # }
  # }
  mat <- mat[, order(colnames(mat))]
  mat_lst[[i]] <- mat
}
save(mat_lst,
     file = "pseudo_bulk_100_lines_data_typexlinextime_mat_lst_05Feb2024.RData")
# names(mat_lst) <- names(type_obj_lsts)

#####
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



# names(adj_mat_list)
adj_mat_lst <-
  mat_lst
names(adj_mat_list)
# [1] "GABA"   "nmglut" "npglut"
covar_table_final <-
  meta_100_lines

cell_cat <- 3L
##  ####
{
  master_df_raw <-
    adj_mat_lst[[cell_cat]]
  master_df_raw <-
    master_df_raw + 2
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
    unique(str_remove(colnames(adj_mat_lst[[cell_cat]]), 
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
  genes_to_plot <-
    c("BDNF", "JUNB",
      "FOSB", "FOSL2",
      "HOMER1", "NR4A2",
      "NPAS4", "MEF2C")
  
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
        each = 8)
  df_2_plot$sex <-
    rep(df_design$covar_table_1hr$sex,
        each = 8)
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
  ggtitle(paste0("GABA 1vs0 hr",
                ", 100 lines"))

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

df_writeout <-
  as.data.frame(t(interim_matrix))
df_writeout$cell_line <-
  cell_line_type_specific
df_writeout$sex <-
  df_design$covar_table_0hr$sex
df_writeout <-
  as.data.frame(t(df_writeout))
colnames(df_writeout) <-
  cell_line_type_specific
  # str_split(string = colnames(df_writeout),
  #           pattern = "_",
  #           simplify = T)[, 1]
write.table(df_writeout,
            file = "npglut_6vs0hr_all_log2FC.txt",
            quote = F, sep = "\t",
            row.names = F, col.names = T)


df_2_plot <-
  melt(interim_matrix)
df_2_plot$cell_line <-
  rep(cell_lines,
      each = 8)
df_2_plot$sex <-
  rep(df_design$covar_table_0hr$sex,
      each = 8)

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
  ggtitle(paste0("GABA 6vs0 hr",
                ", 100 lines"))

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

df_test <- voom_list_post_sva[['6hr']]$E
df_test[rownames(df_test) == "BDNF", ]
