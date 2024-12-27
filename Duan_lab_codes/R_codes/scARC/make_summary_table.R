# 5 line summary

library(Seurat)
library(Signac)

library(ggplot2)

DefaultAssay(obj_complete) <- "RNA"

summary(Idents(obj_complete))

obj_complete@assays$ATAC

cell_type <- c("NEFM+/CUX2- glut","NEFM-/CUX2+ glut","GABA", "NPC")

time_ident <- unique(obj_complete$time.ident)
cell_line <- unique(obj_complete$)

i <- 1
j <- 1

for (j in 1:length(time_ident)) {
  for (i in 1:length(cell_type)) {
    print(paste(i, j))
    if (i == 1 & j == 1) {
      count_df <- subset(x = obj_complete,
                         subset = ((time.ident == time_ident[j]) &
                                     (fine.cell.type == cell_type[i])))
      count_df <- ncol(count_df@assays$ATAC@counts)
      # 
      colnames_vector <- paste(cell_type[i], time_ident[j], sep = "_")
      summary_df <- data.frame(count_df,
                               stringsAsFactors = F)
      colnames(summary_df) <- colnames_vector
    } else {
      count_df <- subset(x = obj_complete,
                         subset = ((time.ident == time_ident[j]) &
                                     (fine.cell.type == cell_type[i])))
      count_df <- ncol(count_df@assays$ATAC@counts)
      summary_df <- data.frame(cbind(summary_df,
                                     count_df),
                               stringsAsFactors = F)
      colnames_vector <- c(colnames_vector,
                           paste(cell_type[i], time_ident[j], sep = "_"))
      colnames(summary_df) <- colnames_vector
      
    }
  }
}

summary(subset(x = obj_complete,
               subset = (time.ident = "0hr" &
                           fine.cell.type = "GABA")))


library(readr)
library(ggplot2)
library(RColorBrewer)


df_to_plot <-
  read_csv("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/cell_counts_per_cell_type_0hr.csv")

df_to_plot <-
  read_csv("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/cell_counts_per_cell_type_1hr.csv")

df_to_plot <-
  read_csv("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/5line_codes_and_rdata/cell_counts_per_cell_type_6hr.csv")


ggplot(df_to_plot, 
       aes(x = cell.line,
           y = counts,
           fill = cell.type)) +
  geom_bar(position = "fill",
           stat = "identity") +
  scale_fill_manual(values = brewer.pal(7, "Dark2")) +
  theme_classic() +
  ggtitle("6 hr")

