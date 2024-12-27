# make peak plot

library(ggplot2)
library(ggrepel)

# by time point ####

ident_list <- c("GABA_1hr", "GABA_0hr",
                "GABA_6hr", "GABA_0hr",
                "GABA_6hr", "GABA_1hr",
                
                "NmCp_glut_1hr", "NmCp_glut_0hr", 
                "NmCp_glut_6hr", "NmCp_glut_0hr", 
                "NmCp_glut_6hr", "NmCp_glut_1hr",
                
                "NpCm_glut_1hr", "NpCm_glut_0hr", 
                "NpCm_glut_6hr", "NpCm_glut_0hr",
                "NpCm_glut_6hr", "NpCm_glut_1hr"
)

# load file
list_of_files <- list.files(path = "/nvmefs/scARC_Duan_018/R_scARC_Duan_018/most_updated_obj_da_peaks_by_time/bed/annotated/annotated_added_peak_id_output/",
                            full.names = T,
                            pattern = "*.txt")

# test
df_to_plot <- read_delim("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/most_updated_obj_da_peaks_by_time/bed/annotated/annotated_added_peak_id_output/NmCp_glut_6hr_NmCp_glut_0hr_uid.txt", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)


df_to_plot$significance <- "nonsig"
df_to_plot$significance[df_to_plot$pval_adj < 0.05 & df_to_plot$logFC > 0] <- "pos"
df_to_plot$significance[df_to_plot$pval_adj < 0.05 & df_to_plot$logFC < 0] <- "neg"
df_to_plot$neg_log_pval <- (0 - log10(df_to_plot$pVal))
df_to_plot$labelling <- ""
for (i in c("NPAS4", "BDNF", "FOS", "VGF", "EGR1", "IGF1", "NRN1", "PNOC")){
  df_to_plot$labelling[df_to_plot$`Gene Name` %in% i] <- i
}
unique(df_to_plot$labelling)

# df_to_plot$labelling <- ifelse(df_to_plot$`Gene Name` %in% c("NPAS4", "BDNF", "FOS", "VGF", "EGR1", "IGF1", "NRN1", "PNOC"),
                               # )

ggplot(data = df_to_plot,
       aes(x = logFC, 
           y = neg_log_pval, 
           color = significance,
           label = labelling)) + 
  geom_point() +
  scale_color_manual(values = colors) +
  theme_minimal() +
  geom_text_repel(box.padding = unit(0.05, 'lines'),
                  min.segment.length = 0,
                  force = 5,
                  max.overlaps = 10000,
                  force_pull = 0.5)


# now write a loop to plot all plots

for (i in 1:length(list_of_files)){
  f <- list_of_files[i]
  # get the correct indices
  ind1 <- 2 * i - 1
  ident1 <- ident_list[ind1]
  ind2 <- 2 * i
  ident2 <- ident_list[ind2]
  
  file_name <- paste0("../new_peak_set_plots/da_peaks_by_time_volcano_plots_labeled/", # pwd: /codes
                      ident1, "_",
                      ident2, "_da_peaks_volcano_plot.pdf"
  )
  
  print(file_name)
  
  # load file
  df_to_plot <- read_delim(f, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # add columns
  df_to_plot$significance <- "nonsig"
  df_to_plot$significance[df_to_plot$pval_adj < 0.05 & df_to_plot$logFC > 0] <- "pos"
  df_to_plot$significance[df_to_plot$pval_adj < 0.05 & df_to_plot$logFC < 0] <- "neg"
  df_to_plot$neg_log_pval <- (0 - log10(df_to_plot$pVal))
  df_to_plot$labelling <- ""
  for (i in c("NPAS4", "BDNF", "FOS", "VGF", "EGR1", "IGF1", "NRN1", "PNOC")){
    df_to_plot$labelling[df_to_plot$`Gene Name` %in% i] <- i
  }
  unique(df_to_plot$labelling)
  
  # plot
  pdf(file = file_name)
  p <- ggplot(data = df_to_plot,
         aes(x = logFC, 
             y = neg_log_pval, 
             color = significance,
             label = labelling)) + 
    geom_point() +
    scale_color_manual(values = colors) +
    theme_minimal() +
    geom_text_repel(box.padding = unit(0.05, 'lines'),
                    min.segment.length = 0,
                    force = 5,
                    max.overlaps = 10000,
                    force_pull = 0.5)
  print(p)
  dev.off()
}


# by cell type ####

ident_list <- c("GABA_0hr", "NmCp_glut_0hr",
                "GABA_1hr", "NmCp_glut_1hr",
                "GABA_6hr", "NmCp_glut_6hr",
                
                "GABA_0hr", "NpCm_glut_0hr",
                "GABA_1hr", "NpCm_glut_1hr",
                "GABA_6hr", "NpCm_glut_6hr",
                
                "NmCp_glut_0hr", "NpCm_glut_0hr",
                "NmCp_glut_1hr", "NpCm_glut_1hr",
                "NmCp_glut_6hr", "NpCm_glut_6hr")

# load file
list_of_files <- list.files(path = "/nvmefs/scARC_Duan_018/R_scARC_Duan_018/most_updated_obj_da_peaks_by_celltype/bed/annotated/annotated_added_peak_id_output/",
                            full.names = T,
                            pattern = "*.txt")


for (i in 1:length(list_of_files)){
  f <- list_of_files[i]
  # get the correct indices
  ind1 <- 2 * i - 1
  ident1 <- ident_list[ind1]
  ind2 <- 2 * i
  ident2 <- ident_list[ind2]
  
  file_name <- paste0("../new_peak_set_plots/da_peaks_by_celltype_volcano_plots_labeled/", # pwd: /codes
                      ident1, "_",
                      ident2, "_da_peaks_volcano_plot.pdf"
  )
  
  print(file_name)
  
  # load file
  df_to_plot <- read_delim(f, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # add columns
  df_to_plot$significance <- "nonsig"
  df_to_plot$significance[df_to_plot$pval_adj < 0.05 & df_to_plot$logFC > 0] <- "pos"
  df_to_plot$significance[df_to_plot$pval_adj < 0.05 & df_to_plot$logFC < 0] <- "neg"
  df_to_plot$neg_log_pval <- (0 - log10(df_to_plot$pVal))
  df_to_plot$labelling <- ""
  for (i in c("NPAS4", "BDNF", "FOS", "VGF", "EGR1", "IGF1", "NRN1", "PNOC")){
    df_to_plot$labelling[df_to_plot$`Gene Name` %in% i] <- i
  }
  unique(df_to_plot$labelling)
  
  # plot
  pdf(file = file_name)
  p <- ggplot(data = df_to_plot,
              aes(x = logFC, 
                  y = neg_log_pval, 
                  color = significance,
                  label = labelling)) + 
    geom_point() +
    scale_color_manual(values = colors) +
    theme_minimal() +
    geom_text_repel(box.padding = unit(0.05, 'lines'),
                    min.segment.length = 0,
                    force = 5,
                    max.overlaps = 10000,
                    force_pull = 0.5)
  print(p)
  dev.off()
}
