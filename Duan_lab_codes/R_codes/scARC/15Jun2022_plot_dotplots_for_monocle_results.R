# Chuxuan Li 06/15/2022
# process the results from process_monocle_res_script: plot dot plots

library(dplyr)
library(stringr)

library(ggplot2)
library(RColorBrewer)
library(ggrepel)


rdatas <- list.files(path = "./monocle_DE/", pattern = "[e|y].RData", full.names = T)
for (d in rdatas){
  load(d)
}
# dotplot ####
df_to_plot_early <- rbind(df_to_plot_GABA_early, df_to_plot_nmglut_early, df_to_plot_npglut_early)
df_to_plot_late <- rbind(df_to_plot_GABA_late, df_to_plot_nmglut_late, df_to_plot_npglut_late)
df_to_plot_late$neg_log_pval[is.infinite(df_to_plot_late$neg_log_pval)] <- 500
ggplot(df_to_plot_early,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = 0 - log10(p_value),
           fill = estimate)) +
  geom_point(shape = 21) +
  scale_size(range = c(2, 5)) +
  ylab("cell type") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "red", "darkred")) +
  facet_grid(cols = vars(gene_short_name)) + 
  theme_bw() +
  labs(size = "-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("single cell differential gene expression\nearly response genes")


ggplot(df_to_plot_late,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = neg_log_pval,
           fill = estimate)) +
  geom_point(shape = 21) +
  scale_size(range = c(3, 6)) +
  ylab("cell type") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "lightsalmon", "red", "darkred")) +
  facet_grid(cols = vars(gene_short_name)) + 
  theme_bw() +
  labs(size = "-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("single cell differential gene expression\nlate response genes")
