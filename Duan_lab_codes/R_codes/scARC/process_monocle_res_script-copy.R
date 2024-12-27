# Chuxuan Li 06/02/2022
# Process the results from Monocle fit_models() for each of the cell types
# Run this script in conda directly

library(monocle3)

library(dplyr)
library(stringr)

library(ggplot2)
library(RColorBrewer)
library(ggrepel)

# # GABA ####
# print("load GABA")
# load("/nvmefs/scARC_Duan_018/GABA_fit_coefs.RData")
# table <- as.data.frame(GABA_fit_coefs)
# rm(GABA_fit_coefs)
# # full list
# terms_1hr <- table %>% filter(term == "time.ident1hr") %>%
#   select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
# write.table(terms_1hr, file = paste0("GABA", "_1v0_full_DEG_res.csv"), quote = F,
#            sep = ",", col.names = T)
#
# terms_6hr <- table %>% filter(term == "time.ident6hr") %>%
#   select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
# write.table(terms_6hr, file = paste0("GABA", "_6v0_full_DEG_res.csv"), quote = F,
#              sep = ",", col.names = T)
#
# # volcano plot on full list
# print("making volcano plots")
# terms_1hr$significance <- "nonsignificant"
# terms_1hr$significance[terms_1hr$q_value < 0.05 & terms_1hr$estimate > 0] <- "up"
# terms_1hr$significance[terms_1hr$q_value < 0.05 & terms_1hr$estimate < 0] <- "down"
# terms_1hr$significance <- factor(terms_1hr$significance,
#                                  levels = c("up", "nonsignificant", "down"))
# terms_1hr$neg_log_pval <- (0 - log2(terms_1hr$p_value))
# terms_1hr$labelling <- ""
# for (i in c("FOS", "NPAS4", "IGF1", "NR4A1", "FOSB", "EGR1", "BDNF", "VGF")){
#   terms_1hr$labelling[terms_1hr$gene_short_name %in% i] <- i
# }
#
# terms_6hr$significance <- "nonsignificant"
# terms_6hr$significance[terms_6hr$q_value < 0.05 & terms_6hr$estimate > 0] <- "up"
# terms_6hr$significance[terms_6hr$q_value < 0.05 & terms_6hr$estimate < 0] <- "down"
# terms_6hr$significance <- factor(terms_6hr$significance, levels = c("up", "nonsignificant", "down"))
# unique(terms_6hr$significance)
# terms_6hr$neg_log_pval <- (0 - log2(terms_6hr$p_value))
# terms_6hr$labelling <- ""
# for (i in c("FOS", "NPAS4", "IGF1", "NR4A1", "FOSB", "EGR1", "BDNF", "VGF")){
#   terms_6hr$labelling[terms_6hr$gene_short_name %in% i] <- i
# }
#
# jpeg("GABA_6v0_volcano_plot.jpeg")
# p <- ggplot(data = terms_6hr,
#             aes(x = estimate,
#                 y = neg_log_pval,
#                 color = significance,
#                 label = labelling)) +
#   geom_point(size = 0.2) +
#   scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
#   theme_minimal() +
#   geom_text_repel(box.padding = unit(0.05, 'lines'),
#                   min.segment.length = 0,
#                   force = 2,
#                   max.overlaps = 8000,
#                   force_pull = 0.5,
#                   show.legend = F) +
#   ggtitle("GABA 6v0hr") +
#   xlim(c(-10, 10))
# print(p)
# dev.off()
#
# jpeg("GABA_1v0_volcano_plot.jpeg")
# p <- ggplot(data = terms_1hr,
#             aes(x = estimate,
#                 y = neg_log_pval,
#                 color = significance,
#                 label = labelling)) +
#   geom_point(size = 0.2) +
#   scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
#   theme_minimal() +
#   geom_text_repel(box.padding = unit(0.05, 'lines'),
#                   min.segment.length = 0,
#                   force = 2,
#                   max.overlaps = 8000,
#                   force_pull = 0.5,
#                   show.legend = F) +
#   ggtitle(paste0("GABA 1v0hr")) +
#   xlim(c(-10, 10))
# print(p)
# dev.off()
#
# # significant, down/up separated
# terms_1hr_up <- table %>% filter(term == "time.ident1hr") %>%
#   filter(q_value < 0.05 & estimate > 0) %>%
#   select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
# write.table(terms_1hr_up, file = paste0("GABA", "_1v0_upregulated_DEGs.csv"), quote = F,
#             sep = ",", col.names = T)
# terms_1hr_do <- table %>% filter(term == "time.ident1hr") %>%
#   filter(q_value < 0.05 & estimate < 0) %>%
#   select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
# write.table(terms_1hr_do, file = paste0("GABA", "_1v0_downregulated_DEGs.csv"), quote = F,
#             sep = ",", col.names = T)
#
# terms_6hr_up <- table %>% filter(term == "time.ident6hr") %>%
#   filter(q_value < 0.05 & estimate > 0) %>%
#   select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
# write.table(terms_6hr_up, file = paste0("GABA", "_6v0_upregulated_DEGs.csv"), quote = F,
#             sep = ",", col.names = T)
# terms_6hr_do <- table %>% filter(term == "time.ident6hr") %>%
#   filter(q_value < 0.05 & estimate < 0) %>%
#   select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
# write.table(terms_6hr_do, file = paste0("GABA", "_6v0_downregulated_DEGs.csv"), quote = F,
#             sep = ",", col.names = T)
#
# terms_1hr_filtered_early <- terms_1hr[terms_1hr$gene_short_name %in%
#                                         c("FOS", "NPAS4", "NR4A1", "FOSB", "EGR1"), ]
# terms_1hr_filtered_late <- terms_1hr[terms_1hr$gene_short_name %in% c("IGF1", "BDNF", "VGF"), ]
#
# terms_6hr_filtered_early <- terms_6hr[terms_6hr$gene_short_name %in%
#                                         c("FOS", "NPAS4", "NR4A1", "FOSB", "EGR1"), ]
# terms_6hr_filtered_late <- terms_6hr[terms_6hr$gene_short_name %in% c("IGF1", "BDNF", "VGF"), ]
#
# df_to_plot_GABA_early <- as.data.frame(rbind(terms_1hr_filtered_early,
#                                              terms_6hr_filtered_early))
# df_to_plot_GABA_early$cell.type <- "GABA"
# df_to_plot_GABA_early$time <- paste0(str_extract(df_to_plot_GABA_early$term, "[0:1:6]"),
#                                      "vs0hr")
# save(df_to_plot_GABA_early, file = "df_to_plot_GABA_early.RData")
#
#
# df_to_plot_GABA_late <- as.data.frame(rbind(terms_1hr_filtered_late,
#                                             terms_6hr_filtered_late))
# df_to_plot_GABA_late$cell.type <- "GABA"
# df_to_plot_GABA_late$time <- paste0(str_extract(df_to_plot_GABA_late$term, "[0:1:6]"),
#                                     "vs0hr")
# save(df_to_plot_GABA_late, file = "df_to_plot_GABA_late.RData")
#
# rm(list = ls())
# gc()

# nmglut ####
print("load nmglut")
load("/nvmefs/scARC_Duan_018/nmglut_fit_coefs.RData")
table <- as.data.frame(nmglut_fit_coefs)
rm(nmglut_fit_coefs)
# full list
terms_1hr <- table %>% filter(term == "time.ident1hr") %>%
  select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
write.table(terms_1hr, file = paste0("nmglut", "_1v0_all_DEGs.csv"), quote = F,
            sep = ",", col.names = T)

terms_6hr <- table %>% filter(term == "time.ident6hr") %>%
  select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
write.table(terms_6hr, file = paste0("nmglut", "_6v0_all_DEGs.csv"), quote = F,
            sep = ",", col.names = T)

# volcano plot on full list
print("volcano plots")
terms_1hr$significance <- "nonsignificant"
terms_1hr$significance[terms_1hr$q_value < 0.05 & terms_1hr$estimate > 0] <- "up"
terms_1hr$significance[terms_1hr$q_value < 0.05 & terms_1hr$estimate < 0] <- "down"
terms_1hr$significance <- factor(terms_1hr$significance,
                                 levels = c("up", "nonsignificant", "down"))
terms_1hr$neg_log_pval <- (0 - log2(terms_1hr$p_value))
terms_1hr$labelling <- ""
for (i in c("FOS", "NPAS4", "IGF1", "NR4A1", "FOSB", "EGR1", "BDNF", "VGF")){
  terms_1hr$labelling[terms_1hr$gene_short_name %in% i] <- i
}

terms_6hr$significance <- "nonsignificant"
terms_6hr$significance[terms_6hr$q_value < 0.05 & terms_6hr$estimate > 0] <- "up"
terms_6hr$significance[terms_6hr$q_value < 0.05 & terms_6hr$estimate < 0] <- "down"
terms_6hr$significance <- factor(terms_6hr$significance, levels = c("up", "nonsignificant", "down"))
unique(terms_6hr$significance)
terms_6hr$neg_log_pval <- (0 - log2(terms_6hr$p_value))
terms_6hr$labelling <- ""
for (i in c("FOS", "NPAS4", "IGF1", "NR4A1", "FOSB", "EGR1", "BDNF", "VGF")){
  terms_6hr$labelling[terms_6hr$gene_short_name %in% i] <- i
}

jpeg("nmglut_6v0_volcano_plot.jpeg")
p <- ggplot(data = terms_6hr,
            aes(x = estimate,
                y = neg_log_pval,
                color = significance,
                label = labelling)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
  theme_minimal() +
  geom_text_repel(box.padding = unit(0.05, 'lines'),
                  min.segment.length = 0,
                  force = 2,
                  max.overlaps = 8000,
                  force_pull = 0.5,
                  show.legend = F) +
  ggtitle("nmglut 6v0hr") +
  xlim(c(-10, 10))
print(p)
dev.off()

jpeg("nmglut_1v0_volcano_plot.jpeg")
p <- ggplot(data = terms_1hr,
            aes(x = estimate,
                y = neg_log_pval,
                color = significance,
                label = labelling)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
  theme_minimal() +
  geom_text_repel(box.padding = unit(0.05, 'lines'),
                  min.segment.length = 0,
                  force = 2,
                  max.overlaps = 8000,
                  force_pull = 0.5,
                  show.legend = F) +
  ggtitle("nmglut 1v0hr") +
  xlim(c(-10, 10))
print(p)
dev.off()

# significant, down/up separated
terms_1hr_up <- table %>% filter(term == "time.ident1hr") %>%
  filter(q_value < 0.05 & estimate > 0) %>%
  select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
write.table(terms_1hr_up, file = paste0("nmglut", "_1v0_upregulated_DEGs.csv"), quote = F,
            sep = ",", col.names = T)
terms_1hr_do <- table %>% filter(term == "time.ident1hr") %>%
  filter(q_value < 0.05 & estimate < 0) %>%
  select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
write.table(terms_1hr_do, file = paste0("nmglut", "_1v0_downregulated_DEGs.csv"), quote = F,
            sep = ",", col.names = T)

terms_6hr_up <- table %>% filter(term == "time.ident6hr") %>%
  filter(q_value < 0.05 & estimate > 0) %>%
  select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
write.table(terms_6hr_up, file = paste0("nmglut", "_6v0_upregulated_DEGs.csv"), quote = F,
            sep = ",", col.names = T)
terms_6hr_do <- table %>% filter(term == "time.ident6hr") %>%
  filter(q_value < 0.05 & estimate < 0) %>%
  select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
write.table(terms_6hr_do, file = paste0("nmglut", "_6v0_downregulated_DEGs.csv"), quote = F,
            sep = ",", col.names = T)

terms_1hr_filtered_early <- terms_1hr[terms_1hr$gene_short_name %in%
                                        c("FOS", "NPAS4", "NR4A1", "FOSB", "EGR1"), ]
terms_1hr_filtered_late <- terms_1hr[terms_1hr$gene_short_name %in% c("IGF1", "BDNF", "VGF"), ]

terms_6hr_filtered_early <- terms_6hr[terms_6hr$gene_short_name %in%
                                        c("FOS", "NPAS4", "NR4A1", "FOSB", "EGR1"), ]
terms_6hr_filtered_late <- terms_6hr[terms_6hr$gene_short_name %in% c("IGF1", "BDNF", "VGF"), ]

df_to_plot_nmglut_early <- as.data.frame(rbind(terms_1hr_filtered_early,
                                               terms_6hr_filtered_early))
df_to_plot_nmglut_early$cell.type <- "nmglut"
df_to_plot_nmglut_early$time <- paste0(str_extract(df_to_plot_nmglut_early$term, "[0:1:6]"),
                                       "vs0hr")
save(df_to_plot_nmglut_early, file = "df_to_plot_nmglut_early.RData")


df_to_plot_nmglut_late <- as.data.frame(rbind(terms_1hr_filtered_late,
                                              terms_6hr_filtered_late))
df_to_plot_nmglut_late$cell.type <- "nmglut"
df_to_plot_nmglut_late$time <- paste0(str_extract(df_to_plot_nmglut_late$term, "[0:1:6]"),
                                      "vs0hr")
save(df_to_plot_nmglut_late, file = "df_to_plot_nmglut_late.RData")

rm(list = ls())
gc()


# npglut ####
print("load npglut")
load("/nvmefs/scARC_Duan_018/npglut_fit_coefs.RData")
table <- as.data.frame(npglut_fit_coefs)
rm(npglut_fit_coefs)
# full list
terms_1hr <- table %>% filter(term == "time.ident1hr") %>%
  select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
write.table(terms_1hr, file = paste0("npglut", "_1v0_all_DEGs.csv"), quote = F,
            sep = ",", col.names = T)

terms_6hr <- table %>% filter(term == "time.ident6hr") %>%
  select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
write.table(terms_6hr, file = paste0("npglut", "_6v0_all_DEGs.csv"), quote = F,
            sep = ",", col.names = T)

# volcano plot on full list
print("volcano plots")
terms_1hr$significance <- "nonsignificant"
terms_1hr$significance[terms_1hr$q_value < 0.05 & terms_1hr$estimate > 0] <- "up"
terms_1hr$significance[terms_1hr$q_value < 0.05 & terms_1hr$estimate < 0] <- "down"
terms_1hr$significance <- factor(terms_1hr$significance,
                                 levels = c("up", "nonsignificant", "down"))
terms_1hr$neg_log_pval <- (0 - log2(terms_1hr$p_value))
terms_1hr$labelling <- ""
for (i in c("FOS", "NPAS4", "IGF1", "NR4A1", "FOSB", "EGR1", "BDNF", "VGF")){
  terms_1hr$labelling[terms_1hr$gene_short_name %in% i] <- i
}

terms_6hr$significance <- "nonsignificant"
terms_6hr$significance[terms_6hr$q_value < 0.05 & terms_6hr$estimate > 0] <- "up"
terms_6hr$significance[terms_6hr$q_value < 0.05 & terms_6hr$estimate < 0] <- "down"
terms_6hr$significance <- factor(terms_6hr$significance, levels = c("up", "nonsignificant", "down"))
unique(terms_6hr$significance)
terms_6hr$neg_log_pval <- (0 - log2(terms_6hr$p_value))
terms_6hr$labelling <- ""
for (i in c("FOS", "NPAS4", "IGF1", "NR4A1", "FOSB", "EGR1", "BDNF", "VGF")){
  terms_6hr$labelling[terms_6hr$gene_short_name %in% i] <- i
}

jpeg("npglut_6v0_volcano_plot.jpeg")
p <- ggplot(data = terms_6hr,
            aes(x = estimate,
                y = neg_log_pval,
                color = significance,
                label = labelling)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
  theme_minimal() +
  geom_text_repel(box.padding = unit(0.05, 'lines'),
                  min.segment.length = 0,
                  force = 2,
                  max.overlaps = 8000,
                  force_pull = 0.5,
                  show.legend = F) +
  ggtitle("npglut 6v0hr") +
  xlim(c(-10, 10))
print(p)
dev.off()

jpeg("npglut_1v0_volcano_plot.jpeg")
p <- ggplot(data = terms_1hr,
            aes(x = estimate,
                y = neg_log_pval,
                color = significance,
                label = labelling)) +
  geom_point(size = 0.2) +
  scale_color_manual(values = c("red3", "grey50", "steelblue3")) +
  theme_minimal() +
  geom_text_repel(box.padding = unit(0.05, 'lines'),
                  min.segment.length = 0,
                  force = 2,
                  max.overlaps = 8000,
                  force_pull = 0.5,
                  show.legend = F) +
  ggtitle(paste0("npglut 1v0hr")) +
  xlim(c(-10, 10))
print(p)
dev.off()

# significant, down/up separated
terms_1hr_up <- table %>% filter(term == "time.ident1hr") %>%
  filter(q_value < 0.05 & estimate > 0) %>%
  select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
write.table(terms_1hr_up, file = paste0("npglut", "_1v0_upregulated_DEGs.csv"), quote = F,
            sep = ",", col.names = T)
terms_1hr_do <- table %>% filter(term == "time.ident1hr") %>%
  filter(q_value < 0.05 & estimate < 0) %>%
  select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
write.table(terms_1hr_do, file = paste0("npglut", "_1v0_downregulated_DEGs.csv"), quote = F,
            sep = ",", col.names = T)

terms_6hr_up <- table %>% filter(term == "time.ident6hr") %>%
  filter(q_value < 0.05 & estimate > 0) %>%
  select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
write.table(terms_6hr_up, file = paste0("npglut", "_6v0_upregulated_DEGs.csv"), quote = F,
            sep = ",", col.names = T)
terms_6hr_do <- table %>% filter(term == "time.ident6hr") %>%
  filter(q_value < 0.05 & estimate < 0) %>%
  select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
write.table(terms_6hr_do, file = paste0("npglut", "_6v0_downregulated_DEGs.csv"), quote = F,
            sep = ",", col.names = T)

terms_1hr_filtered_early <- terms_1hr[terms_1hr$gene_short_name %in%
                                        c("FOS", "NPAS4", "NR4A1", "FOSB", "EGR1"), ]
terms_1hr_filtered_late <- terms_1hr[terms_1hr$gene_short_name %in% c("IGF1", "BDNF", "VGF"), ]

terms_6hr_filtered_early <- terms_6hr[terms_6hr$gene_short_name %in%
                                        c("FOS", "NPAS4", "NR4A1", "FOSB", "EGR1"), ]
terms_6hr_filtered_late <- terms_6hr[terms_6hr$gene_short_name %in% c("IGF1", "BDNF", "VGF"), ]

df_to_plot_npglut_early <- as.data.frame(rbind(terms_1hr_filtered_early,
                                               terms_6hr_filtered_early))
df_to_plot_npglut_early$cell.type <- "npglut"
df_to_plot_npglut_early$time <- paste0(str_extract(df_to_plot_npglut_early$term, "[0:1:6]"),
                                       "vs0hr")
save(df_to_plot_npglut_early, file = "df_to_plot_npglut_early.RData")


df_to_plot_npglut_late <- as.data.frame(rbind(terms_1hr_filtered_late,
                                              terms_6hr_filtered_late))
df_to_plot_npglut_late$cell.type <- "npglut"
df_to_plot_npglut_late$time <- paste0(str_extract(df_to_plot_npglut_late$term, "[0:1:6]"),
                                      "vs0hr")
save(df_to_plot_npglut_late, file = "df_to_plot_npglut_late.RData")

