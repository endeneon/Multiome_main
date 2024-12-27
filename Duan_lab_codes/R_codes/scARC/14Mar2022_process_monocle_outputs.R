# Chuxuan Li 03/14/2022
# Process the results from Monocle fit_models() for each of the cell types


library(monocle3)

library(dplyr)
library(stringr)

library(ggplot2)
library(RColorBrewer)
library(ggrepel)

library(future)

types <- c("GABA", "nmglut", "npglut", "NPC")
titles <- c("GABA", "NEFM- glut", "NEFM+ glut", "NPC")

setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v5_monocle/5line_output")

i <- 0
for (table in list(GABA_fit_coefs, nmglut_fit_coefs, npglut_fit_coefs, NPC_fit_coefs)){
  i <- i + 1
  table <- as.data.frame(table)
  # full list
  terms_1hr <- table %>% filter(term == "time.ident1hr") %>% 
    select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
  write.table(terms_1hr, file = paste0(types[i], "_1v0_all_DEGs.csv"), quote = F, 
              sep = ",", col.names = T)
  
  terms_6hr <- table %>% filter(term == "time.ident6hr") %>% 
    select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
  write.table(terms_6hr, file = paste0(types[i], "_6v0_all_DEGs.csv"), quote = F, 
              sep = ",", col.names = T)
  
  # significant, down/up separated
  terms_1hr_up <- table %>% filter(term == "time.ident1hr") %>% filter(q_value < 0.05 & estimate > 0) %>%
    select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
  write.table(terms_1hr_up, file = paste0(types[i], "_1v0_upregulated_DEGs.csv"), quote = F, 
              sep = ",", col.names = T)
  terms_1hr_do <- table %>% filter(term == "time.ident1hr") %>% filter(q_value < 0.05 & estimate < 0) %>%
    select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
  write.table(terms_1hr_do, file = paste0(types[i], "_1v0_downregulated_DEGs.csv"), quote = F, 
              sep = ",", col.names = T)
  
  terms_6hr_up <- table %>% filter(term == "time.ident6hr") %>% filter(q_value < 0.05 & estimate > 0) %>%
    select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
  write.table(terms_6hr_up, file = paste0(types[i], "_6v0_upregulated_DEGs.csv"), quote = F, 
              sep = ",", col.names = T)
  terms_6hr_do <- table %>% filter(term == "time.ident6hr") %>% filter(q_value < 0.05 & estimate < 0) %>%
    select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
  write.table(terms_6hr_do, file = paste0(types[i], "_6v0_downregulated_DEGs.csv"), quote = F, 
              sep = ",", col.names = T)
  
  terms_1hr_filtered <- terms_1hr[terms_1hr$gene_short_name %in% 
                                    c("FOS", "NPAS4", "BDNF", "VGF"), ]
  
  terms_6hr_filtered <- terms_6hr[terms_6hr$gene_short_name %in% 
                                    c("FOS", "NPAS4", "BDNF", "VGF"), ]
  print(i)
  if (i == 1){
    df_to_plot <- as.data.frame(rbind(terms_1hr_filtered, terms_6hr_filtered))
    df_to_plot$cell.type <- types[i]
    df_to_plot$time <- paste0(str_extract(df_to_plot$term, "[0:1:6]"), "vs0hr")
  } else {
    df_to_append <- as.data.frame(rbind(terms_1hr_filtered, terms_6hr_filtered))
    df_to_append$cell.type <- types[i]
    df_to_append$time <- paste0(str_extract(df_to_append$term, "[0:1:6]"), "vs0hr")
    df_to_plot <- rbind(df_to_plot, df_to_append)
  }
}

# dotplot ####
df_to_plot$p_value_for_plot <- df_to_plot$p_value + 1e-307
ggplot(df_to_plot,
       aes(x = as.factor(time),
           y = as.factor(cell.type),
           size = 0 - log10(p_value_for_plot),
           fill = estimate)) +
  geom_point(shape = 21) +
  ylab("cell type") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "darkred")) +
  facet_grid(cols = vars(gene_short_name)) + 
  theme_bw() +
  labs(size = "-log10(p-value)")

# volcano plots ####
j <- 0
for (table in list(GABA_fit_coefs, nmglut_fit_coefs, npglut_fit_coefs, NPC_fit_coefs)){
  j <- j + 1
  # full list
  terms_1hr <- table %>% filter(term == "time.ident1hr") %>% 
    select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
  terms_1hr$significance <- "nonsignificant"
  terms_1hr$significance[terms_1hr$q_value < 0.05 & terms_1hr$estimate > 0] <- "up"
  terms_1hr$significance[terms_1hr$q_value < 0.05 & terms_1hr$estimate < 0] <- "down"
  terms_1hr$significance <- factor(terms_1hr$significance, levels = c("up", "nonsignificant", "down"))
  terms_1hr$neg_log_pval <- (0 - log2(terms_1hr$p_value))
  terms_1hr$labelling <- ""
  for (i in c("FOS", "NPAS4", "IGF1", "NR4A1", "FOSB", "EGR1", "BDNF")){
    terms_1hr$labelling[terms_1hr$gene_short_name %in% i] <- i
  }
  
  terms_6hr <- table %>% filter(term == "time.ident6hr") %>% 
    select(gene_short_name, term, p_value, q_value, estimate) %>% arrange(estimate)
  terms_6hr$significance <- "nonsignificant"
  terms_6hr$significance[terms_6hr$q_value < 0.05 & terms_6hr$estimate > 0] <- "up"
  terms_6hr$significance[terms_6hr$q_value < 0.05 & terms_6hr$estimate < 0] <- "down"
  terms_6hr$significance <- factor(terms_6hr$significance, levels = c("up", "nonsignificant", "down"))
  unique(terms_6hr$significance)
  terms_6hr$neg_log_pval <- (0 - log2(terms_6hr$p_value))
  terms_6hr$labelling <- ""
  for (i in c("FOS", "NPAS4", "IGF1", "NR4A1", "FOSB", "EGR1", "BDNF")){
    terms_6hr$labelling[terms_6hr$gene_short_name %in% i] <- i
  }
  
  jpeg(paste0(types[j], "_6v0_volcano_plot.jpeg"))
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
    ggtitle(paste0(types[j], " 6v0hr")) +
    xlim(c(-10, 10))
  print(p)
  dev.off()
  
  jpeg(paste0(types[j], "_1v0_volcano_plot.jpeg"))
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
    ggtitle(paste0(types[j], " 1v0hr")) +
    xlim(c(-10, 10))
  print(p)
  dev.off()
}
