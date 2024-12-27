# 23 Feb 2019
# Siwei 

library(ggplot2)

# master_intersect_table <- as.data.frame(rbind(intersect_table_25_types, intersect_table_fb))
# master_summary_table <- as.data.frame(rbind(summary_25_types, summary_fb))

colnames(master_intersect_table) <- c("cell_type", "annotation_type", "length")
colnames(master_summary_table) <- c("annotation_type", "total_feature_length")
colnames(narrowpeaks_length) <- c("cell_type", "total_peak_length")

CN_plot <- log10(
  master_intersect_table[master_intersect_table$cell_type == "CN", 3] /
    (narrowpeaks_length$total_peak_length[1]*master_summary_table$total_feature_length/3234830000)
)

CN_bar <- as.data.frame(cbind(master_summary_table$annotation_type, CN_plot))
colnames(CN_bar) <- c("annotation_type", "CN_plot")

ggplot(CN_bar, aes(x = annotation_type, y = CN_plot)) +
  geom_col(fill = "darkblue") +
  scale_x_discrete(limits = 
                     c("TssA", "PromU", "PromD1", "PromD2", 
                       "Tx5", "Tx", "Tx3", "TxWk", 
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", 
                       "EnhAc", "DNase", "ZNF_Rpts", 
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh")) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  xlab("Enrichment of differentially accessible peaks in epigenetically annotated regions \nGlutamatergic neuron") +
  ylab("log10(fold enrichment)")

###
DN_plot <- log10(
  master_intersect_table[master_intersect_table$cell_type == "DN", 3] /
    (narrowpeaks_length$total_peak_length[2]*master_summary_table$total_feature_length/3234830000)
)

DN_bar <- as.data.frame(cbind(master_summary_table$annotation_type, DN_plot))
colnames(DN_bar) <- c("annotation_type", "DN_plot")

ggplot(DN_bar, aes(x = annotation_type, y = DN_plot)) +
  geom_col(fill = "darkgreen") +
  scale_x_discrete(limits = 
                     c("TssA", "PromU", "PromD1", "PromD2", 
                       "Tx5", "Tx", "Tx3", "TxWk", 
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", 
                       "EnhAc", "DNase", "ZNF_Rpts", 
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh")) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  xlab("Enrichment of differentially accessible peaks in epigenetically annotated regions \nDopaminergic neuron") +
  ylab("log10(fold enrichment)")

###
ips_plot <- log10(
  master_intersect_table[master_intersect_table$cell_type == "ips", 3] /
    (narrowpeaks_length$total_peak_length[4] * master_summary_table$total_feature_length/3234830000)
)


ips_bar <- as.data.frame(cbind(master_summary_table$annotation_type, ips_plot))
colnames(ips_bar) <- c("annotation_type", "ips_plot")

ggplot(ips_bar, aes(x = annotation_type, y = ips_plot)) +
  geom_col(fill = "darkred") +
  scale_x_discrete(limits = 
                     c("TssA", "PromU", "PromD1", "PromD2", 
                       "Tx5", "Tx", "Tx3", "TxWk", 
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", 
                       "EnhAc", "DNase", "ZNF_Rpts", 
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh")) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  xlab("Enrichment of differentially accessible peaks in epigenetically annotated regions \niPS cell") +
  ylab("log10(fold enrichment)")


#####
GA_plot <- log10(
  master_intersect_table[master_intersect_table$cell_type == "GA", 3] /
    (narrowpeaks_length$total_peak_length[3]*master_summary_table$total_feature_length/3234830000)
)

GA_bar <- as.data.frame(cbind(master_summary_table$annotation_type, GA_plot))
colnames(GA_bar) <- c("annotation_type", "GA_plot")

ggplot(GA_bar, aes(x = annotation_type, y = GA_plot)) +
  geom_col(fill = "orange4") +
  scale_x_discrete(limits = 
                     c("TssA", "PromU", "PromD1", "PromD2", 
                       "Tx5", "Tx", "Tx3", "TxWk", 
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", 
                       "EnhAc", "DNase", "ZNF_Rpts", 
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh")) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  xlab("Enrichment of differentially accessible peaks in epigenetically annotated regions \nGABAergic neuron") +
  ylab("log10(fold enrichment)")

#####
NPC_plot <- log10(
  master_intersect_table[master_intersect_table$cell_type == "NSC", 3] /
    (narrowpeaks_length$total_peak_length[3]*master_summary_table$total_feature_length/3234830000)
)

NPC_bar <- as.data.frame(cbind(master_summary_table$annotation_type, NPC_plot))
colnames(NPC_bar) <- c("annotation_type", "NPC_plot")

ggplot(NPC_bar, aes(x = annotation_type, y = NPC_plot)) +
  geom_col(fill = "slategray4") +
  scale_x_discrete(limits = 
                     c("TssA", "PromU", "PromD1", "PromD2", 
                       "Tx5", "Tx", "Tx3", "TxWk", 
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", 
                       "EnhAc", "DNase", "ZNF_Rpts", 
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh")) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  xlab("Enrichment of differentially accessible peaks in epigenetically annotated regions \nNeural progenitor cell") +
  ylab("log10(fold enrichment)")
