#######

######
SNP_unique_plot <- log10(
  main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.05", 3] /
    (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000)
)$length
# 
# SNP_unique_plot[32] <- 
#   log10(sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[9:19]) /
#           (SNP_length$total_SNP_length[1] * 
#              sum(master_summary_table$total_feature_length[9:19])/3234830000))#, 


SNP_unique_bar <- data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot), stringsAsFactors = F)
# SNP_unique_bar$`master_summary_table$annotation_type` <- as.character(as.vector(SNP_unique_bar$`master_summary_table$annotation_type`))
SNP_unique_bar[32, ] <- c("Enh_All", 
                          log10(sum(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.05", 3]$length[9:19]) /
                                  (SNP_length$total_SNP_length[1] * 
                                     sum(master_summary_table$total_feature_length[9:19])/3234830000))#, 
                          # "all_Enh_Added"
)
colnames(SNP_unique_bar) <- c("annotation_type", "SNP_unique_plot")

ggplot(SNP_unique_bar, aes(x = annotation_type, 
                           y = as.numeric(SNP_unique_bar$SNP_unique_plot))) +
  geom_col(fill = "royalblue3") +
  scale_x_discrete(limits = 
                     c("TssA", "PromU", "PromD1", "PromD2", 
                       "Tx5", "Tx", "Tx3", "TxWk", 
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", 
                       "EnhAc", "DNase", "ZNF_Rpts", 
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh", 
                       "H3K27acF", "H3K4me2F", "H3K27acO", "H3K4me2O", 
                       "Enh_All")) +
  geom_text(aes(label = as.vector(c(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.05", 3]$length, 
                                    sum(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.05", 3]$length[9:19])))),
            hjust = -0.25,
            # y = log10(c(
            #   (main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3] /
            #     (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000))$length, 6) + 0.1),
            
            # hjust = "inward",
            angle = 90) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  xlab("Enrichment of iPS-neural shared ASoC SNPs (FDR=0.05) in epigenetically annotated regions") +
  ylab("log10(fold enrichment)")

######

SNP_length$cell_type
SNP_unique_plot <- log10(
  main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.1", 3] /
    (SNP_length$total_SNP_length[2]*master_summary_table$total_feature_length/3234830000)
)$length
# 
# SNP_unique_plot[32] <- 
#   log10(sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[9:19]) /
#           (SNP_length$total_SNP_length[1] * 
#              sum(master_summary_table$total_feature_length[9:19])/3234830000))#, 


SNP_unique_bar <- data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot), stringsAsFactors = F)
# SNP_unique_bar$`master_summary_table$annotation_type` <- as.character(as.vector(SNP_unique_bar$`master_summary_table$annotation_type`))
SNP_unique_bar[32, ] <- c("Enh_All", 
                          log10(sum(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.1", 3]$length[9:19]) /
                                  (SNP_length$total_SNP_length[2] * 
                                     sum(master_summary_table$total_feature_length[9:19])/3234830000))#, 
                          # "all_Enh_Added"
)
colnames(SNP_unique_bar) <- c("annotation_type", "SNP_unique_plot")

ggplot(SNP_unique_bar, aes(x = annotation_type, 
                           y = as.numeric(SNP_unique_bar$SNP_unique_plot))) +
  geom_col(fill = "royalblue3") +
  scale_x_discrete(limits = 
                     c("TssA", "PromU", "PromD1", "PromD2", 
                       "Tx5", "Tx", "Tx3", "TxWk", 
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", 
                       "EnhAc", "DNase", "ZNF_Rpts", 
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh", 
                       "H3K27acF", "H3K4me2F", "H3K27acO", "H3K4me2O", 
                       "Enh_All")) +
  geom_text(aes(label = as.vector(c(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.1", 3]$length, 
                                    sum(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.1", 3]$length[9:19])))),
            hjust = -0.25,
            # y = log10(c(
            #   (main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3] /
            #     (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000))$length, 6) + 0.1),
            
            # hjust = "inward",
            angle = 90) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  xlab("Enrichment of iPS-neural shared ASoC SNPs (FDR=0.1) in epigenetically annotated regions") +
  ylab("log10(fold enrichment)")

######

SNP_length$cell_type
SNP_unique_plot <- log10(
  main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.2", 3] /
    (SNP_length$total_SNP_length[3]*master_summary_table$total_feature_length/3234830000)
)$length
# 
# SNP_unique_plot[32] <- 
#   log10(sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[9:19]) /
#           (SNP_length$total_SNP_length[1] * 
#              sum(master_summary_table$total_feature_length[9:19])/3234830000))#, 


SNP_unique_bar <- data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot), stringsAsFactors = F)
# SNP_unique_bar$`master_summary_table$annotation_type` <- as.character(as.vector(SNP_unique_bar$`master_summary_table$annotation_type`))
SNP_unique_bar[32, ] <- c("Enh_All", 
                          log10(sum(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.2", 3]$length[9:19]) /
                                  (SNP_length$total_SNP_length[3] * 
                                     sum(master_summary_table$total_feature_length[9:19])/3234830000))#, 
                          # "all_Enh_Added"
)
colnames(SNP_unique_bar) <- c("annotation_type", "SNP_unique_plot")

ggplot(SNP_unique_bar, aes(x = annotation_type, 
                           y = as.numeric(SNP_unique_bar$SNP_unique_plot))) +
  geom_col(fill = "royalblue3") +
  scale_x_discrete(limits = 
                     c("TssA", "PromU", "PromD1", "PromD2", 
                       "Tx5", "Tx", "Tx3", "TxWk", 
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", 
                       "EnhAc", "DNase", "ZNF_Rpts", 
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh", 
                       "H3K27acF", "H3K4me2F", "H3K27acO", "H3K4me2O", 
                       "Enh_All")) +
  geom_text(aes(label = as.vector(c(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.2", 3]$length, 
                                    sum(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.2", 3]$length[9:19])))),
            hjust = -0.25,
            # y = log10(c(
            #   (main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3] /
            #     (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000))$length, 6) + 0.1),
            
            # hjust = "inward",
            angle = 90) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  xlab("Enrichment of iPS-neural shared ASoC SNPs (FDR=0.2) in epigenetically annotated regions") +
  ylab("log10(fold enrichment)")

######

SNP_length$cell_type
SNP_unique_plot <- log10(
  main_SNP_table[main_SNP_table$cell_type == "iPS_only_specific", 3] /
    (SNP_length$total_SNP_length[4]*master_summary_table$total_feature_length/3234830000)
)$length
# 
# SNP_unique_plot[32] <- 
#   log10(sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[9:19]) /
#           (SNP_length$total_SNP_length[1] * 
#              sum(master_summary_table$total_feature_length[9:19])/3234830000))#, 


SNP_unique_bar <- data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot), stringsAsFactors = F)
# SNP_unique_bar$`master_summary_table$annotation_type` <- as.character(as.vector(SNP_unique_bar$`master_summary_table$annotation_type`))
SNP_unique_bar[32, ] <- c("Enh_All", 
                          log10(sum(main_SNP_table[main_SNP_table$cell_type == "iPS_only_specific", 3]$length[9:19]) /
                                  (SNP_length$total_SNP_length[4] * 
                                     sum(master_summary_table$total_feature_length[9:19])/3234830000))#, 
                          # "all_Enh_Added"
)
colnames(SNP_unique_bar) <- c("annotation_type", "SNP_unique_plot")

ggplot(SNP_unique_bar, aes(x = annotation_type, 
                           y = as.numeric(SNP_unique_bar$SNP_unique_plot))) +
  geom_col(fill = "royalblue3") +
  scale_x_discrete(limits = 
                     c("TssA", "PromU", "PromD1", "PromD2", 
                       "Tx5", "Tx", "Tx3", "TxWk", 
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", 
                       "EnhAc", "DNase", "ZNF_Rpts", 
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh", 
                       "H3K27acF", "H3K4me2F", "H3K27acO", "H3K4me2O", 
                       "Enh_All")) +
  geom_text(aes(label = as.vector(c(main_SNP_table[main_SNP_table$cell_type == "iPS_only_specific", 3]$length, 
                                    sum(main_SNP_table[main_SNP_table$cell_type == "iPS_only_specific", 3]$length[9:19])))),
            hjust = -0.25,
            # y = log10(c(
            #   (main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3] /
            #     (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000))$length, 6) + 0.1),
            
            # hjust = "inward",
            angle = 90) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  xlab("Enrichment of iPS-specific shared ASoC SNPs in epigenetically annotated regions") +
  ylab("log10(fold enrichment)")

######

SNP_length$cell_type
SNP_unique_plot <- log10(
  main_SNP_table[main_SNP_table$cell_type == "neuronal_4_cell_types", 3] /
    (SNP_length$total_SNP_length[5]*master_summary_table$total_feature_length/3234830000)
)$length
# 
# SNP_unique_plot[32] <- 
#   log10(sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[9:19]) /
#           (SNP_length$total_SNP_length[1] * 
#              sum(master_summary_table$total_feature_length[9:19])/3234830000))#, 


SNP_unique_bar <- data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot), stringsAsFactors = F)
# SNP_unique_bar$`master_summary_table$annotation_type` <- as.character(as.vector(SNP_unique_bar$`master_summary_table$annotation_type`))
SNP_unique_bar[32, ] <- c("Enh_All", 
                          log10(sum(main_SNP_table[main_SNP_table$cell_type == "neuronal_4_cell_types", 3]$length[9:19]) /
                                  (SNP_length$total_SNP_length[5] * 
                                     sum(master_summary_table$total_feature_length[9:19])/3234830000))#, 
                          # "all_Enh_Added"
)
colnames(SNP_unique_bar) <- c("annotation_type", "SNP_unique_plot")

ggplot(SNP_unique_bar, aes(x = annotation_type, 
                           y = as.numeric(SNP_unique_bar$SNP_unique_plot))) +
  geom_col(fill = "royalblue3") +
  scale_x_discrete(limits = 
                     c("TssA", "PromU", "PromD1", "PromD2", 
                       "Tx5", "Tx", "Tx3", "TxWk", 
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", 
                       "EnhAc", "DNase", "ZNF_Rpts", 
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh", 
                       "H3K27acF", "H3K4me2F", "H3K27acO", "H3K4me2O", 
                       "Enh_All")) +
  geom_text(aes(label = as.vector(c(main_SNP_table[main_SNP_table$cell_type == "neuronal_4_cell_types", 3]$length, 
                                    sum(main_SNP_table[main_SNP_table$cell_type == "neuronal_4_cell_types", 3]$length[9:19])))),
            hjust = -0.25,
            # y = log10(c(
            #   (main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3] /
            #     (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000))$length, 6) + 0.1),
            
            # hjust = "inward",
            angle = 90) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  xlab("Enrichment of neuronal-specific shared ASoC SNPs in epigenetically annotated regions") +
  ylab("log10(fold enrichment)")

