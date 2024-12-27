#######

######
SNP_unique_plot <- log10(
  main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.05", 3] /
    (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000)
)$length
# 
# SNP_unique_plot[32] <- 
#   log10(sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[13:18]) /
#           (SNP_length$total_SNP_length[1] * 
#              sum(master_summary_table$total_feature_length[13:18])/3234830000))#, 


SNP_unique_bar <- data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot), stringsAsFactors = F)
# SNP_unique_bar$`master_summary_table$annotation_type` <- as.character(as.vector(SNP_unique_bar$`master_summary_table$annotation_type`))
SNP_unique_bar[32, ] <- c("Enh_All", 
                          log10(sum(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.05", 3]$length[13:18]) /
                                  (SNP_length$total_SNP_length[1] * 
                                     sum(master_summary_table$total_feature_length[13:18])/3234830000))#, 
                          # "all_Enh_Added"
)
colnames(SNP_unique_bar) <- c("annotation_type", "SNP_unique_plot")

SNP_unique_bar$GREAT_p <- 0

for (i in 1:nrow(SNP_unique_bar-1)) {
  SNP_unique_bar$GREAT_p[i] <- (0-log10(binom.test(x = main_SNP_table$length[main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.05"]],
                                                      n = 4411, #SNP_unique_table$n_total_features_length[i],
                                                      p = SNP_unique_table$total_feature_length[i]/3234830000,
                                                      alternative = "greater")$p.value))
}

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
                                    sum(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.05", 3]$length[13:18])))),
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
#   log10(sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[13:18]) /
#           (SNP_length$total_SNP_length[1] * 
#              sum(master_summary_table$total_feature_length[13:18])/3234830000))#, 


SNP_unique_bar <- data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot), stringsAsFactors = F)
# SNP_unique_bar$`master_summary_table$annotation_type` <- as.character(as.vector(SNP_unique_bar$`master_summary_table$annotation_type`))
SNP_unique_bar[32, ] <- c("Enh_All", 
                          log10(sum(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.1", 3]$length[13:18]) /
                                  (SNP_length$total_SNP_length[2] * 
                                     sum(master_summary_table$total_feature_length[13:18])/3234830000))#, 
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
                                    sum(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.1", 3]$length[13:18])))),
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
#   log10(sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[13:18]) /
#           (SNP_length$total_SNP_length[1] * 
#              sum(master_summary_table$total_feature_length[13:18])/3234830000))#, 


SNP_unique_bar <- data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot), stringsAsFactors = F)
# SNP_unique_bar$`master_summary_table$annotation_type` <- as.character(as.vector(SNP_unique_bar$`master_summary_table$annotation_type`))
SNP_unique_bar[32, ] <- c("Enh_All", 
                          log10(sum(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.2", 3]$length[13:18]) /
                                  (SNP_length$total_SNP_length[3] * 
                                     sum(master_summary_table$total_feature_length[13:18])/3234830000))#, 
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
                                    sum(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.2", 3]$length[13:18])))),
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
#   log10(sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[13:18]) /
#           (SNP_length$total_SNP_length[1] * 
#              sum(master_summary_table$total_feature_length[13:18])/3234830000))#, 


SNP_unique_bar <- data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot), stringsAsFactors = F)
# SNP_unique_bar$`master_summary_table$annotation_type` <- as.character(as.vector(SNP_unique_bar$`master_summary_table$annotation_type`))
SNP_unique_bar[32, ] <- c("Enh_All", 
                          log10(sum(main_SNP_table[main_SNP_table$cell_type == "iPS_only_specific", 3]$length[13:18]) /
                                  (SNP_length$total_SNP_length[4] * 
                                     sum(master_summary_table$total_feature_length[13:18])/3234830000))#, 
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
                                    sum(main_SNP_table[main_SNP_table$cell_type == "iPS_only_specific", 3]$length[13:18])))),
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
#   log10(sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[13:18]) /
#           (SNP_length$total_SNP_length[1] * 
#              sum(master_summary_table$total_feature_length[13:18])/3234830000))#, 


SNP_unique_bar <- data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot), stringsAsFactors = F)
# SNP_unique_bar$`master_summary_table$annotation_type` <- as.character(as.vector(SNP_unique_bar$`master_summary_table$annotation_type`))
SNP_unique_bar[32, ] <- c("Enh_All", 
                          log10(sum(main_SNP_table[main_SNP_table$cell_type == "neuronal_4_cell_types", 3]$length[13:18]) /
                                  (SNP_length$total_SNP_length[5] * 
                                     sum(master_summary_table$total_feature_length[13:18])/3234830000))#, 
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
                                    sum(main_SNP_table[main_SNP_table$cell_type == "neuronal_4_cell_types", 3]$length[13:18])))),
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

####

SNP_unique_table <- as.data.frame(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.2", 1:3])
SNP_unique_table[32, ] <- c("iPS_neural_FDR_shared_q_0.2", "Enh_All", as.numeric(sum(SNP_unique_table$length[13:18])))
SNP_unique_table$length <- as.numeric(SNP_unique_table$length)

master_summary_table[32, ] <- c("Enh_All", sum(as.numeric(master_summary_table$total_feature_length[13:18])))
SNP_unique_table$total_feature_length <- master_summary_table$total_feature_length
# SNP_unique_table$total_feature_length <- SNP_length$total_SNP_length[3]
# SNP_length$total_SNP_length[3]
SNP_unique_table$GREATenrichment <- 1
i <- 1
# j <- 1
for (i in 1:nrow(SNP_unique_table)) {
  print(paste("i=", i, collapse = ""))
  print(paste(SNP_unique_table$annotation_type[i]))
  SNP_unique_table$GREATenrichment[i] <- (0 - log10(binom.test(x = as.numeric(SNP_unique_table$length[i]),
                                                             n = sum(SNP_unique_table$length), #357, #4411, #357# #SNP_unique_table$n_total_features_length[i],
                                                             p = as.numeric(as.character(SNP_unique_table$total_feature_length[i]))/3234830000,
                                                             alternative = "greater")$p.value))
}

ggplot(SNP_unique_table, aes(x = annotation_type, y = GREATenrichment)) +
  geom_col(fill = "turquoise4") +
  # scale_x_discrete(limits =
  #                    c("TssA", "PromU", "PromD1", "PromD2",
  #                      "Tx5", "Tx", "Tx3", "TxWk",
  #                      "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2",
  #                      "EnhAc", "DNase", "ZNF_Rpts",
  #                      "Het", "PromP", "PromBiv", "ReprPC", "Quies",
  #                      "forebrain_enh", "non_forebrain_enh")) +
  scale_x_discrete(limits = SNP_unique_table$annotation_type) +
  # geom_text(aes(label = SNP_unique_table$kpi_SNP_count), 
  #           vjust = -0.25) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  ylim(c(0, 300)) +
  # ylim(c(0, max(SNP_unique_table$GREATenrichment))) +
  xlab("GREAT Enrichment of iPS-neural shared ASoC SNPs (FDR=0.2) \nin epigenetically annotated regions") +
  # xlab("GREAT Enrichment of cell-type-specific ASoC SNPs \nin epigenetically annotated regions") +
  ylab("-log10(Pvalue) fold enrichment")

####

SNP_unique_table <- as.data.frame(main_SNP_table[main_SNP_table$cell_type == "neuronal_4_cell_types", 1:3])
SNP_unique_table[32, ] <- c("neuronal_4_cell_types", "Enh_All", as.numeric(sum(SNP_unique_table$length[13:18])))
SNP_unique_table$length <- as.numeric(SNP_unique_table$length)

master_summary_table[32, ] <- c("Enh_All", sum(as.numeric(master_summary_table$total_feature_length[13:18])))
SNP_unique_table$total_feature_length <- master_summary_table$total_feature_length
# SNP_unique_table$total_feature_length <- SNP_length$total_SNP_length[3]
# SNP_length$total_SNP_length[3]
SNP_unique_table$GREATenrichment <- 1
i <- 1
# j <- 1
for (i in 1:nrow(SNP_unique_table)) {
  print(paste("i=", i, collapse = ""))
  print(paste(SNP_unique_table$annotation_type[i]))
  SNP_unique_table$GREATenrichment[i] <- (0 - log10(binom.test(x = as.numeric(SNP_unique_table$length[i]),
                                                               n = sum(SNP_unique_table$length), #357, #4411, #357# #SNP_unique_table$n_total_features_length[i],
                                                               p = as.numeric(as.character(SNP_unique_table$total_feature_length[i]))/3234830000,
                                                               alternative = "greater")$p.value))
}

ggplot(SNP_unique_table, aes(x = annotation_type, y = GREATenrichment)) +
  geom_col(fill = "turquoise4") +
  # scale_x_discrete(limits =
  #                    c("TssA", "PromU", "PromD1", "PromD2",
  #                      "Tx5", "Tx", "Tx3", "TxWk",
  #                      "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2",
  #                      "EnhAc", "DNase", "ZNF_Rpts",
  #                      "Het", "PromP", "PromBiv", "ReprPC", "Quies",
  #                      "forebrain_enh", "non_forebrain_enh")) +
  scale_x_discrete(limits = SNP_unique_table$annotation_type) +
  # geom_text(aes(label = SNP_unique_table$kpi_SNP_count), 
  #           vjust = -0.25) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  ylim(c(0, 300)) +
  # ylim(c(0, max(SNP_unique_table$GREATenrichment))) +
  xlab("GREAT Enrichment of neuronal-specific ASoC SNPs \nin epigenetically annotated regions") +
  # xlab("GREAT Enrichment of cell-type-specific ASoC SNPs \nin epigenetically annotated regions") +
  ylab("-log10(Pvalue) fold enrichment")

####

SNP_unique_table <- as.data.frame(main_SNP_table[main_SNP_table$cell_type == "iPS_neural_FDR_shared_q_0.1", 1:3])
SNP_unique_table[32, ] <- c("iPS_neural_FDR_shared_q_0.1", "Enh_All", as.numeric(sum(SNP_unique_table$length[13:18])))
SNP_unique_table$length <- as.numeric(SNP_unique_table$length)

master_summary_table[32, ] <- c("Enh_All", sum(as.numeric(master_summary_table$total_feature_length[13:18])))
SNP_unique_table$total_feature_length <- master_summary_table$total_feature_length
# SNP_unique_table$total_feature_length <- SNP_length$total_SNP_length[3]
# SNP_length$total_SNP_length[3]
SNP_unique_table$GREATenrichment <- 1
i <- 1
# j <- 1
for (i in 1:nrow(SNP_unique_table)) {
  print(paste("i=", i, collapse = ""))
  print(paste(SNP_unique_table$annotation_type[i]))
  SNP_unique_table$GREATenrichment[i] <- (0 - log10(binom.test(x = as.numeric(SNP_unique_table$length[i]),
                                                               n = sum(SNP_unique_table$length), #357, #4411, #357# #SNP_unique_table$n_total_features_length[i],
                                                               p = as.numeric(as.character(SNP_unique_table$total_feature_length[i]))/3234830000,
                                                               alternative = "greater")$p.value))
  SNP_unique_table$p_GREAT[i] <- binom.test(x = as.numeric(SNP_unique_table$length[i]),
                                                               n = sum(SNP_unique_table$length), #357, #4411, #357# #SNP_unique_table$n_total_features_length[i],
                                                               p = as.numeric(as.character(SNP_unique_table$total_feature_length[i]))/3234830000,
                                                               alternative = "greater")$p.value
}

p.adjust(SNP_unique_table$p_GREAT, method = "fdr")


ggplot(SNP_unique_table, aes(x = annotation_type, y = GREATenrichment)) +
  geom_col(fill = "turquoise4") +
  # scale_x_discrete(limits =
  #                    c("TssA", "PromU", "PromD1", "PromD2",
  #                      "Tx5", "Tx", "Tx3", "TxWk",
  #                      "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2",
  #                      "EnhAc", "DNase", "ZNF_Rpts",
  #                      "Het", "PromP", "PromBiv", "ReprPC", "Quies",
  #                      "forebrain_enh", "non_forebrain_enh")) +
  scale_x_discrete(limits = SNP_unique_table$annotation_type) +
  # geom_text(aes(label = SNP_unique_table$kpi_SNP_count), 
  #           vjust = -0.25) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  ylim(c(0, 300)) +
  # ylim(c(0, max(SNP_unique_table$GREATenrichment))) +
  xlab("GREAT Enrichment of iPS-neural shared ASoC SNPs (FDR=0.1) \nin epigenetically annotated regions") +
  # xlab("GREAT Enrichment of cell-type-specific ASoC SNPs \nin epigenetically annotated regions") +
  ylab("-log10(Pvalue) fold enrichment")

####

SNP_unique_table <- as.data.frame(main_SNP_table[main_SNP_table$cell_type == "iPS_only_specific", 1:3])
SNP_unique_table[32, ] <- c("iPS_only_specific", "Enh_All", as.numeric(sum(SNP_unique_table$length[13:18])))
SNP_unique_table$length <- as.numeric(SNP_unique_table$length)

master_summary_table[32, ] <- c("Enh_All", sum(as.numeric(master_summary_table$total_feature_length[13:18])))
SNP_unique_table$total_feature_length <- master_summary_table$total_feature_length
# SNP_unique_table$total_feature_length <- SNP_length$total_SNP_length[3]
# SNP_length$total_SNP_length[3]
SNP_unique_table$GREATenrichment <- 1
i <- 1
# j <- 1
for (i in 1:nrow(SNP_unique_table)) {
  print(paste("i=", i, collapse = ""))
  print(paste(SNP_unique_table$annotation_type[i]))
  SNP_unique_table$GREATenrichment[i] <- (0 - log10(binom.test(x = as.numeric(SNP_unique_table$length[i]),
                                                               n = sum(SNP_unique_table$length), #357, #4411, #357# #SNP_unique_table$n_total_features_length[i],
                                                               p = as.numeric(as.character(SNP_unique_table$total_feature_length[i]))/3234830000,
                                                               alternative = "greater")$p.value))
  SNP_unique_table$p_GREAT[i] <- binom.test(x = as.numeric(SNP_unique_table$length[i]),
                                            n = sum(SNP_unique_table$length), #357, #4411, #357# #SNP_unique_table$n_total_features_length[i],
                                            p = as.numeric(as.character(SNP_unique_table$total_feature_length[i]))/3234830000,
                                            alternative = "greater")$p.value
}

p.adjust(SNP_unique_table$p_GREAT, method = "fdr")


ggplot(SNP_unique_table, aes(x = annotation_type, y = GREATenrichment)) +
  geom_col(fill = "turquoise4") +
  # scale_x_discrete(limits =
  #                    c("TssA", "PromU", "PromD1", "PromD2",
  #                      "Tx5", "Tx", "Tx3", "TxWk",
  #                      "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2",
  #                      "EnhAc", "DNase", "ZNF_Rpts",
  #                      "Het", "PromP", "PromBiv", "ReprPC", "Quies",
  #                      "forebrain_enh", "non_forebrain_enh")) +
  scale_x_discrete(limits = SNP_unique_table$annotation_type) +
  # geom_text(aes(label = SNP_unique_table$kpi_SNP_count), 
  #           vjust = -0.25) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  ylim(c(0, 300)) +
  # ylim(c(0, max(SNP_unique_table$GREATenrichment))) +
  xlab("GREAT Enrichment of iPS-specific ASoC SNPs \nin epigenetically annotated regions") +
  # xlab("GREAT Enrichment of cell-type-specific ASoC SNPs \nin epigenetically annotated regions") +
  ylab("-log10(Pvalue) fold enrichment")

####

SNP_unique_table <- as.data.frame(main_SNP_table[main_SNP_table$cell_type == "iPS_only_specific", 1:3])
SNP_unique_table[32, ] <- c("iPS_only_specific", "Enh_All", as.numeric(sum(SNP_unique_table$length[13:18])))
SNP_unique_table$length <- as.numeric(SNP_unique_table$length)

master_summary_table[32, ] <- c("Enh_All", sum(as.numeric(master_summary_table$total_feature_length[13:18])))
SNP_unique_table$total_feature_length <- master_summary_table$total_feature_length
# SNP_unique_table$total_feature_length <- SNP_length$total_SNP_length[3]
# SNP_length$total_SNP_length[3]
SNP_unique_table$GREATenrichment <- 1
i <- 1
# j <- 1
for (i in 1:nrow(SNP_unique_table)) {
  print(paste("i=", i, collapse = ""))
  print(paste(SNP_unique_table$annotation_type[i]))
  SNP_unique_table$GREATenrichment[i] <- (0 - log10(binom.test(x = as.numeric(SNP_unique_table$length[i]),
                                                               n = sum(SNP_unique_table$length), #357, #4411, #357# #SNP_unique_table$n_total_features_length[i],
                                                               p = as.numeric(as.character(SNP_unique_table$total_feature_length[i]))/3234830000,
                                                               alternative = "greater")$p.value))
  SNP_unique_table$p_GREAT[i] <- binom.test(x = as.numeric(SNP_unique_table$length[i]),
                                            n = sum(SNP_unique_table$length), #357, #4411, #357# #SNP_unique_table$n_total_features_length[i],
                                            p = as.numeric(as.character(SNP_unique_table$total_feature_length[i]))/3234830000,
                                            alternative = "greater")$p.value
}

p.adjust(SNP_unique_table$p_GREAT, method = "fdr")


ggplot(SNP_unique_table, aes(x = annotation_type, y = GREATenrichment)) +
  geom_col(fill = "turquoise4") +
  # scale_x_discrete(limits =
  #                    c("TssA", "PromU", "PromD1", "PromD2",
  #                      "Tx5", "Tx", "Tx3", "TxWk",
  #                      "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2",
  #                      "EnhAc", "DNase", "ZNF_Rpts",
  #                      "Het", "PromP", "PromBiv", "ReprPC", "Quies",
  #                      "forebrain_enh", "non_forebrain_enh")) +
  scale_x_discrete(limits = SNP_unique_table$annotation_type) +
  # geom_text(aes(label = SNP_unique_table$kpi_SNP_count), 
  #           vjust = -0.25) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  ylim(c(0, 300)) +
  # ylim(c(0, max(SNP_unique_table$GREATenrichment))) +
  xlab("GREAT Enrichment of iPS-specific ASoC SNPs \nin epigenetically annotated regions") +
  # xlab("GREAT Enrichment of cell-type-specific ASoC SNPs \nin epigenetically annotated regions") +
  ylab("-log10(Pvalue) fold enrichment")


####

main_SNP_table <- main_SNP_table_old_09May2019


SNP_unique_table <- as.data.frame(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 1:3])
SNP_unique_table[32, ] <- c("cell_type_specific", "Enh_All", as.numeric(sum(SNP_unique_table$length[4:9])))
SNP_unique_table$length <- as.numeric(SNP_unique_table$length)

master_summary_table[32, ] <- c("Enh_All", sum(as.numeric(master_summary_table$total_feature_length[4:9])))
SNP_unique_table$total_feature_length <- master_summary_table$total_feature_length
# SNP_unique_table$total_feature_length <- SNP_length$total_SNP_length[3]
# SNP_length$total_SNP_length[3]
SNP_unique_table$GREATenrichment <- 1
i <- 1
# j <- 1
for (i in 1:nrow(SNP_unique_table)) {
  print(paste("i=", i, as.numeric(as.character(SNP_unique_table$total_feature_length[i])), collapse = ""))
  print(paste("x=", SNP_unique_table$length[i], 
              "n=", sum(SNP_unique_table$length), 
              "p=", as.numeric(as.character(SNP_unique_table$total_feature_length[i]))/3234830000))
  # print(paste(SNP_unique_table$annotation_type[i]))
  SNP_unique_table$GREATenrichment[i] <- (0 - log10(binom.test(x = as.numeric(SNP_unique_table$length[i]),
                                                               n = sum(SNP_unique_table$length), #357, #4411, #357# #SNP_unique_table$n_total_features_length[i],
                                                               p = as.numeric(as.character(SNP_unique_table$total_feature_length[i]))/3234830000,
                                                               alternative = "greater")$p.value))
  SNP_unique_table$p_GREAT[i] <- binom.test(x = SNP_unique_table$length[i],
                                            n = sum(SNP_unique_table$length), #357, #4411, #357# #SNP_unique_table$n_total_features_length[i],
                                            p = as.numeric(as.character(SNP_unique_table$total_feature_length[i]))/3234830000,
                                            alternative = "greater")$p.value
}

p.adjust(SNP_unique_table$p_GREAT, method = "fdr")


ggplot(SNP_unique_table, aes(x = annotation_type, y = GREATenrichment)) +
  geom_col(fill = "turquoise4") +
  # scale_x_discrete(limits =
  #                    c("TssA", "PromU", "PromD1", "PromD2",
  #                      "Tx5", "Tx", "Tx3", "TxWk",
  #                      "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2",
  #                      "EnhAc", "DNase", "ZNF_Rpts",
  #                      "Het", "PromP", "PromBiv", "ReprPC", "Quies",
  #                      "forebrain_enh", "non_forebrain_enh")) +
  scale_x_discrete(limits = SNP_unique_table$annotation_type) +
  # geom_text(aes(label = SNP_unique_table$kpi_SNP_count), 
  #           vjust = -0.25) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  ylim(c(0, 300)) +
  # ylim(c(0, max(SNP_unique_table$GREATenrichment))) +
  xlab("GREAT Enrichment of iPS-specific ASoC SNPs \nin epigenetically annotated regions") +
  # xlab("GREAT Enrichment of cell-type-specific ASoC SNPs \nin epigenetically annotated regions") +
  ylab("-log10(Pvalue) fold enrichment")
