# 14 May 2019
#############
# Enrichment
#############
# n is the total number of genomic regions in the input set;
# pπ is the a priori probability of selecting a base pair annotated with π
#     when selecting a single base pair uniformly from all non–assembly gap base pairs in the genome;
# kπ is the number of genomic regions in the input set that cause annotation π to be selected.
#############
# special library for very large number calculation required
library(Rmpfr)
library(Brobdingnag)
library(ggplot2)
library(ggpubr)
#############

# cell-type-specific
SNP_unique_table <- as.data.frame(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 1:3])
# SNPs shared by at least 3 SNPs
SNP_unique_table <- as.data.frame(main_SNP_table[main_SNP_table$cell_type == "shared_by_three", 1:3])

SNP_unique_table$annotation_type[26:27] <- c("forebrain_enh", "non_forebrain_enh")

SNP_unique_table <- as.data.frame(cbind(master_summary_table$total_feature_length, SNP_unique_table))
annotations_features_count$annotation_type[26:27] <- c("forebrain_enh", "non_forebrain_enh")
annotations_features_count$order[26:27] <- c("forebrain_enh", "non_forebrain_enh")
SNP_unique_table <- as.data.frame(merge(SNP_unique_table, annotations_features_count, 
                                        by = "annotation_type", all.x = T))
SNP_unique_table$order[8] <- "d_forebrain_enh"
SNP_unique_table$order[14] <- "d_non_forebrain_enh"
SNP_unique_table <- SNP_unique_table[order(SNP_unique_table$order), ]
# SNP_unique_table$enrichment <- 0
colnames(SNP_unique_table) <- c("annotation_type", "total_feature_length", "cell_type", 
                                "kpi_SNP_count", "n_total_features_length", "order")
# SNP_unique_table <- merge(SNP_unique_table, master_summary_table,
#                           by = "annotation_type")
SNP_unique_table$GREATenrichment <- 1
# SNP_unique_table$enrichment <- NULL
# binomP <- 1
i <- 1
# j <- 1
for(i in 1:nrow(SNP_unique_table)) {
  print(paste("i=", i, collapse = ""))
  print(paste(SNP_unique_table$annotation_type[i]))
  # print((chooseMpfr(n = SNP_unique_table$kpi_SNP_count[i], 
  #                   a = SNP_unique_table$n_total_features_length[i], 
  #                   p = SNP_unique_table$total_feature_length[i]/3234830000,
  #                   alternative = "greater")$p.value))
  # binomP <- mpfr(0, precBits = 3500)
  SNP_unique_table$GREATenrichment[i] <- (0-log10(binom.test(x = SNP_unique_table$kpi_SNP_count[i],
                                                          n = 4411, #SNP_unique_table$n_total_features_length[i],
                                                          p = SNP_unique_table$total_feature_length[i]/3234830000,
                                                          alternative = "greater")$p.value))
  # binomP <- mpfr(0, precBits = 3500)
  # for(j in SNP_unique_table$kpi_SNP_count[i]:SNP_unique_table$n_total_features_length[i]) {
  #   # print(paste("j=", j, collapse = ""))
  #   # binomP <- binomP + 
  #   #   chooseMpfr(a = SNP_unique_table$n_total_features_length[i], # use log to calculate approx value
  #   #                n = SNP_unique_table$kpi_SNP_count[i]) *
  #   #   mpfr((1 - 1/3234830000)^(SNP_unique_table$n_total_features_length - j), precBits = 3500) *
  #   #   mpfr((1/3234830000)^j, precBits = 3500) 
  #   # binomP <- binomP + pbinom()
  #   binomP <- binom.test(x = SNP_unique_table$kpi_SNP_count[])
  #   print(paste("binominal P=", binomP, collapse = ""))
  # }

  
}




# SNP_unique_plot <- log10(
#   main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3] /
#     (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000)
# # )
# 
# SNP_unique_bar <- as.data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot))
# colnames(SNP_unique_bar) <- c("annotation_type", "SNP_unique_plot")
ggplot(SNP_unique_table, aes(x = annotation_type, y = GREATenrichment)) +
  geom_col(fill = "royalblue3") +
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
  ylim(c(0, 250)) +
  # ylim(c(0, max(SNP_unique_table$GREATenrichment))) +
  xlab("GREAT Enrichment of cell-type-specific ASoC SNPs \nin epigenetically annotated regions") +
  ylab("-log10(Pvalue) fold enrichment")

###
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
  geom_text(aes(label = SNP_unique_table$kpi_SNP_count), 
            vjust = -0.25) +
  # ylim(c(0, 250)) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  # ylim(c(0, max(SNP_unique_table$GREATenrichment))) +
  xlab("GREAT Enrichment of ASoC SNPs shared by minimum 3 cell types \nin epigenetically annotated regions") +
  ylab("-log10(Pvalue) GREAT enrichment")

