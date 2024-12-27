# 23 Feb 2019
# Siwei 
# Use SNPs to make epigenetic 

CN_8 <- CN_8[, 1:3]
DN_8 <- DN_8[, 1:3]
NPC_8 <- NPC_8[, 1:3]
iPS_8 <- iPS_8[, 1:3]
GA_8 <- GA_8[, 1:3]

colnames(GA_8) <- colnames(CN_8)


SNP_lookup_table <- as.data.frame(rbind(CN_8, DN_8, GA_8, iPS_8, NPC_8))
SNP_lookup_table <- SNP_lookup_table[!duplicated(SNP_lookup_table$rsID), ]
SNP_lookup_table$END <- SNP_lookup_table$POS + 1
SNP_lookup_table <- SNP_lookup_table[, c(1,2,4,3)]

specific_coord <- SNP_lookup_table[SNP_lookup_table$rsID %in% specific_ASoC_all$X1, ]
shared_2_coord <- SNP_lookup_table[SNP_lookup_table$rsID %in% shared_ASoC_all$X1, ]
shared_3_coord <- SNP_lookup_table[SNP_lookup_table$rsID %in% at_least_3_types$x, ]

write.table(specific_coord, file = "cell_type_specific_SNP.narrowPeak", 
            quote = F, row.names = F, col.names = F, sep = "\t")
write.table(shared_2_coord, file = "shared_by_2_SNP.narrowPeak", 
            quote = F, row.names = F, col.names = F, sep = "\t")
write.table(shared_3_coord, file = "shared_by_3_SNP.narrowPeak", 
            quote = F, row.names = F, col.names = F, sep = "\t")


#########

colnames(main_SNP_table) <- c("cell_type", "annotation_type", "length")
colnames(SNP_length) <- c("cell_type", "total_SNP_length")


######
main_SNP_table <- main_SNP_table_old_09May2019
SNP_length <- SNP_length_old_09May2019
master_summary_table <- master_summary_table_old


SNP_unique_plot <- log10(as.numeric(c(
  (main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3] /
    (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000))$length,
  sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[4:9]) /
    (SNP_length$total_SNP_length[1]*sum(master_summary_table$total_feature_length[4:9])/3234830000)
)))



# SNP_unique_plot <- log10(SNP_unique_plot)

# master_summary_table[32, ] <- c("Enh A1+A2", sum(as.numeric(master_summary_table$total_feature_length[4:5])))

SNP_unique_bar <- data.frame(cbind(c(master_summary_table$annotation_type, "Enh_All"),
                                      SNP_unique_plot), 
                             stringsAsFactors = F)
# SNP_unique_bar[32, ] <- c("all_Enh_Added", sum())
colnames(SNP_unique_bar) <- c("annotation_type", "SNP_unique_plot")
SNP_unique_bar$order <- SNP_unique_bar$annotation_type


ggplot(SNP_unique_bar, aes(x = annotation_type, y = as.numeric(SNP_unique_bar$SNP_unique_plot))) +
  geom_col(fill = "royalblue3") +
  scale_x_discrete(limits =
                     c("TssA", "PromU", "PromD1", "PromD2",
                       "Tx5", "Tx", "Tx3", "TxWk",
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2",
                       "EnhAc", "DNase", "ZNF_Rpts",
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh", 
                       "H3K27acF", "H3K27acO", "H3K4me2F", "H3K4me2O", 
                       "Enh_All")) +
  geom_text(aes(label = as.vector(c(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length, 
                                    sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[4:9])))),
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
  xlab("Enrichment of cell-type-specific ASoC SNPs in epigenetically annotated regions") +
  ylab("log10(fold enrichment)")

######
# [1] "cell_type_specific" "shared_by_two"      "shared_by_three"   

shared_by_2_plot <- log10(
  main_SNP_table[main_SNP_table$cell_type == "shared_by_two", 3] /
    (SNP_length$total_SNP_length[2]*master_summary_table$total_feature_length/3234830000)
)

shared_by_2_bar <- as.data.frame(cbind(master_summary_table$annotation_type, shared_by_2_plot))
colnames(shared_by_2_bar) <- c("annotation_type", "shared_by_2_plot")

ggplot(shared_by_2_bar, aes(x = annotation_type, y = shared_by_2_plot)) +
  geom_col(fill = "salmon4") +
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
  xlab("Enrichment of ASoC SNPs shared by minimum 2 cell types \nin epigenetically annotated regions") +
  ylab("log10(fold enrichment)")

######
# [1] "cell_type_specific" "shared_by_two"      "shared_by_three"   

shared_by_3_plot <- log10(
  main_SNP_table[main_SNP_table$cell_type == "shared_by_three", 3] /
    (SNP_length$total_SNP_length[3]*master_summary_table$total_feature_length/3234830000)
)

shared_by_3_bar <- as.data.frame(cbind(master_summary_table$annotation_type, shared_by_3_plot))
colnames(shared_by_3_bar) <- c("annotation_type", "shared_by_3_plot")

ggplot(shared_by_3_bar, aes(x = annotation_type, y = shared_by_3_plot)) +
  geom_col(fill = "turquoise4") +
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
  xlab("Enrichment of ASoC SNPs shared by minimum 3 cell types \nin epigenetically annotated regions") +
  ylab("log10(fold enrichment)")

#############
# Enrichment
#############

SNP_unique_table <- as.data.frame(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 1:3])
# SNP_unique_table <- as.data.frame(cbind(master_summary_table$annotation_type, SNP_unique_table))
# colnames(SNP_unique_table) <- c("annotation_type", "SNP_unique_plot")
SNP_unique_table <- cbind(SNP_unique_table, master_summary_table)
SNP_unique_table$binomP <- 1
i <- 1
for(i in 1:nrow(SNP_unique_table)) {
  SNP_unique_table$binomP[i] <- 0-log10(binom.test(x = SNP_unique_table$SNP_unique_plot[i], 
                                                   n = SNP_unique_table$total_feature_length[i], 
                                                   p = SNP_unique_table$total_feature_length[i]/3234830000, 
                                                   alternative = "two.sided")$p.value)
}




SNP_unique_plot <- log10(
  main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3] /
    (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000)
)

SNP_unique_bar <- as.data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot))
colnames(SNP_unique_bar) <- c("annotation_type", "SNP_unique_plot")

ggplot(SNP_unique_bar, aes(x = annotation_type, y = SNP_unique_plot)) +
  geom_col(fill = "royalblue3") +
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
  xlab("Enrichment of cell-type-specific ASoC SNPs in epigenetically annotated regions") +
  ylab("log10(fold enrichment)")


######


SNP_unique_plot <- log10(as.numeric(c(
  (main_SNP_table[main_SNP_table$cell_type == "shared_by_three", 3] /
     (SNP_length$total_SNP_length[3]*master_summary_table$total_feature_length/3234830000))$length,
  sum(main_SNP_table[main_SNP_table$cell_type == "shared_by_three", 3]$length[4:9]) /
    (SNP_length$total_SNP_length[3]*sum(master_summary_table$total_feature_length[4:9])/3234830000)
)))



# SNP_unique_plot <- log10(SNP_unique_plot)

# master_summary_table[32, ] <- c("Enh A1+A2", sum(as.numeric(master_summary_table$total_feature_length[4:5])))

SNP_unique_bar <- data.frame(cbind(c(master_summary_table$annotation_type, "Enh_All"),
                                   SNP_unique_plot), 
                             stringsAsFactors = F)
# SNP_unique_bar[32, ] <- c("all_Enh_Added", sum())
colnames(SNP_unique_bar) <- c("annotation_type", "SNP_unique_plot")
SNP_unique_bar$order <- SNP_unique_bar$annotation_type


ggplot(SNP_unique_bar, aes(x = annotation_type, y = as.numeric(SNP_unique_bar$SNP_unique_plot))) +
  geom_col(fill = "turquoise4") +
  scale_x_discrete(limits =
                     c("TssA", "PromU", "PromD1", "PromD2",
                       "Tx5", "Tx", "Tx3", "TxWk",
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2",
                       "EnhAc", "DNase", "ZNF_Rpts",
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh", 
                       "H3K27acF", "H3K27acO", "H3K4me2F", "H3K4me2O", 
                       "Enh_All")) +
  theme_classic() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  xlab("Enrichment of ASoC SNPs shared by minimum 3 cell types \nin epigenetically annotated regions") +
  ylab("log10(fold enrichment)")

