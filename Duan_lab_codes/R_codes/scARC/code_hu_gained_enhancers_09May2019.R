# 09 May 2019
# Siwei 

library(ggplot2)
library(ggpubr)

#########

main_SNP_table_25state <- main_SNP_table
master_summary_table_25state <- master_summary_table

# main_SNP_table <- rbind(main_SNP_table_25state, main_SNP_table)
# master_summary_table <- rbind(master_summary_table_25state, master_summary_table)

main_SNP_table[82:93, 1] <- c(rep("cell_type_specific", 4), 
                              rep("shared_by_two", 4), 
                              rep("shared_by_three", 4))
main_SNP_table[82:93, 2] <- rep(c("H3K27acF", "H3K4me2F", "H3K27acO", "H3K4me2O"), 3)

master_summary_table[28:31, 1] <- c("H3K27acF", "H3K4me2F", "H3K27acO", "H3K4me2O")

colnames(main_SNP_table) <- c("cell_type", "annotation_type", "length")
colnames(SNP_length) <- c("cell_type", "total_SNP_length")
colnames(master_summary_table) <- colnames(master_summary_table_25state)
########



######
SNP_unique_plot <- log10(
  main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3] /
    (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000)
)$length

SNP_unique_plot[32] <- 
  log10(sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[4:9]) /
          (SNP_length$total_SNP_length[1] * 
             sum(master_summary_table$total_feature_length[4:9])/3234830000))#, 


SNP_unique_bar <- data.frame(cbind(master_summary_table$annotation_type, SNP_unique_plot), stringsAsFactors = F)
# SNP_unique_bar$`master_summary_table$annotation_type` <- as.character(as.vector(SNP_unique_bar$`master_summary_table$annotation_type`))
SNP_unique_bar[32, ] <- c("Enh_All", 
                          log10(sum(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length[4:9]) /
                                  (SNP_length$total_SNP_length[1] * 
                                     sum(master_summary_table$total_feature_length[4:9])/3234830000))#, 
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
   geom_text(aes(label = as.vector(c(main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3]$length, "896"))),
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


SNP_unique_plot <- log10(
  main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3] /
    (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000)
)$length

shared_by_3_plot <- 
  main_SNP_table[main_SNP_table$cell_type == "shared_by_three", 3] /
    (SNP_length$total_SNP_length[3]*master_summary_table$total_feature_length/3234830000)

k <- 
  sum(shared_by_3_plot$length[4:9]) /
    (SNP_length$total_SNP_length[3]*sum(master_summary_table$total_feature_length[4:9])/3234830000)

shared_by_3_plot <- as.data.frame(log10(c(as.vector(shared_by_3_plot$length), k)))


# shared_by_3_plot$length[32] <- k

shared_by_3_bar <- as.data.frame(cbind(c(master_summary_table$annotation_type, "Enh_All"), shared_by_3_plot))
colnames(shared_by_3_bar) <- c("annotation_type", "shared_by_3_plot")

ggplot(shared_by_3_bar, aes(x = annotation_type, y = shared_by_3_plot)) +
  geom_col(fill = "turquoise4") +
  scale_x_discrete(limits = 
                     c("TssA", "PromU", "PromD1", "PromD2", 
                       "Tx5", "Tx", "Tx3", "TxWk", 
                       "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", 
                       "EnhAc", "DNase", "ZNF_Rpts", 
                       "Het", "PromP", "PromBiv", "ReprPC", "Quies",
                       "forebrain_enh", "non_forebrain_enh",
                       "H3K27acF", "H3K4me2F", "H3K27acO", "H3K4me2O", 
                       "Enh_All")) +
  geom_text(aes(label = c(main_SNP_table[main_SNP_table$cell_type == "shared_by_three", 3]$length, 
                          sum(main_SNP_table[main_SNP_table$cell_type == "shared_by_three", 3]$length[4:9]))), 
            # hjust = -0.25, 
            # y = 2,
            # y = (log10(
            #   main_SNP_table[main_SNP_table$cell_type == "cell_type_specific", 3] /
            #     (SNP_length$total_SNP_length[1]*master_summary_table$total_feature_length/3234830000)))$length 
            # + 0.1,
            # position = position_stack((vjust = 1)),
            hjust = "inward",
            angle = 90) +
  theme_classic() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  ggpubr::rotate_x_text() +
  xlab("Enrichment of ASoC SNPs shared by minimum 3 cell types \nin epigenetically annotated regions") +
  ylab("log10(fold enrichment)") #+
  ylim(c(-1, 2.5))
