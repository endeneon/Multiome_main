# Siwei 15 Nov 2022
# plot the published data on the allele balance of significant SNPs

library(readxl)
library(ggplot2)
library(readr)



excel_sheets("bulk_ATAC/aay3983_Suppl_Excel_seq1_v5.xlsx")

# [1] "Table summary"                   "S1_Annotated_elements"           "S2_20 subjects"                 
# [4] "S3_Next-gen sequencing stats"    "S4_peak overlap PsychENCODE"     "S5_iPSC_8_ASoC"                 
# [7] "S6_NPC_8_ASoC"                   "S7_Glut_8_ASoC"                  "S8_GA_8_ASoC"                   
# [10] "S9_DN_8_ASoC"                    "S10_ASoC variant Numbers"        "S11_SNP_NPC_20"                 
# [13] "S12_SNP_Glut_20_QTL_enrich"      "S13_iNGlut_ASoC_MPRA"            " S14_NPC_ASoC_MPRA"             
# [16] "S15_reporter gene assay in NPC"  "S16_ASoC_QTL_HiC"                "S17ASTB_vs.ASoC"                
# [19] "S18_Homer_motif_ASoC"            "S19_Homer_motif_OCRs"            "S20_ASoC_SZ_Glut20"             
# [22] "S21_ASoC_SZ_NPC20"               "S22_GWAS  for TORUS"             "S23_Multiplex _CROPSeq_gRNAs"   
# [25] "S24_cis_gene_DE_500kb"           "S25_primers_PCR assay"           "S26_SZ_loci_CROP_HiC"           
# [28] "S27_Summary_PIP_ATAC-seq_ASoC"   "S28_VPS45_editing_gRNA_primers"  "S29_BCL11B_editing_gRNA_primers"

df_to_plot <-
  read_excel("bulk_ATAC/aay3983_Suppl_Excel_seq1_v5.xlsx", 
             sheet = "S12_SNP_Glut_20_QTL_enrich", 
             skip = 1)

df_to_plot <-
  read_excel("bulk_ATAC/aay3983_Suppl_Excel_seq1_v5.xlsx", 
             sheet = "S5_iPSC_8_ASoC", 
             skip = 1)

df_to_plot <-
  read_excel("bulk_ATAC/Glut_NGN2_neuron_20_all_1ATAC-seq_ASoC_63066_FDRsig_8205.xlsx")

df_to_plot$REF_ratio <-
  df_to_plot$REF / (df_to_plot$REF + df_to_plot$ALT)
df_to_plot$`-log10FDR` <- 
  0 - log10(df_to_plot$FDR)

ggplot(df_to_plot,
       aes(x = REF_ratio,
           y = `-log10FDR`,
           colour = ifelse(test = FDR < 0.05,
                           yes = "FDR < 0.05",
                           no = "FDR > 0.05"))) +
  geom_point(size = 0.2) +
  scale_colour_manual(values = c("red", "black")) +
  labs(colour = "") +
  ylim(0, 50) +
  theme_classic() +
  ggtitle(label = paste("iPS 8 lines,\n N(FDR < 0.05) = ",
                        sum(df_to_plot$FDR < 0.05),
                        ",\n",
                        "N(FDR < 0.05 & REF/DP < 0.5) = ",
                        sum(df_to_plot$FDR < 0.05 & df_to_plot$REF_ratio < 0.5),
                        ",\n",
                        "N(FDR < 0.05 & REF/DP > 0.5) = ",
                        sum(df_to_plot$FDR < 0.05 & df_to_plot$REF_ratio > 0.5),
                        sep = ""))

sum(df_to_plot$FDR < 0.05)
sum(df_to_plot$FDR < 0.05 & df_to_plot$REF_ratio < 0.5)
sum(df_to_plot$FDR < 0.05 & df_to_plot$REF_ratio > 0.5)


######
df_to_plot <- iPS_8_REF
df_to_plot <- DN_core_8_REF_ALT_30Aug2018
df_to_plot <- new_GA_SNP_list


df_to_plot <-
  read_delim("/nvmefs/MG_17_lines_analysis/R_MG_17/DP20_data_files/RGlut_41_lines_500bp_SNP_262K_06Jun2022.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

df_to_plot <-
  read_delim("/nvmefs/MG_17_lines_analysis/R_MG_17/DP20_data_files/GABA_34_lines_SNP_307K_06Jun2022.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

df_to_plot <-
  read_delim("/nvmefs/MG_17_lines_analysis/R_MG_17/DP20_data_files/MG_28_lines_SNP_361KK_06Jun2022.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

df_to_plot <-
  read_delim("/nvmefs/MG_17_lines_analysis/R_MG_17/DP20_data_files/AST_12_lines_SNP_235K_06Jun2022.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

df_to_plot <-
  read_delim("/nvmefs/MG_17_lines_analysis/R_MG_17/DP20_data_files/outdated_files/DN_20_lines_SNP_174K_22Jun2021.txt", 
             delim = "\t", escape_double = FALSE, 
             trim_ws = TRUE)

df_to_plot$REF_ratio <-
  df_to_plot$REF_count / (df_to_plot$REF_count + df_to_plot$ALT_count)

df_to_plot$REF_ratio <-
  df_to_plot$REF_C / (df_to_plot$REF_C + df_to_plot$ALT_C)
df_to_plot$`-log10FDR` <- 
  0 - log10(df_to_plot$FDR)

df_to_plot$REF_ratio <-
  df_to_plot$REF_N / (df_to_plot$REF_N + df_to_plot$ALT_N)
df_to_plot$`-log10FDR` <- 
  0 - log10(df_to_plot$FDR)

ggplot(df_to_plot,
       aes(x = REF_ratio,
           y = `-log10FDR`,
           colour = ifelse(test = FDR < 0.05,
                           yes = "FDR < 0.05",
                           no = "FDR > 0.05"))) +
  geom_point(size = 0.2) +
  scale_colour_manual(values = c("red", "black")) +
  labs(colour = "") +
  ylim(0, 25) +
  theme_classic() +
  ggtitle(label = paste("GABA 8 lines,\nN(FDR < 0.05) = ",
                        sum(df_to_plot$FDR < 0.05),
                        ",\n",
                        "N(FDR < 0.05 & REF/DP < 0.5) = ",
                        sum(df_to_plot$FDR < 0.05 & df_to_plot$REF_ratio < 0.5),
                        ",\n",
                        "N(FDR < 0.05 & REF/DP > 0.5) = ",
                        sum(df_to_plot$FDR < 0.05 & df_to_plot$REF_ratio > 0.5),
                        sep = ""))
