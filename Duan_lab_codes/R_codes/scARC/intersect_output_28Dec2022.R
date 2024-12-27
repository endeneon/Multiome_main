# Siwei 28 Dec 2022
# Intersect 03 Oct results with mgs list

# init
library(readr)
library(vcfR)
library(stringr)
library(ggplot2)

library(parallel)

library(MASS)

library(RColorBrewer)
library(grDevices)

# outputNEFM_neg_0hr <- 
#   read_delim("BQSR_no_interval_23Sept2022/outputNEFM_neg_0hr.avinput", 
#              delim = "\t", escape_double = FALSE, 
#              col_names = FALSE, trim_ws = TRUE)
# colnames(outputNEFM_neg_0hr) <-
#   c("CHR", "START", "END", "REF", "ALT", "rsID",
#     "NCALLED", "REF_C", "ALT_C", "DP", "PVal", "FDR")

outputNEFM_neg_0hr <-
  read.vcfR(file = "/home/zhangs3/NVME/scARC_Duan_024_raw_fastq/GRCh38_output/hr_0/known_genotypes_neg_pt1/vcf_output_all_mgs/output/NEFM_neg_0hr_merged_SNP_VF.vcf",
            limit = 5e+8,
            verbose = T)


mgs_all_intervals <-
  read.vcfR(file = "~/Data/Databases/GWAS/mgs_new_imputation_guo/output_vcf/Duan_Project_024_samples.vcf",
            limit = 1e+8,
            verbose = T)

mgs_SNP_list <-
  as.data.frame(mgs_all_intervals@fix)
mgs_SNP_list$rsID <-
  str_split(string = mgs_SNP_list$ID,
            pattern = "\\:",
            simplify = T)[, 1]
mgs_SNP_list <-
  mgs_SNP_list[str_detect(string = mgs_SNP_list$rsID,
                          pattern = "^rs.*"), 
               ]

df_outputNEFM_neg_0hr <-
  as.data.frame(outputNEFM_neg_0hr@fix)
df_outputNEFM_neg_0hr$uuid <-
  str_c(df_outputNEFM_neg_0hr$CHROM,
        df_outputNEFM_neg_0hr$POS,
        sep = "_")

df_gt_df <-
  extract.gt(outputNEFM_neg_0hr,
             element = "DF",
             as.numeric = T)
df_gt_sum <-
  rowSums(x = df_gt_df,
          na.rm = T)

df_gt_dr <-
  extract.gt(outputNEFM_neg_0hr,
             element = "DR",
             as.numeric = T)
df_gt_dr_sum <-
  rowSums(x = df_gt_dr,
          na.rm = T)

df_count_outputNEFM_neg_0hr <-
  data.frame(uuid = df_outputNEFM_neg_0hr$uuid,
             # REF = df_outputNEFM_neg_0hr$REF,
             DF = df_gt_sum,
             DR = df_gt_dr_sum)

mgs_SNP_list$uuid <-
  str_c(mgs_SNP_list$CHROM,
        mgs_SNP_list$POS,
        sep = "_")

df_count_outputNEFM_neg_0hr <-
  merge(x = df_count_outputNEFM_neg_0hr,
        y = mgs_SNP_list,
        by = "uuid")
df_count_outputNEFM_neg_0hr$DP <-
  df_count_outputNEFM_neg_0hr$DF +
  df_count_outputNEFM_neg_0hr$DR



df_to_plot <- df_count_outputNEFM_neg_0hr
df_to_plot <- 
  df_to_plot[df_to_plot$DP > 19, ]

df_to_plot <- df_to_plot[df_to_plot$DP > 19, ]
df_to_plot <- 
  df_to_plot[((df_to_plot$DF > 1) &
                df_to_plot$DR > 1), ]
# df_to_plot <- df_to_plot[df_to_plot$DP < 1000, ]

clust_df_plot_use <- makeCluster(type = "FORK", 6)
clusterExport(clust_df_plot_use, "df_to_plot")
df_to_plot$pVal <- 
  parApply(cl = clust_df_plot_use, X = df_to_plot, MARGIN = 1,
           FUN = function(x)(binom.test(x = as.numeric(x[2]),
                                        n = as.numeric(x[ncol(df_to_plot)]),
                                        p = 0.5,
                                        alternative = "t")$p.value))
stopCluster(clust_df_plot_use)
rm(clust_df_plot_use)

df_to_plot$FDR <- 
  p.adjust(p = df_to_plot$pVal,
           method = "fdr")

df_to_plot <-
  df_to_plot[order(df_to_plot$pVal), ]

ggplot(df_to_plot,
       aes(x = (DF/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle(paste(#"",
    "bwa_align_direct_count_wo_HaplotypeCaller\n",
    "NEFM_neg_0hr, DP >= 20, minAllele >= 2, 18 lines\n",
    # "SNP lists from 23 Sept results, removed all NCALLED >= 14\n",
    nrow(df_to_plot),
    "SNPs, of",
    sum(df_to_plot$FDR < 0.05),
    "FDR < 0.05\n",
    "FDR < 0.05 & REF/DP < 0.5 = ", 
    sum((df_to_plot$DF/df_to_plot$DP < 0.5) & (df_to_plot$FDR < 0.05)),
    ";\n FDR < 0.05 & REF/DP > 0.5 = ", "",
    sum((df_to_plot$DF/df_to_plot$DP > 0.5) & (df_to_plot$FDR < 0.05)), "\n",
    "rs2027349: REF_C =",
    df_to_plot[df_to_plot$rsID %in% "rs2027349", ]$DF,
    ", DP = ",
    df_to_plot[df_to_plot$rsID %in% "rs2027349", ]$DP, 
    ",\n Pval =",
    format(signif(df_to_plot[df_to_plot$rsID %in% "rs2027349", ]$pVal, digits = 2),
           nsmall = 2), 
    ", FDR =",
    format(signif(df_to_plot[df_to_plot$rsID %in% "rs2027349", ]$FDR, digits = 2),
           nsmall = 2))) +
  ylab("-log10P") +
  ylim(0, 75) +
  theme_classic()

sum(df_to_plot$DF/df_to_plot$DP < 0.5)
sum(df_to_plot$DF/df_to_plot$DP > 0.5)

write.table(df_to_plot,
            file = "bwa_align_direct_count_wo_HaplotypeCaller_29Dec2022.txt",
            quote = F, sep = "\t",
            row.names = F, col.names = T)

View(df_to_plot[df_to_plot$rsID %in% "rs472498", ])
View(df_to_plot[df_to_plot$rsID %in% "rs12406605", ])


### generate BED file for homer enrichment
# deduplicate df_to_plot by uuid
df_to_plot <-
  df_to_plot[!duplicated(df_to_plot$uuid), ]

df_to_bed <-
  df_to_plot[, c("CHROM", "POS", "rsID", "FDR")]

df_to_bed$strand_info <-
  df_to_plot$DF/df_to_plot$DP
df_to_bed$POS <-
  as.numeric(df_to_bed$POS)

### case 1, DF/DP < 0.5 -> (+) strand
df_to_bed_enrich_final <-
  df_to_bed[df_to_bed$FDR < 0.05, ]
df_to_bed_enrich_final$strand <-
  ifelse(df_to_bed_enrich_final$strand_info < 0.5, 
         yes = "+",
         no = "-")

df_to_bed_enrich_final$`POS-1` <-
  df_to_bed_enrich_final$POS - 1

df_to_bed_enrich_final <-
  df_to_bed_enrich_final[, c(1, 7, 2, 3, 4, 6)]


df_to_bed_background_final <-
  df_to_bed[df_to_bed$FDR > 0.05, ]
df_to_bed_background_final <-
  df_to_bed_background_final[sample(x = nrow(df_to_bed_background_final),
                                    size = nrow(df_to_bed_enrich_final),
                                    replace = F),
                             ]
df_to_bed_background_final$strand <-
  ifelse(df_to_bed_background_final$strand_info < 0.5, 
         yes = "+",
         no = "-")

df_to_bed_background_final$`POS-1` <-
  df_to_bed_background_final$POS - 1

df_to_bed_background_final <-
  df_to_bed_background_final[, c(1, 7, 2, 3, 4, 6)]

write.table(df_to_bed_enrich_final,
            file = "bed_4_homer/bed_homer_enrich_DF_DP_less_005_plus_strand.bed",
            quote = F, sep = "\t", 
            row.names = F, col.names = F)
write.table(df_to_bed_background_final,
            file = "bed_4_homer/bed_homer_bkgrnd_DF_DP_less_005_plus_strand.bed",
            quote = F, sep = "\t", 
            row.names = F, col.names = F)


### case 2, DF/DP < 0.5 -> (+) strand
df_to_bed_enrich_final <-
  df_to_bed[df_to_bed$FDR < 0.05, ]
df_to_bed_enrich_final$strand <-
  ifelse(df_to_bed_enrich_final$strand_info > 0.5, 
         yes = "+",
         no = "-")

df_to_bed_enrich_final$`POS-1` <-
  df_to_bed_enrich_final$POS - 1

df_to_bed_enrich_final <-
  df_to_bed_enrich_final[, c(1, 7, 2, 3, 4, 6)]


df_to_bed_background_final <-
  df_to_bed[df_to_bed$FDR > 0.05, ]
df_to_bed_background_final <-
  df_to_bed_background_final[sample(x = nrow(df_to_bed_background_final),
                                    size = nrow(df_to_bed_enrich_final),
                                    replace = F),
  ]
df_to_bed_background_final$strand <-
  ifelse(df_to_bed_background_final$strand_info > 0.5, 
         yes = "+",
         no = "-")

df_to_bed_background_final$`POS-1` <-
  df_to_bed_background_final$POS - 1

df_to_bed_background_final <-
  df_to_bed_background_final[, c(1, 7, 2, 3, 4, 6)]

write.table(df_to_bed_enrich_final,
            file = "bed_4_homer/bed_homer_enrich_DF_DP_more_005_plus_strand.bed",
            quote = F, sep = "\t", 
            row.names = F, col.names = F)
write.table(df_to_bed_background_final,
            file = "bed_4_homer/bed_homer_bkgrnd_DF_DP_more_005_plus_strand.bed",
            quote = F, sep = "\t", 
            row.names = F, col.names = F)
