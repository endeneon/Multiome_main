# Siwei 05 Aug 2022
# calculate ASoC ratio from the aggregation

# init
library(ggplot2)
library(ggrepel)
library(readr)

library(stringr)

library(RColorBrewer)

# df_raw <- read_delim("bwa_call_test/NEFM_pos_glut_scATAC_0hr_rededup_no_VQSR_MAPQ_30_05Aug2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)
# 
# df_raw <- read_delim("bwa_call_test/NEFM_pos_glut_scATAC_6hr_merged_10x_dedup_SNP_05Aug2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)
# 
# df_raw <- read_delim("bwa_call_test/NEFM_pos_glut_scATAC_0hr_no_VQSR_MAPQ_30_05Aug2022_4_R_5x.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)
# 
# df_raw <- read_delim("pre_WASP_test/GABA_0hr_scATAC_no_VQSR_MAPQ_30_pre_WASP_10Aug2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)
# 
# df_raw <- read_delim("pre_WASP_test/NEFM_neg_0hr_scATAC_no_VQSR_MAPQ_30_pre_WASP_10Aug2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)
# 
# df_raw <- read_delim("pre_WASP_test/NEFM_pos_0hr_scATAC_no_VQSR_MAPQ_30_pre_WASP_10Aug2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)
# 
# df_raw <- read_delim("bowtie2_test/GABA_hr_0_scATAC_no_VQSR_MAPQ_30_bowtie2_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)
# 
# df_raw <- read_delim("bowtie2_test/NEFM_neg_hr_0_scATAC_no_VQSR_MAPQ_30_bowtie2_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)

# df_raw <- read_delim("bowtie2_test/NEFM_pos_hr_0_scATAC_no_VQSR_MAPQ_30_bowtie2_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)

df_raw <- read_delim("NotDuplicateReadFilter_15Aug2022/GABA_0hr_scATAC_NotDuplicatedReadFilter_15Aug2022_4_R.txt",
                     delim = "\t", escape_double = FALSE,
                     trim_ws = TRUE)
df_raw <- read_delim("NotDuplicateReadFilter_15Aug2022/NEFM_neg_0hr_scATAC_NotDuplicatedReadFilter_15Aug2022_4_R.txt",
                     delim = "\t", escape_double = FALSE,
                     trim_ws = TRUE)
df_raw <- read_delim("NotDuplicateReadFilter_15Aug2022/NEFM_pos_0hr_scATAC_NotDuplicatedReadFilter_15Aug2022_4_R.txt",
                     delim = "\t", escape_double = FALSE,
                     trim_ws = TRUE)

# df_raw <- read_delim("NotDuplicateReadFilter_15Aug2022/GABA_1hr_scATAC_NotDuplicatedReadFilter_15Aug2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)
# df_raw <- read_delim("NotDuplicateReadFilter_15Aug2022/NEFM_neg_1hr_scATAC_NotDuplicatedReadFilter_15Aug2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)
# df_raw <- read_delim("NotDuplicateReadFilter_15Aug2022/NEFM_pos_1hr_scATAC_NotDuplicatedReadFilter_15Aug2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)

# df_raw <- read_delim("NotDuplicateReadFilter_15Aug2022/GABA_6hr_scATAC_NotDuplicatedReadFilter_15Aug2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)
# df_raw <- read_delim("NotDuplicateReadFilter_15Aug2022/NEFM_neg_6hr_scATAC_NotDuplicatedReadFilter_15Aug2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)
df_raw <- read_delim("NotDuplicateReadFilter_15Aug2022/NEFM_pos_6hr_scATAC_NotDuplicatedReadFilter_15Aug2022_4_R.txt",
                     delim = "\t", escape_double = FALSE,
                     trim_ws = TRUE)

df_raw <- read_delim("0hr_CD_63_GABA_WASP_and_dedup_all_by_barcode_23Sept2022.txt",
                     delim = "\t", escape_double = FALSE,
                     trim_ws = TRUE)


df_to_plot <- df_raw
df_to_plot <- df_to_plot[df_to_plot$ID != ".", ]

df_to_plot$REF_C <- 
  apply(df_to_plot, 1, 
        function(x)(sum(as.numeric(x[seq(from = 7,
                                         to = ncol(df_to_plot) - 1,
                                         by = 2)]))))
df_to_plot$ALT_C <- 
  apply(df_to_plot, 1, 
        function(x)(sum(as.numeric(x[seq(from = 8,
                                         to = ncol(df_to_plot),
                                         by = 2)]))))


df_to_plot <- df_to_plot[df_to_plot$REF_C > 1, ]
df_to_plot <- df_to_plot[df_to_plot$ALT_C > 1, ]
df_to_plot$DP <- df_to_plot$REF_C + df_to_plot$ALT_C

# sum(df_to_plot$DP > 19) 

df_to_plot <- df_to_plot[df_to_plot$DP > 8, ]


df_to_plot$pVal <- apply(df_to_plot, 1, 
                         function(x)(binom.test(x = as.numeric(x[ncol(df_to_plot) - 1]),
                                                n = sum(as.numeric(x[ncol(df_to_plot) - 1]), 
                                                        as.numeric(x[ncol(df_to_plot) - 2])),
                                                p = 0.5,
                                                alternative = "t")$p.value))
df_to_plot$FDR <- p.adjust(p = df_to_plot$pVal,
                           method = "fdr")

# df_writeout <- df_to_plot[c(1:6, 43:47)]

# View(df_to_plot[df_to_plot$ID %in% "rs2027349", ])

# df_to_plot[df_to_plot$ID %in% "rs2027349", ]$
# 
# df_to_plot$SNP <- NA
# df_to_plot$SNP[df_to_plot$ID %in% "rs2027349"] <- "rs2027349"


## !! assign time and cell type !!
ggplot(df_to_plot,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)),
           label = "SNP")) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle(paste("hr_6, NEFM_pos_glut, DP >= 20, 18 lines\n",
                # "Duplicates from picard MarkDuplicates\n",
                # "Duplicates marked by 10x Cell Ranger\n",
                nrow(df_to_plot),
                "SNPs, of",
                sum(df_to_plot$FDR < 0.05),
                "FDR < 0.05\n",
                sep = " ")) +
  ylab("-log10P") +
  theme_classic() # +
# geom_text_repel(max.time = 0.5, 
#                 na.rm = T,
#                 verbose = T)
## !! assign time and cell type !!
ggplot(df_to_plot,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)),
           label = "SNP")) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle(paste("hr_6, NEFM_pos_glut, DP >= 20, 18 lines\n",
                # "Duplicates from picard MarkDuplicates\n",
                # "Duplicates marked by 10x Cell Ranger\n",
                nrow(df_to_plot),
                "SNPs, of",
                sum(df_to_plot$FDR < 0.05),
                "FDR < 0.05\n",
                "rs2027349: REF_C =",
                df_to_plot[df_to_plot$ID %in% "rs2027349", ]$REF_C,
                ", DP = ",
                df_to_plot[df_to_plot$ID %in% "rs2027349", ]$DP, 
                ",\n Pval =",
                format(signif(df_to_plot[df_to_plot$ID %in% "rs2027349", ]$pVal, digits = 2),
                       nsmall = 2), 
                ", FDR =",
                format(signif(df_to_plot[df_to_plot$ID %in% "rs2027349", ]$FDR, digits = 2),
                       nsmall = 2),
                sep = " ")) +
  ylab("-log10P") +
  theme_classic() # +
  # geom_text_repel(max.time = 0.5, 
  #                 na.rm = T,
  #                 verbose = T)



sum(df_to_plot$REF_C/df_to_plot$DP < 0.5) #29685
sum(df_to_plot$REF_C/df_to_plot$DP > 0.5) #26342


sum((df_to_plot$REF_C/df_to_plot$DP < 0.5) & df_to_plot$FDR < 0.05)
sum((df_to_plot$REF_C/df_to_plot$DP > 0.5) & df_to_plot$FDR < 0.05)

##### make a loop to send output

input_file_list <-
  list.files(path = "NotDuplicateReadFilter_15Aug2022",
             pattern = "*15Aug2022_4_R.txt$", 
             full.names = T,
             include.dirs = F)

i <- 1

for (i in 1:length(input_file_list)) {
  file_name_output <-
    str_replace_all(string = input_file_list[i],
                    pattern = "_4_R\\.txt$",
                    replacement = "_results\\.tsv")
  print("calc")
  
  ### read data
  df_raw <- read_delim(input_file_list[i],
                       delim = "\t", escape_double = FALSE,
                       trim_ws = TRUE)
  
  df_to_plot <- df_raw
  df_to_plot <- df_to_plot[df_to_plot$ID != ".", ]
  
  df_to_plot$REF_C <- 
    apply(df_to_plot, 1, 
          function(x)(sum(as.numeric(x[seq(from = 7,
                                           to = ncol(df_to_plot) - 1,
                                           by = 2)]))))
  df_to_plot$ALT_C <- 
    apply(df_to_plot, 1, 
          function(x)(sum(as.numeric(x[seq(from = 8,
                                           to = ncol(df_to_plot),
                                           by = 2)]))))
  
  
  df_to_plot <- df_to_plot[df_to_plot$REF_C > 1, ]
  df_to_plot <- df_to_plot[df_to_plot$ALT_C > 1, ]
  df_to_plot$DP <- df_to_plot$REF_C + df_to_plot$ALT_C

  df_to_plot <- df_to_plot[df_to_plot$DP > 19, ]
  
  df_to_plot$pVal <- apply(df_to_plot, 1, 
                           function(x)(binom.test(x = as.numeric(x[ncol(df_to_plot) - 1]),
                                                  n = sum(as.numeric(x[ncol(df_to_plot) - 1]), 
                                                          as.numeric(x[ncol(df_to_plot) - 2])),
                                                  p = 0.5,
                                                  alternative = "t")$p.value))
  df_to_plot$FDR <- p.adjust(p = df_to_plot$pVal,
                             method = "fdr")
  
  df_writeout <- df_to_plot[c(1:6, 43:47)]
  
  print(paste(i, 
              ncol(df_writeout),
              sum(df_writeout$FDR < 0.05),
              sep = ","))
  
  write.table(df_writeout,
              file = file_name_output,
              quote = F, sep = "\t",
              row.names = F, col.names = T)
}
