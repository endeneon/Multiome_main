# Siwei 03 Jan 2022
# calculate ASoC ratio from the aggregation of 5 lines

# init
library(ggplot2)
library(readr)

## make a df_raw section for general use
df_raw <- read_delim("ASoC_18_line_raw_28Jul2022/GABA_scATAC_0hr_merged_SNP_28Jul2022_4_R.txt",
                     delim = "\t", escape_double = FALSE,
                     trim_ws = TRUE)
# 
# df_raw <- 
#   read_delim("ASoC_18_line_raw_28Jul2022/NEFM_neg_scATAC_0hr_merged_SNP_28Jul2022_4_R.txt", 
#              delim = "\t", escape_double = FALSE, 
#              trim_ws = TRUE)

df_raw <- read_delim("ASoC_18_line_raw_28Jul2022/NEFM_pos_scATAC_0hr_merged_SNP_28Jul2022_4_R.txt",
                     delim = "\t", escape_double = FALSE,
                     trim_ws = TRUE)

# df_raw <- read_delim("ASoC_18_line_raw_28Jul2022/GABA_scATAC_1hr_merged_SNP_28Jul2022_4_R.txt", 
#                      delim = "\t", escape_double = FALSE, 
#                      trim_ws = TRUE)

df_raw <- read_delim("ASoC_18_line_raw_28Jul2022/NEFM_neg_scATAC_1hr_merged_SNP_28Jul2022_4_R.txt",
                     delim = "\t", escape_double = FALSE,
                     trim_ws = TRUE)

# df_raw <- read_delim("ASoC_18_line_raw_28Jul2022/NEFM_pos_scATAC_1hr_merged_SNP_28Jul2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)

# df_raw <- read_delim("ASoC_18_line_raw_28Jul2022/NEFM_neg_scATAC_6hr_merged_SNP_28Jul2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)
# 
# df_raw <- read_delim("ASoC_18_line_raw_28Jul2022/NEFM_pos_scATAC_6hr_merged_SNP_28Jul2022_4_R.txt",
#                      delim = "\t", escape_double = FALSE,
#                      trim_ws = TRUE)

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

df_to_plot <- df_to_plot[df_to_plot$DP > 19, ]


df_to_plot$pVal <- apply(df_to_plot, 1, 
                     function(x)(binom.test(x = as.numeric(x[ncol(df_to_plot) - 1]),
                                            n = sum(as.numeric(x[ncol(df_to_plot) - 1]), 
                                                    as.numeric(x[ncol(df_to_plot) - 2])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_to_plot$FDR <- p.adjust(p = df_to_plot$pVal,
                       method = "fdr")

View(df_to_plot[df_to_plot$ID %in% "rs2027349", ])

## !! assign time and cell type !!
ggplot(df_to_plot,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle(paste("GABA, 0hr, DP >= 20,\n",
                nrow(df_to_plot),
                "SNPs, of",
                sum(df_to_plot$FDR < 0.05),
                "FDR < 0.05")) +
  ylab("-log10P") +
  theme_classic()




##### GABA #####
## make 0 hr group
df_0hr <- read_delim("ASoC_18_line_raw_28Jul2022/GABA_scATAC_0hr_merged_SNP_28Jul2022_4_R.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_0hr <- df_0hr[df_0hr$ID != ".", ]

df_0hr$REF_C <- 
  apply(df_0hr, 1, 
        function(x)(sum(as.numeric(x[seq(from = 7,
                                         to = ncol(df_0hr) - 1,
                                         by = 2)]))))
df_0hr$ALT_C <- 
  apply(df_0hr, 1, 
        function(x)(sum(as.numeric(x[seq(from = 8,
                                         to = ncol(df_0hr),
                                         by = 2)]))))


df_0hr <- df_0hr[df_0hr$REF_C > 1, ]
df_0hr <- df_0hr[df_0hr$ALT_C > 1, ]
df_0hr$DP <- df_0hr$REF_C + df_0hr$ALT_C

sum(df_0hr$DP > 19) # 10161

df_0hr <- df_0hr[df_0hr$DP > 19, ]


df_0hr$pVal <- apply(df_0hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[ncol(df_0hr) - 1]),
                                            n = sum(as.numeric(x[ncol(df_0hr) - 1]), 
                                                    as.numeric(x[ncol(df_0hr) - 2])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_0hr$FDR <- p.adjust(p = df_0hr$pVal,
                       method = "fdr")

sum((df_0hr$REF_C/df_0hr$DP) > 0.5)
sum((df_0hr$REF_C/df_0hr$DP) < 0.5)

sum(df_0hr$FDR < 0.05)

ggplot(df_0hr,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle("GABA, 0hr, DP >= 20, 93727 SNPs,\nof 3127 FDR < 0.05") +
  ylab("-log10P") +
  theme_classic()

## make 1 hr group
df_1hr <- read_delim("raw_import/GABA__1hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_1hr$REF_C <- apply(df_1hr, 1, 
                    function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_1hr$ALT_C <- apply(df_1hr, 1, 
                    function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_1hr <- df_1hr[df_1hr$REF_C > 1, ]
df_1hr <- df_1hr[df_1hr$ALT_C > 1, ]

df_1hr$DP <- df_1hr$REF_C + df_1hr$ALT_C

df_1hr <- df_1hr[df_1hr$DP > 19, ]


df_1hr$pVal <- apply(df_1hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[16]),
                                            n = sum(as.numeric(x[16]), 
                                                    as.numeric(x[17])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_1hr$FDR <- p.adjust(p = df_1hr$pVal,
                       method = "fdr")

## make 6 hr group
df_6hr <- read_delim("raw_import/GABA__6hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_6hr$REF_C <- apply(df_6hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_6hr$ALT_C <- apply(df_6hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_6hr <- df_6hr[df_6hr$REF_C > 1, ]
df_6hr <- df_6hr[df_6hr$ALT_C > 1, ]

df_6hr$DP <- df_6hr$REF_C + df_6hr$ALT_C

df_6hr <- df_6hr[df_6hr$DP > 19, ]


df_6hr$pVal <- apply(df_6hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[16]),
                                            n = sum(as.numeric(x[16]), 
                                                    as.numeric(x[17])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_6hr$FDR <- p.adjust(p = df_6hr$pVal,
                       method = "fdr")
df_6hr[df_6hr$ID %in% "rs11606396", ]

### see if any different ASoC SNPs at the three time points
sum(df_0hr$FDR < 0.05) # 212
sum(df_1hr$FDR < 0.05) # 60
sum(df_6hr$FDR < 0.05) # 352

median(df_0hr$DP)
median(df_1hr$DP)
median(df_6hr$DP)

Reduce(intersect, list(df_0hr$FDR[df_0hr$FDR < 0.05],
                       df_6hr$FDR[df_1hr$FDR < 0.05]))

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in%
      df_1hr$ID[df_1hr$FDR < 0.05]) # 2
df_0hr$ID[df_0hr$FDR < 0.05][(df_0hr$ID[df_0hr$FDR < 0.05] %in%
      df_1hr$ID[df_1hr$FDR < 0.05])] # 2


sum(df_0hr$ID[df_0hr$FDR < 0.05] %in%
      df_6hr$ID[df_6hr$FDR < 0.05]) # 10
df_0hr$ID[df_0hr$FDR < 0.05] [(df_0hr$ID[df_0hr$FDR < 0.05] %in%
      df_6hr$ID[df_6hr$FDR < 0.05])] # 10

sum(df_1hr$ID[df_1hr$FDR < 0.05] %in%
      df_6hr$ID[df_6hr$FDR < 0.05]) # 8
df_1hr$ID[df_1hr$FDR < 0.05][(df_1hr$ID[df_1hr$FDR < 0.05] %in%
      df_6hr$ID[df_6hr$FDR < 0.05])]

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in%
      df_1hr$ID[df_1hr$FDR < 0.05] %in%
      df_6hr$ID[df_6hr$FDR < 0.05]) # 


sum(df_0hr$ID[df_0hr$pVal < 0.01] %in%
      df_1hr$ID[df_1hr$pVal < 0.01]) # 9

df_1hr[df_1hr$ID %in% "rs1010837814", ]
View(df_6hr[df_6hr$ID %in% "rs1010837814", ])

ggplot(df_0hr,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle("GABA, 0hr, DP >= 20, 12866 SNPs,\nof 352 FDR < 0.05") +
  ylab("-log10P") +
  theme_classic()


temp_writeout_table <- df_0hr
temp_writeout_table <- df_1hr
temp_writeout_table <- df_6hr
output_table <- data.frame(cbind(temp_writeout_table$CHROM,
                                 temp_writeout_table$POS,
                                 temp_writeout_table$POS,
                                 temp_writeout_table$REF,
                                 temp_writeout_table$ALT,
                                 temp_writeout_table$ID,
                                 temp_writeout_table$REF_C,
                                 temp_writeout_table$ALT_C,
                                 temp_writeout_table$DP,
                                 temp_writeout_table$pVal,
                                 temp_writeout_table$FDR))
# output_table <- output_table[order(output_table$X8)]

write.table(output_table,
            file = "txt_4_annovar/GABA_0hr.avinput",
            quote = F,
            sep = "\t", 
            col.names = T, row.names = F)






##### NEFMm_CUXp #####
## make 0 hr group
df_0hr <- read_delim("raw_import/NEMFm__0hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_0hr$REF_C <- apply(df_0hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_0hr$ALT_C <- apply(df_0hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_0hr <- df_0hr[df_0hr$REF_C > 1, ]
df_0hr <- df_0hr[df_0hr$ALT_C > 1, ]
df_0hr$DP <- df_0hr$REF_C + df_0hr$ALT_C

sum(df_0hr$DP > 9) # 2086

df_0hr <- df_0hr[df_0hr$DP > 9, ]


df_0hr$pVal <- apply(df_0hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[17]),
                                            n = sum(as.numeric(x[17]), 
                                                    as.numeric(x[16])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_0hr$FDR <- p.adjust(p = df_0hr$pVal,
                       method = "fdr")
df_0hr <- df_0hr[order(df_0hr$pVal), ]

sum(df_0hr$FDR < 0.05) # 1

# df_0hr_FDR <- df_0hr[df_0hr$FDR < 0.05, ]

## make 1 hr group
df_1hr <- read_delim("raw_import/NEMFm__1hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)


df_1hr$REF_C <- apply(df_1hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_1hr$ALT_C <- apply(df_1hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_1hr <- df_1hr[df_1hr$REF_C > 1, ]
df_1hr <- df_1hr[df_1hr$ALT_C > 1, ]
df_1hr$DP <- df_1hr$REF_C + df_1hr$ALT_C

sum(df_1hr$DP > 9) # 2794

df_1hr <- df_1hr[df_1hr$DP > 9, ]


df_1hr$pVal <- apply(df_1hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[17]),
                                            n = sum(as.numeric(x[17]), 
                                                    as.numeric(x[16])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_1hr$FDR <- p.adjust(p = df_1hr$pVal,
                       method = "fdr")
df_1hr <- df_1hr[order(df_1hr$pVal), ]

sum(df_1hr$FDR < 0.05) # 4

## make 6 hr group
df_6hr <- read_delim("raw_import/NEMFm__6hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_6hr$REF_C <- apply(df_6hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_6hr$ALT_C <- apply(df_6hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_6hr <- df_6hr[df_6hr$REF_C > 1, ]
df_6hr <- df_6hr[df_6hr$ALT_C > 1, ]
df_6hr$DP <- df_6hr$REF_C + df_6hr$ALT_C

sum(df_6hr$DP > 9) # 7069

df_6hr <- df_6hr[df_6hr$DP > 9, ]


df_6hr$pVal <- apply(df_6hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[16]),
                                            n = sum(as.numeric(x[16]), 
                                                    as.numeric(x[17])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_6hr$FDR <- p.adjust(p = df_6hr$pVal,
                       method = "fdr")

### see if any different ASoC SNPs at the three time points
sum(df_0hr$FDR < 0.05) # 1
sum(df_1hr$FDR < 0.05) # 4
sum(df_6hr$FDR < 0.05) # 2

Reduce(intersect, list(df_0hr$FDR[df_0hr$FDR < 0.05],
                       df_6hr$FDR[df_1hr$FDR < 0.05]))

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in%
      df_1hr$ID[df_1hr$FDR < 0.05]) # 0

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in%
      df_6hr$ID[df_6hr$FDR < 0.05]) # 5

sum(df_1hr$ID[df_1hr$FDR < 0.05] %in%
      df_6hr$ID[df_6hr$FDR < 0.05]) # 3

sum(df_0hr$ID[df_0hr$pVal < 0.01] %in%
      df_1hr$ID[df_1hr$pVal < 0.01]) # 9

df_1hr[df_1hr$ID %in% "rs1010837814", ]
View(df_6hr[df_6hr$ID %in% "rs1010837814", ])

ggplot(df_6hr,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle("NEFMm_CUXp_Glut, 6hr, DP >= 10,\n1910 SNPs of 2 FDR < 0.05") +
  ylab("-log10P") +
  theme_classic()


##### NEFMp_CUXm #####
## make 0 hr group
df_0hr <- read_delim("raw_import/NEFMp__0hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_0hr$REF <- apply(df_0hr, 1, 
                    function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12)]))))
df_0hr$ALT <- apply(df_0hr, 1, 
                    function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 1]))))

df_0hr <- df_0hr[df_0hr$REF > 1, ]
df_0hr <- df_0hr[df_0hr$ALT > 1, ]
df_0hr$DP <- df_0hr$REF + df_0hr$ALT

df_0hr <- df_0hr[df_0hr$DP > 9, ]


df_0hr$pVal <- apply(df_0hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[14]),
                                            n = sum(as.numeric(x[14]), 
                                                    as.numeric(x[15])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_0hr$FDR <- p.adjust(p = df_0hr$pVal,
                       method = "fdr")

## make 1 hr group
df_1hr <- read_delim("raw_import/NEFMp__1hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_1hr$REF <- apply(df_1hr, 1, 
                    function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12)]))))
df_1hr$ALT <- apply(df_1hr, 1, 
                    function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 1]))))

df_1hr <- df_1hr[df_1hr$REF > 1, ]
df_1hr <- df_1hr[df_1hr$ALT > 1, ]

df_1hr$DP <- df_1hr$REF + df_1hr$ALT

hist(df_1hr$DP, breaks = 100)
df_1hr <- df_1hr[df_1hr$DP > 9, ]


df_1hr$pVal <- apply(df_1hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[14]),
                                            n = sum(as.numeric(x[14]), 
                                                    as.numeric(x[15])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_1hr$FDR <- p.adjust(p = df_1hr$pVal,
                       method = "fdr")

## make 6 hr group
df_6hr <- read_delim("raw_import/NEFMp__6hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_6hr$REF <- apply(df_6hr, 1, 
                    function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12)]))))
df_6hr$ALT <- apply(df_6hr, 1, 
                    function(x)(sum(as.numeric(x[c(4, 6, 8, 10, 12) + 1]))))

df_6hr <- df_6hr[df_6hr$REF > 1, ]
df_6hr <- df_6hr[df_6hr$ALT > 1, ]

df_6hr$DP <- df_6hr$REF + df_6hr$ALT

df_6hr <- df_6hr[df_6hr$DP > 9, ]

df_6hr$pVal <- apply(df_6hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[14]),
                                            n = sum(as.numeric(x[14]), 
                                                    as.numeric(x[15])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_6hr$FDR <- p.adjust(p = df_6hr$pVal,
                       method = "fdr")

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in% 
      df_1hr$ID[df_1hr$FDR < 0.05]) # 2

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in% 
      df_6hr$ID[df_6hr$FDR < 0.05]) # 10

sum(df_1hr$ID[df_1hr$FDR < 0.05] %in% 
      df_6hr$ID[df_6hr$FDR < 0.05]) # 8

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in% 
      df_1hr$ID[df_1hr$FDR < 0.05] %in%
      df_6hr$ID[df_6hr$FDR < 0.05]) # 0

### see if any different ASoC SNPs at the three time points
sum(df_0hr$FDR < 0.05) # 129
sum(df_1hr$FDR < 0.05) # 54
sum(df_6hr$FDR < 0.05) # 145

Reduce(intersect, list(df_0hr$FDR[df_0hr$FDR < 0.05],
                       df_6hr$FDR[df_1hr$FDR < 0.05]))

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in%
      df_1hr$ID[df_1hr$FDR < 0.05]) # 5

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in%
      df_6hr$ID[df_6hr$FDR < 0.05]) # 3

sum(df_1hr$ID[df_1hr$FDR < 0.05] %in%
      df_6hr$ID[df_6hr$FDR < 0.05]) # 2

sum(df_0hr$ID[df_0hr$pVal < 0.01] %in%
      df_1hr$ID[df_1hr$pVal < 0.01]) # 

df_1hr[df_1hr$ID %in% "rs1010837814", ]
View(df_6hr[df_6hr$ID %in% "rs1010837814", ])

ggplot(df_0hr,
       aes(x = (REF/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle("NEFMp_CUXm_Glut, 0hr,\nDP >= 10, 27523 SNPs") +
  ylab("-log10P") +
  theme_classic()


##### NPC #####
## make 0 hr group
df_0hr <- read_delim("raw_import/NPC__0hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_0hr$REF_C <- apply(df_0hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_0hr$ALT_C <- apply(df_0hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_0hr <- df_0hr[df_0hr$REF_C > 1, ]
df_0hr <- df_0hr[df_0hr$ALT_C > 1, ]
df_0hr$DP <- df_0hr$REF_C + df_0hr$ALT_C

sum(df_0hr$DP > 9) # 10161

df_0hr <- df_0hr[df_0hr$DP > 9, ]


df_0hr$pVal <- apply(df_0hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[16]),
                                            n = sum(as.numeric(x[16]), 
                                                    as.numeric(x[17])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_0hr$FDR <- p.adjust(p = df_0hr$pVal,
                       method = "fdr")
df_0hr <- df_0hr[order(df_0hr$pVal), ]

sum(df_0hr$FDR < 0.05) # 34

## make 1 hr group
df_1hr <- read_delim("raw_import/NPC__1hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_1hr$REF_C <- apply(df_1hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_1hr$ALT_C <- apply(df_1hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_1hr <- df_1hr[df_1hr$REF_C > 1, ]
df_1hr <- df_1hr[df_1hr$ALT_C > 1, ]
df_1hr$DP <- df_1hr$REF_C + df_1hr$ALT_C

sum(df_1hr$DP > 9) # 7069

df_1hr <- df_1hr[df_1hr$DP > 9, ]


df_1hr$pVal <- apply(df_1hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[16]),
                                            n = sum(as.numeric(x[16]), 
                                                    as.numeric(x[17])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_1hr$FDR <- p.adjust(p = df_1hr$pVal,
                       method = "fdr")
df_1hr <- df_1hr[order(df_1hr$pVal), ]

sum(df_1hr$FDR < 0.05) # 4


## make 6 hr group
df_6hr <- read_delim("raw_import/NPC__6hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_6hr$REF_C <- apply(df_6hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_6hr$ALT_C <- apply(df_6hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_6hr <- df_6hr[df_6hr$REF_C > 1, ]
df_6hr <- df_6hr[df_6hr$ALT_C > 1, ]
df_6hr$DP <- df_6hr$REF_C + df_6hr$ALT_C

sum(df_6hr$DP > 9) # 4866

df_6hr <- df_6hr[df_6hr$DP > 9, ]


df_6hr$pVal <- apply(df_6hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[16]),
                                            n = sum(as.numeric(x[16]), 
                                                    as.numeric(x[17])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_6hr$FDR <- p.adjust(p = df_6hr$pVal,
                       method = "fdr")
df_6hr <- df_6hr[order(df_6hr$pVal), ]

sum(df_6hr$FDR < 0.05) # 22

### see if any different ASoC SNPs at the three time points
sum(df_0hr$FDR < 0.05) # 34
sum(df_1hr$FDR < 0.05) # 4
sum(df_6hr$FDR < 0.05) # 22

max(df_0hr$DP)
max(df_1hr$DP)
max(df_6hr$DP)

median(df_0hr$DP)
median(df_1hr$DP)
median(df_6hr$DP)

Reduce(intersect, list(df_0hr$FDR[df_0hr$FDR < 0.05],
                       df_6hr$FDR[df_1hr$FDR < 0.05]))

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in%
      df_1hr$ID[df_1hr$FDR < 0.05]) # 0

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in%
      df_6hr$ID[df_6hr$FDR < 0.05]) # 0

sum(df_1hr$ID[df_1hr$FDR < 0.05] %in%
      df_6hr$ID[df_6hr$FDR < 0.05]) # 0

sum(df_0hr$ID[df_0hr$pVal < 0.01] %in%
      df_1hr$ID[df_1hr$pVal < 0.01]) # 3

sum(df_0hr$ID[df_0hr$pVal < 0.01] %in%
      df_6hr$ID[df_6hr$pVal < 0.01]) # 9

df_1hr[df_1hr$ID %in% "rs1010837814", ]
View(df_6hr[df_6hr$ID %in% "rs1010837814", ])

ggplot(df_6hr,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle("NPC, 6hr, DP >= 10, \n4866 SNPs of 22 FDR < 0.05") +
  ylab("-log10P") +
  theme_classic()


temp_writeout_table <- df_0hr
temp_writeout_table <- df_1hr
temp_writeout_table <- df_6hr
output_table <- data.frame(cbind(temp_writeout_table$CHROM,
                                 temp_writeout_table$POS,
                                 temp_writeout_table$POS,
                                 temp_writeout_table$REF,
                                 temp_writeout_table$ALT,
                                 temp_writeout_table$ID,
                                 temp_writeout_table$REF_C,
                                 temp_writeout_table$ALT_C,
                                 temp_writeout_table$DP,
                                 temp_writeout_table$pVal,
                                 temp_writeout_table$FDR))
# output_table <- output_table[order(output_table$X8)]

write.table(output_table,
            file = "txt_4_annovar/NPC_0hr.avinput",
            quote = F,
            sep = "\t", 
            col.names = F, row.names = F)
