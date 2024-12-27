# Siwei 25 Jan 2022
# calculate ASoC ratio from the aggregation of 5 lines
# use direct count results

# init
library(ggplot2)
library(readr)

##### GABA #####
## make 0 hr group
df_0hr <- read_delim("direct_count_raw_import/GABA_fromGVCF__0hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_0hr <- read_delim("raw_import/GABA_fromGVCF__0hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

#colnames(df_0hr)[6:7] <- c("REF_C", "ALT_C")

df_0hr$REF_C <- apply(df_0hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))
df_0hr$ALT_C <- apply(df_0hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 2]))))

df_0hr <- df_0hr[df_0hr$REF_C > 1, ]
df_0hr <- df_0hr[df_0hr$ALT_C > 1, ]
df_0hr$DP <- df_0hr$REF_C + df_0hr$ALT_C

sum(df_0hr$DP > 19) # 10161

max(df_0hr$DP)
df_0hr <- df_0hr[df_0hr$DP > 19, ]


df_0hr$pVal <- apply(df_0hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[16]),
                                            n = sum(as.numeric(x[16]), 
                                                    as.numeric(x[17])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_0hr$FDR <- p.adjust(p = df_0hr$pVal,
                       method = "fdr")
sum(df_0hr$FDR < 0.05) # 432

ggplot(df_0hr,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle("GABA, 0hr, DP >= 20, 12866 SNPs,\nof 432 FDR < 0.05") +
  ylab("-log10P") +
  theme_classic()

sum((df_0hr$REF_C/df_0hr$DP) > 0.5)
sum((df_0hr$REF_C/df_0hr$DP) < 0.5)

#####
## make 0 hr group
df_0hr <- read_delim("direct_count_raw_import/CD_54_0hr_NEFMp_CUX2m_glut_Duan_018_2-0_ATAC_FL_S1_L004_WASPed_scARC_het.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

colnames(df_0hr)[6:7] <- c("REF_C", "ALT_C")

# df_0hr$REF_C <- apply(df_0hr, 1, 
#                       function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
# df_0hr$ALT_C <- apply(df_0hr, 1, 
#                       function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_0hr <- df_0hr[df_0hr$REF_C > 1, ]
df_0hr <- df_0hr[df_0hr$ALT_C > 1, ]
df_0hr$DP <- df_0hr$REF_C + df_0hr$ALT_C

sum(df_0hr$DP > 10) # 10161

max(df_0hr$DP)
df_0hr <- df_0hr[df_0hr$DP > 10, ]

df_0hr$pVal <- apply(df_0hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[6]),
                                            n = sum(as.numeric(x[6]), 
                                                    as.numeric(x[7])),
                                            p = 0.5,
                                            alternative = "t")$p.value))

# df_0hr$pVal <- apply(df_0hr, 1, 
#                      function(x)(binom.test(x = as.numeric(x[16]),
#                                             n = sum(as.numeric(x[16]), 
#                                                     as.numeric(x[17])),
#                                             p = 0.5,
#                                             alternative = "t")$p.value))
df_0hr$FDR <- p.adjust(p = df_0hr$pVal,
                       method = "fdr")
sum(df_0hr$FDR < 0.05) # 432

ggplot(df_0hr,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle("GABA, 0hr, DP >= 20, 12866 SNPs,\nof 432 FDR < 0.05") +
  ylab("-log10P") +
  theme_classic()

sum((df_0hr$REF_C/df_0hr$DP) > 0.5)
sum((df_0hr$REF_C/df_0hr$DP) < 0.5)
