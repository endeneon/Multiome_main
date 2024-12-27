# Siwei 19 Jan 2022
# calculate ASoC ratio from the aggregation of 5 lines
# Use Glut from NEFMp

# init
library(ggplot2)
library(readr)

##### Glut, NEFMp #####
## make 0 hr group
df_0hr <- read_delim("raw_import/NEFMp__0hr_new.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_0hr <- read_delim("raw_import/NEFMp_CUX2m_glut_fromGVCF__0hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_0hr$REF_C <- apply(df_0hr, 1, 
                    function(x)(sum(as.numeric(x[c(7, 9, 11, 13, 15)]))))
df_0hr$ALT_C <- apply(df_0hr, 1, 
                    function(x)(sum(as.numeric(x[c(8, 10, 12, 14, 16)]))))

df_0hr <- df_0hr[df_0hr$REF_C > 1, ]
df_0hr <- df_0hr[df_0hr$ALT_C > 1, ]
df_0hr$DP_GF <- df_0hr$REF_C + df_0hr$ALT_C

sum(df_0hr$DP > 19) # 10161
sum(df_0hr$DP_GF > 19) 
df_0hr <- df_0hr[df_0hr$DP > 19, ]
max(df_0hr$DP_GF)

df_0hr$pVal <- apply(df_0hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[17]),
                                            n = sum(as.numeric(x[17]), 
                                                    as.numeric(x[18])),
                                            p = 0.5,
                                            alternative = "t")$p.value))

df_0hr$pVal <- apply(df_0hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[17]),
                                            n = as.numeric(x[19]),
                                            p = 0.5,
                                            alternative = "t")$p.value))

df_0hr$FDR <- p.adjust(p = df_0hr$pVal,
                       method = "fdr")
df_0hr <- df_0hr[order(df_0hr$pVal), ]

sum(df_0hr$FDR < 0.05) # 247

# df_0hr_FDR <- df_0hr[df_0hr$FDR < 0.05, ]
ggplot(df_0hr,
       aes(x = (ALT_C/DP_GF),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  # ggtitle("NEFMp_CUXm_Glut, 6hr, DP >= 20,\n6921 SNPs of 267 FDR < 0.05") +
  ylab("-log10P") +
  theme_classic()

## make 1 hr group
df_1hr <- read_delim("raw_import/NEFMp__1hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_1hr$REF_C <- apply(df_1hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_1hr$ALT_C <- apply(df_1hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_1hr <- df_1hr[df_1hr$REF_C > 1, ]
df_1hr <- df_1hr[df_1hr$ALT_C > 1, ]
df_1hr$DP <- df_1hr$REF_C + df_1hr$ALT_C

sum(df_1hr$DP > 19) # 7069

df_1hr <- df_1hr[df_1hr$DP > 19, ]


df_1hr$pVal <- apply(df_1hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[16]),
                                            n = sum(as.numeric(x[16]), 
                                                    as.numeric(x[17])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_1hr$FDR <- p.adjust(p = df_1hr$pVal,
                       method = "fdr")
df_1hr <- df_1hr[order(df_1hr$pVal), ]

sum(df_1hr$FDR < 0.05) # 145

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in% df_1hr$ID[df_1hr$FDR < 0.05])
df_0hr$ID[df_0hr$FDR < 0.05][df_0hr$ID[df_0hr$FDR < 0.05] %in% 
                               df_1hr$ID[df_1hr$FDR < 0.05]]



## make 6 hr group
df_6hr <- read_delim("raw_import/NEFMp__6hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_6hr$REF_C <- apply(df_6hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_6hr$ALT_C <- apply(df_6hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_6hr <- df_6hr[df_6hr$REF_C > 1, ]
df_6hr <- df_6hr[df_6hr$ALT_C > 1, ]
df_6hr$DP <- df_6hr$REF_C + df_6hr$ALT_C

sum(df_6hr$DP > 19) # 7069

df_6hr <- df_6hr[df_6hr$DP > 19, ]


df_6hr$pVal <- apply(df_6hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[16]),
                                            n = sum(as.numeric(x[16]), 
                                                    as.numeric(x[17])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_6hr$FDR <- p.adjust(p = df_6hr$pVal,
                       method = "fdr")
df_6hr <- df_6hr[order(df_6hr$pVal), ]

sum(df_6hr$FDR < 0.05) # 267

sum(df_0hr$ID[df_0hr$FDR < 0.05] %in% df_6hr$ID[df_6hr$FDR < 0.05])
df_0hr$ID[df_0hr$FDR < 0.05][df_0hr$ID[df_0hr$FDR < 0.05] %in% 
                               df_6hr$ID[df_6hr$FDR < 0.05]]

ggplot(df_0hr,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  ggtitle("NEFMp_CUXm_Glut, 6hr, DP >= 20,\n6921 SNPs of 267 FDR < 0.05") +
  ylab("-log10P") +
  theme_classic()
