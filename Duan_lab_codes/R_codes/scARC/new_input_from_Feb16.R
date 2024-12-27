# Siwei 03 Jan 2022
# calculate ASoC ratio from the aggregation of 5 lines

# init
library(ggplot2)
library(readr)

##### GABA #####
## make 0 hr group
df_0hr <- read_delim("Feb16_raw_input/GABA__0hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_0hr$REF_C <- apply(df_0hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_0hr$ALT_C <- apply(df_0hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_0hr <- df_0hr[df_0hr$REF_C > 1, ]
df_0hr <- df_0hr[df_0hr$ALT_C > 1, ]
df_0hr$DP <- df_0hr$REF_C + df_0hr$ALT_C

sum(df_0hr$DP > 19) # 10161

df_0hr <- df_0hr[df_0hr$DP > 59, ]

max(df_0hr$DP)
df_0hr$pVal <- apply(df_0hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[16]),
                                            n = sum(as.numeric(x[16]), 
                                                    as.numeric(x[17])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_0hr$FDR <- p.adjust(p = df_0hr$pVal,
                       method = "fdr")

sum((df_0hr$REF_C/df_0hr$DP) > 0.5)
sum((df_0hr$REF_C/df_0hr$DP) < 0.5)

sum(df_0hr$FDR < 0.05)

## make 1 hr group
df_1hr <- read_delim("Feb16_raw_input/GABA__1hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_1hr$REF_C <- apply(df_1hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_1hr$ALT_C <- apply(df_1hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_1hr <- df_1hr[df_1hr$REF_C > 1, ]
df_1hr <- df_1hr[df_1hr$ALT_C > 1, ]

df_1hr$DP <- df_1hr$REF_C + df_1hr$ALT_C
max(df_1hr$DP)
df_1hr <- df_1hr[df_1hr$DP > 59, ]


df_1hr$pVal <- apply(df_1hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[16]),
                                            n = sum(as.numeric(x[16]), 
                                                    as.numeric(x[17])),
                                            p = 0.5,
                                            alternative = "t")$p.value))
df_1hr$FDR <- p.adjust(p = df_1hr$pVal,
                       method = "fdr")

## make 6 hr group
df_6hr <- read_delim("Feb16_raw_input/GABA__6hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)

df_6hr$REF_C <- apply(df_6hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14)]))))
df_6hr$ALT_C <- apply(df_6hr, 1, 
                      function(x)(sum(as.numeric(x[c(6, 8, 10, 12, 14) + 1]))))

df_6hr <- df_6hr[df_6hr$REF_C > 1, ]
df_6hr <- df_6hr[df_6hr$ALT_C > 1, ]

df_6hr$DP <- df_6hr$REF_C + df_6hr$ALT_C
max(df_6hr$DP)
df_6hr <- df_6hr[df_6hr$DP > 59, ]


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

ggplot(df_6hr,
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
