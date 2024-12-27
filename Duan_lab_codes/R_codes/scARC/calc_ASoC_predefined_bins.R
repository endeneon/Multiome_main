# Siwei 18 Jan 2022
# calculate ASoC ratio from the aggregation of 5 lines
# use direct counting from predefined bins only



# init
library(ggplot2)
library(readr)

##### GABA #####
## make 0 hr group
df_0hr <- read_delim("raw_import/GABA_fromGVCF__0hr.table", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)
colnames(df_0hr)[c(7, 8)] <- c("REF_C", "ALT_C")

df_0hr <- df_0hr[df_0hr$REF_C > 1, ]
df_0hr <- df_0hr[df_0hr$ALT_C > 1, ]
df_0hr$DP_GF <- df_0hr$REF_C + df_0hr$ALT_C

sum(df_0hr$DP > 19, na.rm = T)
sum(df_0hr$DP_GF > 19, na.rm = T)# 532
max(df_0hr$DP)

qdf_0hr <- df_0hr[df_0hr$DP > 10, ]

df_0hr$pVal <- apply(df_0hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[17]),
                                            n = as.numeric(x[19]),
                                            p = 0.5,
                                            alternative = "t")$p.value))

df_0hr$pVal <- apply(df_0hr, 1, 
                     function(x)(binom.test(x = as.numeric(x[17]),
                                            n = as.numeric(x[6]),
                                            p = 0.5,
                                            alternative = "t")$p.value))

# df_0hr$pVal <- apply(df_0hr, 1, 
#                      function(x)(binom.test(x = as.numeric(x[6]),
#                                             n = sum(as.numeric(x[6]), 
#                                                     as.numeric(x[7])),
#                                             p = 0.5,
#                                             alternative = "t")$p.value))
df_0hr$FDR <- p.adjust(p = df_0hr$pVal,
                       method = "fdr")
df_0hr <- df_0hr[order(df_0hr$pVal), ]

sum(df_0hr$FDR < 0.05) # 247



ggplot(df_0hr,
       aes(x = (REF_C/DP_GF),
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
