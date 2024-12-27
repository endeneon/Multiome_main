# Siwei 13 Jan 2022
# calculate ASoC ratio from the aggregation of CD27 aggregated call

# init
library(ggplot2)
library(readr)

CD_27_df <- read_delim("CD_27_all_time_merged_WASPed_scARC_het.txt", 
                       delim = "\t", escape_double = FALSE,
                       col_names = F,
                       trim_ws = TRUE)
colnames(CD_27_df) <- c("CHR", "POS", "ID", "REF_C", "ALT_C")

CD_27_df$DP <- CD_27_df$REF_C + CD_27_df$ALT_C
CD_27_df$pVal <- apply(CD_27_df, 1, 
                       function(x)(binom.test(x = as.numeric(x[4]),
                                              n = sum(as.numeric(x[4]), 
                                                      as.numeric(x[5])),
                                              p = 0.5,
                                              alternative = "t")$p.value))
CD_27_df$FDR <- p.adjust(p = CD_27_df$pVal,
                         method = "fdr")

ggplot(CD_27_df,
       aes(x = (REF_C/DP),
           y = (0 - log10(pVal)))) +
  geom_point(size = 0.5,
             aes(colour = factor(ifelse(FDR < 0.05, 
                                        "FDR < 0.05", 
                                        "FDR > 0.05")))) +
  scale_color_manual(values = c("red", "black")) +
  labs(colour = "FDR value") +
  # ggtitle("GABA, 0hr, DP >= 10, 10161 SNPs") +
  ylab("-log10P") +
  theme_classic()
