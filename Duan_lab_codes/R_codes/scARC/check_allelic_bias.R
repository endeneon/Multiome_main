# check the allelic bias of WASPed VCF file (per line)
# Siwei 09 Nov 2021

# init
library(ggplot2)
library(readr)

# load data
CD_27_input <- read_delim("CD_27_input.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

CD_27_input <- read_delim("CD_27_bowtie2_DP20.table", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

CD_27_input <- read_delim("CD_54_bowtie2_DP20.table", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

CD_27_input <- read_delim("CD_27_2_0_ASoC.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

colnames(CD_27_input)[4:5] <- c("REF_C", "ALT_C")

CD_27_input <- 
  CD_27_input[(CD_27_input$REF_C > 1) & (CD_27_input$ALT_C > 1), ]

CD_27_input$ALT_C <- as.integer(CD_27_input$ALT_C)
CD_27_input <- CD_27_input[!is.na(CD_27_input$ALT_C), ]

CD_27_input$DP <- CD_27_input$REF_C + CD_27_input$ALT_C

CD_27_input <- 
  CD_27_input[CD_27_input$DP > 19, ]

count_matrix <- CD_27_input[, 4:6]
count_matrix <- as.matrix(count_matrix)


binom.apply <- function(x) {
  value_to_return <-
    binom.test(x = x[1], 
               n = x[3],
               p = 0.5,
               alternative = "t")$p.value
  return(value_to_return)
}

CD_27_input$pBinom <-
  apply(count_matrix, 
        MARGIN = 1,
        FUN = binom.apply)
CD_27_input$FDR <-
  p.adjust(CD_27_input$pBinom,
           method = "fdr")

CD_27_input$REF_ratio <- 
  CD_27_input$REF_C / CD_27_input$DP
CD_27_input$`-logP` <-
  0 - log10(CD_27_input$pBinom)

CD_27_input$signif <- "FDR > 0.05"
CD_27_input$signif[CD_27_input$FDR < 0.5] <- "FDR < 0.05"


ggplot(CD_27_input,
       aes(x = REF_ratio,
           y = `-logP`,
           colour = signif)) +
  geom_point(size = 0.5) +
  scale_colour_manual(values = c("red", "black"), 
                        name = "FDR") +
  # guides()
  ylim(0, 5) +
  ylab("-log10P") +
  ggtitle(paste("CD_27_0hr", 
                "REF_count = 8143, ALT_count = 8087",
                sep = "\n")) +
  theme_classic()

sum(CD_27_input$REF_ratio < 0.5)
sum(CD_27_input$REF_ratio > 0.5)
