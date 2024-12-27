# Chuxuan Li 02/16/2022
# Take a GO analysis result .txt file, plot bar graphs showing the top 10
#enriched cellular processes and enrichment score

# init
library(ggplot2)
library(readr)


# load data
setwd("/nvmefs/scARC_Duan_018/R_scARC_Duan_018/Analysis_RNAseq_v3_normalize_by_6libs_integrate_by_Seurat/20line_DEG_case_control_by_time")
analysis_0hr_GABA_do <- read.delim("./analysis_0hr_GABA_do.txt")
col_names <- c("term", "REF", "#", "expected", "over_under", 
               "fold_enrichment", "raw_p_value", "FDR")
colnames(analysis_0hr_GABA_do) <- col_names
analysis_0hr_GABA_do$term <- gsub(pattern = '\\(GO:[0-9]+\\)', 
                                  replacement = "", 
                                  x = analysis_0hr_GABA_do$term)
analysis_0hr_GABA_do$fold_enrichment[analysis_0hr_GABA_do$fold_enrichment == " > 100"] <- 100
analysis_0hr_GABA_do$fold_enrichment[analysis_0hr_GABA_do$fold_enrichment == " < 0.01"] <- 0
analysis_0hr_GABA_do$fold_enrichment <- as.numeric(analysis_0hr_GABA_do$fold_enrichment)

# bar graph
ggplot(analysis_0hr_GABA_do[1:10, ], 
       aes(term, log2(fold_enrichment))) +
  geom_col() + 
  geom_text(aes(label = analysis_0hr_GABA_do[1:10, ]$term, hjust = -0.05)) + 
  theme(axis.text.y = element_blank()) +
  scale_y_continuous(limits = c(0, 10)) +
  ggtitle("GO terms enriched in genes with higher expression in Case cell lines \n in 0hr GABA cells") +
  coord_flip()


# GABA 0hr do
analysis_0hr_GABA_do <- read.delim("./analysis_0hr_GABA_do.txt")
colnames(analysis_0hr_GABA_do) <- col_names
analysis_0hr_GABA_do$term <- gsub(pattern = '\\(GO:[0-9]+\\)', 
                                  replacement = "", 
                                  x = analysis_0hr_GABA_do$term)
analysis_0hr_GABA_do$fold_enrichment[analysis_0hr_GABA_do$fold_enrichment == " > 100"] <- 100
analysis_0hr_GABA_do$fold_enrichment[analysis_0hr_GABA_do$fold_enrichment == " < 0.01"] <- 0
analysis_0hr_GABA_do$fold_enrichment <- as.numeric(analysis_0hr_GABA_do$fold_enrichment)

# bar graph
ggplot(analysis_0hr_GABA_do[1:10, ], 
       aes(term, log2(fold_enrichment))) +
  geom_col() + 
  geom_text(aes(label = analysis_0hr_GABA_do[1:10, ]$term, hjust = -0.05)) + 
  theme(axis.text.y = element_blank()) +
  scale_y_continuous(limits = c(0, 10)) +
  ggtitle("GO terms enriched in genes with higher expression in Control cell lines \n in 0hr GABA cells") +
  coord_flip()


# GABA 1hr up
analysis_1hr_GABA_up <- read.delim("./analysis_1hr_GABA_up.txt")
colnames(analysis_1hr_GABA_up) <- col_names
analysis_1hr_GABA_up$term <- gsub(pattern = '\\(GO:[0-9]+\\)', 
                                  replacement = "", 
                                  x = analysis_1hr_GABA_up$term)
analysis_1hr_GABA_up$fold_enrichment[analysis_1hr_GABA_up$fold_enrichment == " > 100"] <- 100
analysis_1hr_GABA_up$fold_enrichment[analysis_1hr_GABA_up$fold_enrichment == " < 0.01"] <- 0
analysis_1hr_GABA_up$fold_enrichment <- as.numeric(analysis_1hr_GABA_up$fold_enrichment)

# bar graph
ggplot(analysis_1hr_GABA_up[1:10, ], 
       aes(term, log2(fold_enrichment))) +
  geom_col() + 
  geom_text(aes(label = analysis_1hr_GABA_up[1:10, ]$term, hjust = -0.05)) + 
  theme(axis.text.y = element_blank()) +
  scale_y_continuous(limits = c(0, 10)) +
  ggtitle("GO terms enriched in genes with higher expression in Case cell lines \n in 1hr GABA cells") +
  coord_flip()


# GABA 1hr do
analysis_1hr_GABA_do <- read.delim("./analysis_1hr_GABA_do.txt")
colnames(analysis_1hr_GABA_do) <- col_names
analysis_1hr_GABA_do$term <- gsub(pattern = '\\(GO:[0-9]+\\)', 
                                  replacement = "", 
                                  x = analysis_1hr_GABA_do$term)
analysis_1hr_GABA_do$fold_enrichment[analysis_1hr_GABA_do$fold_enrichment == " > 100"] <- 100
analysis_1hr_GABA_do$fold_enrichment[analysis_1hr_GABA_do$fold_enrichment == " < 0.01"] <- 0
analysis_1hr_GABA_do$fold_enrichment <- as.numeric(analysis_1hr_GABA_do$fold_enrichment)

# bar graph
ggplot(analysis_1hr_GABA_do[1:10, ], 
       aes(term, log2(fold_enrichment))) +
  geom_col() + 
  geom_text(aes(label = analysis_1hr_GABA_do[1:10, ]$term, hjust = -0.05)) + 
  theme(axis.text.y = element_blank()) +
  scale_y_continuous(limits = c(0, 10)) +
  ggtitle("GO terms enriched in genes with higher expression in Control cell lines \n in 1hr GABA cells") +
  coord_flip()


# GABA 6hr up
analysis_6hr_GABA_up <- read.delim("./analysis_6hr_GABA_up.txt")
colnames(analysis_6hr_GABA_up) <- col_names
analysis_6hr_GABA_up$term <- gsub(pattern = '\\(GO:[0-9]+\\)', 
                                  replacement = "", 
                                  x = analysis_6hr_GABA_up$term)
analysis_6hr_GABA_up$fold_enrichment[analysis_6hr_GABA_up$fold_enrichment == " > 100"] <- 100
analysis_6hr_GABA_up$fold_enrichment[analysis_6hr_GABA_up$fold_enrichment == " < 0.01"] <- 0
analysis_6hr_GABA_up$fold_enrichment <- as.numeric(analysis_6hr_GABA_up$fold_enrichment)

# bar graph
ggplot(analysis_6hr_GABA_up[1:10, ], 
       aes(term, log2(fold_enrichment))) +
  geom_col() + 
  geom_text(aes(label = analysis_6hr_GABA_up[1:10, ]$term, hjust = -0.05)) + 
  theme(axis.text.y = element_blank()) +
  scale_y_continuous(limits = c(0, 10)) +
  ggtitle("GO terms enriched in genes with higher expression in Case cell lines \n in 6hr GABA cells") +
  coord_flip()


# GABA 6hr do
analysis_6hr_GABA_do <- read.delim("./analysis_6hr_GABA_do.txt")
colnames(analysis_6hr_GABA_do) <- col_names
analysis_6hr_GABA_do$term <- gsub(pattern = '\\(GO:[0-9]+\\)', 
                                  replacement = "", 
                                  x = analysis_6hr_GABA_do$term)
analysis_6hr_GABA_do$fold_enrichment[analysis_6hr_GABA_do$fold_enrichment == " > 100"] <- 100
analysis_6hr_GABA_do$fold_enrichment[analysis_6hr_GABA_do$fold_enrichment == " < 0.01"] <- 0
analysis_6hr_GABA_do$fold_enrichment <- as.numeric(analysis_6hr_GABA_do$fold_enrichment)

# bar graph
ggplot(analysis_6hr_GABA_do[1:10, ], 
       aes(term, log2(fold_enrichment))) +
  geom_col() + 
  geom_text(aes(label = analysis_6hr_GABA_do[1:10, ]$term, hjust = -0.05)) + 
  theme(axis.text.y = element_blank()) +
  scale_y_continuous(limits = c(0, 10)) +
  ggtitle("GO terms enriched in genes with higher expression in Control cell lines \n in 6hr GABA cells") +
  coord_flip()
