# Siwei 30 Aug 2024
# count ASoC enrichment in 25-way regions
# generate necessary bed6 files for bash counting

# init
{
  library(readr)
  # library(vcfR)
  library(stringr)
  library(ggplot2)
  
  library(scales)
  
  library(parallel)
  library(future)
  
  library(MASS)
  
  library(RColorBrewer)
  library(grDevices)
  
  library(ggvenn)
}

# param
plan("multisession", workers = 12)
options(expressions = 20000)
options(future.globals.maxSize = 207374182400)
set.seed(42)

#####


# length("a" %in% c("a", "b","c"))
# length(c("a", "b","c") %in% "a")

hg38_genome_size <- 3234830000

df_raw_data <-
  read.table(file = "ASoC_25_ways_sum_dir/intersect_table.txt",
             header = F, sep = "\t",
             quote = "")
colnames(df_raw_data) <-
  c("Sample_SNP", "SNP_coverage", "Interval_region", "Interval_coverage", "Overlap_coverage")
# df_raw_data <-
#   as.data.frame(rbind(df_raw_data,
#                       c("Z_Enh_all",
#                         )))

df_2_plot <-
  df_raw_data
df_2_plot$Enrichment <-
  (df_2_plot$Overlap_coverage / hg38_genome_size) / 
  ((df_2_plot$Interval_coverage / hg38_genome_size) * (df_2_plot$SNP_coverage / hg38_genome_size))
df_2_plot$log10_enrichment <-
  log10(df_2_plot$Enrichment)
df_2_plot$log10_enrichment[is.infinite(df_2_plot$log10_enrichment)] <- 0

df_2_plot$Cell_type <-
  str_split(df_2_plot$Sample_SNP,
            pattern = "_",
            n = 2,
            simplify = T)[ , 1]
df_2_plot$Cell_time <-
  str_split(df_2_plot$Sample_SNP,
            pattern = "_",
            n = 2,
            simplify = T)[ , 2]

df_2_plot$Region <-
  str_split(df_2_plot$Interval_region,
            pattern = "_",
            n = 2,
            simplify = T)[, 2]
df_2_plot$Region <-
  str_replace_all(df_2_plot$Region,
                  pattern = "^12F",
                  replacement = "Me_F")
df_2_plot$Region <-
  str_replace_all(df_2_plot$Region,
                  pattern = "^12O",
                  replacement = "Me_O")

df_2_plot$Cell_type <-
  factor(df_2_plot$Cell_type,
         levels = c("GABA",
       "nmglut",
       "npglut",
       "atLeastInOneCellType",
       "sharedByAllCellTypes"))

# remove atLeastInOneCellType, shared in all
df_2_plot <-
  df_2_plot[!(df_2_plot$Cell_type %in% "atLeastInOneCellType"), ]
df_2_plot <-
  df_2_plot[!(df_2_plot$Cell_type %in% "sharedByAllCellTypes"), ]
df_2_plot <-
  df_2_plot[!(df_2_plot$Cell_time %in% "at_least_one"), ]

# df_2_plot_backup <- df_2_plot
df_2_plot <- df_2_plot_backup

df_2_plot$Region <-
  str_remove_all(string = df_2_plot$Region,
                 pattern = "^[0-2][0-9]_")
df_2_plot$Region <-
  str_remove_all(string = df_2_plot$Region,
                 pattern = "^hg38_")
df_2_plot$Region <-
  str_remove_all(string = df_2_plot$Region,
                 pattern = "_Hu_gain_hg38$")

df_2_plot$Region[df_2_plot$Region %in% "zEnh_all"] <- "Enh_all"
df_2_plot$Region[df_2_plot$Region %in% "zProm"] <- "Prom"

df_2_plot$Region[df_2_plot$Region %in% "Me_Fpcw_ac"] <- "H3K27ac_F"
df_2_plot$Region[df_2_plot$Region %in% "Me_Fpcw_me2"] <- "H3K4Me2_F"
df_2_plot$Region[df_2_plot$Region %in% "Me_Opcw_ac"] <- "H3K27ac_O"
df_2_plot$Region[df_2_plot$Region %in% "Me_Opcw_me2"] <- "H3K4Me2_O"

df_2_plot$fill <- "black"
df_2_plot$fill[df_2_plot$Region %in% 
                 c("TssA", "PromU","PromD1", "PromD2", "Prom")] <-
  brewer.pal(n = 8,
             name = "Dark2")[1]
df_2_plot$fill[df_2_plot$Region %in% 
                 c("Tx", "Tx5", "Tx3", "TxWk", 
                   "TxReg", "TxEnh5", "TxEnh3", "TxEnhW")] <-
  brewer.pal(n = 8,
             name = "Dark2")[2]
df_2_plot$fill[df_2_plot$Region %in% 
                 c("EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", "EnhAc", "Enh_all")] <-
  brewer.pal(n = 8,
             name = "Dark2")[3]
df_2_plot$fill[df_2_plot$Region %in% 
                 c("Het", "ReprPC", "Quies")] <-
  brewer.pal(n = 8,
             name = "Dark2")[4]
df_2_plot$fill[df_2_plot$Region %in% 
                 c("forebrain_enh", "non_forebrain_enh")] <-
  brewer.pal(n = 8,
             name = "Dark2")[5]
df_2_plot$fill[df_2_plot$Region %in% 
                 c("H3K27ac_F", "H3K4Me2_F", "H3K27ac_O", "H3K4Me2_O")] <-
  brewer.pal(n = 8,
             name = "Dark2")[6]
df_2_plot$fill[df_2_plot$Region %in% 
                 c("DNase", "ZNF_Rpts", "PromP","PromBiv")] <-
  brewer.pal(n = 8,
             name = "Dark2")[7]


df_2_plot$Region <-
  factor(df_2_plot$Region,
         levels = c("TssA", "PromU","PromD1", "PromD2", "Prom",
                    "Tx", "Tx5", "Tx3", "TxWk", 
                    "TxReg", "TxEnh5", "TxEnh3", "TxEnhW", 
                    "EnhA1", "EnhA2", "EnhAF", "EnhW1", "EnhW2", "EnhAc", "Enh_all", 
                    "Het", "ReprPC", "Quies",            
                    "forebrain_enh", "non_forebrain_enh",
                    "H3K27ac_F", "H3K4Me2_F", "H3K27ac_O", "H3K4Me2_O",
                    "DNase", "ZNF_Rpts",         
                    "PromP","PromBiv"))



ggplot(df_2_plot,
       aes(x = Region,
           y = log10_enrichment)) +
  geom_bar(fill = df_2_plot$fill,
           stat = "identity",
           linewidth = 0.2) +
  geom_text(aes(label = Overlap_coverage),
            # colour = df_2_plot$fill,
            angle = 270,
            inherit.aes = T,
            vjust = 0.5,
            hjust = 0,
            colour = "black") +
  theme_classic() +
  # theme_minimal() +
  # theme_bw() +
  theme(axis.text.x = element_text(angle = 270,
                                   vjust = 0.5,
                                   hjust = 0,
                                   colour = "black")) +
  facet_wrap(Cell_type ~ Cell_time)

df_2_plot$GREAT_p <-
  unlist(apply(X = as.matrix(df_2_plot), 
               MARGIN = 1, 
               FUN = function(x) {
                 0 - log10(binom.test(x = as.numeric(x[5]),
                                      n = as.numeric(x[2]),
                                      p = as.numeric(x[4]) / hg38_genome_size,
                                      alternative = "g")$p.value)
               }))

# unlist(apply(X = as.matrix(df_2_plot), 
#              MARGIN = 1, 
#              FUN = function(x) {
#                return(x[2])
#              }))
df_2_plot$GREAT_p[is.infinite(df_2_plot$GREAT_p)] <- 500

ggplot(df_2_plot,
       aes(x = Region,
           y = GREAT_p)) +
  geom_bar(aes(fill = factor(ifelse(test = as.numeric(factor(df_2_plot$Region)) %% 2,
                                    yes = 0,
                                    no = 1)),
               colour = factor(as.numeric(factor(df_2_plot$Cell_type)))),
           stat = "identity",
           linewidth = 0.2) +
  geom_text(aes(label = df_2_plot$Overlap_coverage),
            angle = 270,
            inherit.aes = T,
            vjust = 0.5,
            hjust = 0,
            colour = "black") +
  scale_fill_manual(values = c("pink", "lightblue"),
                    # values = brewer.pal(n = 3,
                    #                     name = "Dark2")[c(1, 2)],
                    guide = "none") +
  scale_colour_manual(values = brewer.pal(n = 5,
                                          name = "Dark2"),
                      guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270,
                                   vjust = 0.5,
                                   hjust = 0,
                                   colour = "black")) +
   facet_grid(Cell_type ~ Cell_time)
  # facet_wrap(~ Cell_time)

ggplot(df_2_plot,
       aes(x = Region,
           y = GREAT_p)) +
  geom_bar(fill = df_2_plot$fill,
           stat = "identity",
           linewidth = 0.2) +
  geom_text(aes(label = Overlap_coverage),
            # colour = df_2_plot$fill,
            angle = 270,
            inherit.aes = T,
            vjust = 0.5,
            hjust = 0,
            colour = "black") +
  theme_classic() +
  # theme_minimal() +
  # theme_bw() +
  theme(axis.text.x = element_text(angle = 270,
                                   vjust = 0.5,
                                   hjust = 0,
                                   colour = "black")) +
  facet_wrap(Cell_type ~ Cell_time)

