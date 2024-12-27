# Chuxuan Li 1/12/2023
# Summarize 025 limma DEG MAGMA enrichment results

# init ####
library(readr)
library(stringr)

library(ggplot2)
library(RColorBrewer)

# read results ####
results_paths <- list.files("./Analysis_part2_GRCh38/magma/gsa_out", full.names = T)
results_dfs <- vector("list", length(results_paths))

for (i in 1:length(results_paths)) {
  if (i %% 7 == 1 | i %% 7 == 0) {
    results_dfs[[i]] <- read_table(results_paths[i], skip = 3)
  } else {
    results_dfs[[i]] <- read_table(results_paths[i], skip = 4)
  }
}

namelist <- str_remove(results_paths, "^\\./Analysis_part2_GRCh38/magma/gsa_out/")
namelist <- str_remove(namelist, "\\.gsa.out$")
celltype_list <- str_extract(namelist, "^[A-Za-z]+")
time_list <- str_extract(namelist, "[1|6]v0")
sig_list <- str_remove(str_extract(namelist, "0_[a-z]+"), "0_")
disease_list <- str_extract(namelist, "[A-Za-z0-9]+$")

names(results_dfs) <- namelist

# make dfs for plotting ####
GABA_dfs <- results_dfs[str_detect(names(results_dfs), "GABA")]
nmglut_dfs <- results_dfs[str_detect(names(results_dfs), "nmglut")]
npglut_dfs <- results_dfs[str_detect(names(results_dfs), "npglut")]
addInfo2Df <- function(df_temp, dfs) {
  df_temp$cell_type <- str_extract(names(dfs)[i], "^[A-Za-z]+")
  df_temp$time_points <- str_extract(names(dfs)[i], "[1|6]v0")
  df_temp$disease <- str_extract(names(dfs)[i], "[A-Za-z0-9]+$")
  df_temp$significance <- str_remove(str_extract(names(dfs)[i], "0_[a-z]+"), "0_")
  df_temp$name <- names(dfs)[i]
  if (df_temp$P > 0.05) {
    df_temp$neglogp <- NA
  } else {
    df_temp$neglogp <- (-1) * log10(GABA_df_temp$P)
  }
  return(df_temp)
}
for (i in 1:length(GABA_dfs)) {
  GABA_df_temp <- GABA_dfs[[i]]
  GABA_df_temp <- addInfo2Df(GABA_df_temp, GABA_dfs)
  nmglut_df_temp <- nmglut_dfs[[i]]
  nmglut_df_temp <- addInfo2Df(nmglut_df_temp, nmglut_dfs)
  npglut_df_temp <- npglut_dfs[[i]]
  npglut_df_temp <- addInfo2Df(npglut_df_temp, npglut_dfs)
  
  if (i == 1) {
    GABA_df <- GABA_df_temp
    nmglut_df <- nmglut_df_temp
    npglut_df <- npglut_df_temp
  } else {
    GABA_df <- rbind(GABA_df, GABA_df_temp)
    nmglut_df <- rbind(nmglut_df, nmglut_df_temp)
    npglut_df <- rbind(npglut_df, npglut_df_temp)
  }
}

# plot dotplots for each cell type ####
limit <- max(abs(GABA_df$BETA)) * c(-1, 1)
ggplot(GABA_df, aes(x = disease, y = significance, size = neglogp, color = BETA)) +
  geom_point() +
  facet_grid(rows = vars(time_points), switch = "y") +
  scale_color_gradientn(colors = rev(brewer.pal(9, "RdBu")), limit = limit) +
  labs(color = "Beta", size = "- log10(p)", x = "Diseases", y = "") +
  ggtitle("GABA")

limit <- max(abs(nmglut_df$BETA)) * c(-1, 1)
ggplot(nmglut_df, aes(x = disease, y = significance, size = neglogp, color = BETA)) +
  geom_point() +
  facet_grid(rows = vars(time_points), switch = "y") +
  scale_color_gradientn(colors = rev(brewer.pal(9, "RdBu")), limit = limit) +
  labs(color = "Beta", size = "- log10(p)", x = "Diseases", y = "") +
  ggtitle("nmglut")

limit <- max(abs(npglut_df$BETA)) * c(-1, 1)
ggplot(npglut_df, aes(x = disease, y = significance, size = neglogp, color = BETA)) +
  geom_point() +
  facet_grid(rows = vars(time_points), switch = "y") +
  scale_color_gradientn(colors = rev(brewer.pal(9, "RdBu")), limit = limit) +
  labs(color = "Beta", size = "- log10(p)", x = "Diseases", y = "") +
  ggtitle("npglut")
