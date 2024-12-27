# Chuxuan Li 01/10/2022
# Using the output of 05Jan2022_time_stimulation_related_peaks_by_celltype.R
#plo dotplots showing for four cell types/three cell types, the expression of BDNF
#and percentage of cells expressing BDNF


# init ####
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(future)

library(readr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(graphics)
library(stringr)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(GenomicFeatures)
library(AnnotationDbi)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)


set.seed(1105)



# BDNF only ####
# create dataframe
df_to_plot_1v0 <- data.frame(cell_type = c("NpCm_glut", "NmCp_glut", "GABA"),
                             p_val = c(2.44e-09, NA, NA),
                             avg_log2FC = c(0.102464771, 0, 0),
                             pct.1 = c(0.03, 0, 0),
                             pct.2 = c(0.013, 0, 0),
                             p_val_adj = c(0.000272767, NA, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("NpCm_glut", "NmCp_glut", "GABA"),
                             p_val = c(1.78e-08, 9.1e-09, NA),
                             avg_log2FC = c(0.093871895, 0.14671314, 0),
                             pct.1 = c(0.029, 0.026, 0),
                             pct.2 = c(0.013, 0.004, 0),
                             p_val_adj = c(0.001983788, 0.00101486, NA))

# make a bubble plot
df_to_plot <- rbind(df_to_plot_1v0,
                    df_to_plot_6v0)
df_to_plot$source <- c(rep_len("1v0hr", length.out = length(df_to_plot_1v0$cell_type)),
                       rep_len("6v0hr", length.out = length(df_to_plot_6v0$cell_type)))


ggplot(df_to_plot,
       aes(x = as.factor(source),
           y = as.factor(cell_type),
           size = pct.1 * 100,
           fill = 2 ^ avg_log2FC * pct.1 / pct.2)) +
  geom_point(shape = 21) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "red", "darkred")) +
  theme_classic() 



# four cell types ####

#ATP2A2
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(NA, NA, NA, NA),
                             pct.1 = c(NA, NA, NA, NA),
                             pct.2 = c(NA, NA, NA, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.027857631, 0.041902501, NA, NA),
                             pct.1 = c(0.009, 0.01, NA, NA),
                             pct.2 = c(0.005, 0.004, NA, NA))

#BDNF 
# chr11 27770354	27771236
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(NA, NA, 0.102464771, NA),
                             pct.1 = c(NA, NA, 0.03, NA),
                             pct.2 = c(NA, NA, 0.013, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(NA, 0.14671314, 0.093871895, NA),
                             pct.1 = c(NA, 0.026, 0.029, NA),
                             pct.2 = c(NA, 0.004, 0.013, NA))

#EGR1
#chr5    138469785       138470701
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.195929101, 0.1808981, 0.137601656, 0.121791281),
                             pct.1 = c(0.059, 0.072, 0.057, 0.092),
                             pct.2 = c(0.018, 0.031, 0.025, 0.061))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(NA, 0.070983828, NA, NA),
                             pct.1 = c(NA, 0.047, NA, NA),
                             pct.2 = c(NA, 0.031, NA, NA))

#FOS
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.056708944, 0.077558196, NA, NA),
                             pct.1 = c(0.162, 0.171, NA, NA),
                             pct.2 = c(0.14, 0.137, NA, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(NA, 0.039158278, NA, NA),
                             pct.1 = c(NA, 0.165, NA, NA),
                             pct.2 = c(NA, 0.137, NA, NA))


#ETF1
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(NA, 0.029931999, 0.048218465, NA),
                             pct.1 = c(NA, 0.023, 0.033, NA),
                             pct.2 = c(NA, 0.018, 0.022, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(NA, NA, NA, NA),
                             pct.1 = c(NA, NA, NA, NA),
                             pct.2 = c(NA, NA, NA, NA))

#MMP1
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.017320692, 0.059959005, 0.035506167, NA),
                             pct.1 = c(0.008, 0.019, 0.028, NA),
                             pct.2 = c(0.005, 0.008, 0.022, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.053316639, 0.10618522, 0.054392112, 0.04109853),
                             pct.1 = c(0.014, 0.027, 0.033, 0.017),
                             pct.2 = c(0.005, 0.008, 0.022, 0.008))

#NPAS4
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.23125368, 0.219222734, 0.222985304, NA),
                             pct.1 = c(0.143, 0.145, 0.147, NA),
                             pct.2 = c(0.071, 0.075, 0.071, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.220993672, 0.169692108, 0.184128655, 0.035536597),
                             pct.1 = c(0.143, 0.13, 0.135, 0.064),
                             pct.2 = c(0.071, 0.075, 0.071, 0.047))
		
		
#NPAS4_2
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.161067779, 0.113855974, 0.142019554, 0.066142998),
                             pct.1 = c(0.063, 0.058, 0.081, 0.095),
                             pct.2 = c(0.025, 0.031, 0.043, 0.076))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.17054478, 0.164151371, 0.110231614, 0.081566271),
                             pct.1 = c(0.067, 0.074, 0.072, 0.105),
                             pct.2 = c(0.025, 0.031, 0.043, 0.076))

#NR4A1
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.367708639, 0.436440371, 0.341603452, 0.260307923),
                             pct.1 = c(0.181, 0.218, 0.199, 0.271),
                             pct.2 = c(0.053, 0.056, 0.065, 0.163))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.289279739, 0.258080532, 0.261958857, 0.217027089),
                             pct.1 = c(0.155, 0.154, 0.165, 0.273),
                             pct.2 = c(0.053, 0.056, 0.065, 0.163))
		
		
#NR4A1_2
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.230019339, 0.27466292, 0.205806565, 0.223923737),
                             pct.1 = c(0.136, 0.151, 0.132, 0.172),
                             pct.2 = c(0.067, 0.062, 0.064, 0.097))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.08166616, 0.109687278, 0.098331862, 0.103014975),
                             pct.1 = c(0.093, 0.101, 0.097, 0.139),
                             pct.2 = c(0.067, 0.062, 0.064, 0.097))

#RGS2
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.034798515, 0.094060914, NA, NA),
                             pct.1 = c(0.118, 0.126, NA, NA),
                             pct.2 = c(0.106, 0.094, NA, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.041126712, 0.055278912, NA, NA),
                             pct.1 = c(0.125, 0.117, NA, NA),
                             pct.2 = c(0.106, 0.094, NA, NA))
                             

#ZSWIM6
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.069003373, 0.055302832, 0.101127731, NA),
                             pct.1 = c(0.019, 0.027, 0.047, NA),
                             pct.2 = c(0.008, 0.016, 0.026, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.061571643, 0.055481837, 0.055672587, 0.038725044),
                             pct.1 = c(0.019, 0.028, 0.037, 0.025),
                             pct.2 = c(0.008, 0.016, 0.026, 0.016))

#ZSWIM6_2
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.091410967, 0.120384368, 0.06051371, NA),
                             pct.1 = c(0.083, 0.091, 0.103, NA),
                             pct.2 = c(0.055, 0.056, 0.083, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.05514496, 0.098148657, NA, NA),
                             pct.1 = c(0.075, 0.086, NA, NA),
                             pct.2 = c(0.055, 0.056, NA, NA))
		
#ZSWIM6_3
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.058853561, 0.055499825, 0.078881739, 0.013375744),
                             pct.1 = c(0.023, 0.02, 0.046, 0.027),
                             pct.2 = c(0.013, 0.011, 0.03, 0.025))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.05269498, 0.065566417, 0.068932811, 0.06125984),
                             pct.1 = c(0.023, 0.023, 0.045, 0.04),
                             pct.2 = c(0.013, 0.011, 0.03, 0.025))
		
#PPP1R13B
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.046755439, NA, 0.042915065, NA),
                             pct.1 = c(0.013, NA, 0.014, NA),
                             pct.2 = c(0.006, NA, 0.008, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.064231542, 0.091007689, 0.051150378, NA),
                             pct.1 = c(0.016, 0.019, 0.016, NA),
                             pct.2 = c(0.006, 0.006, 0.008, NA))		
		
		
# make a bubble plot
df_to_plot <- rbind(df_to_plot_1v0,
                    df_to_plot_6v0)
df_to_plot$source <- c(rep_len("1v0hr", length.out = length(df_to_plot_1v0$cell_type)),
                       rep_len("6v0hr", length.out = length(df_to_plot_6v0$cell_type)))
df_to_plot$cell_type <- factor(df_to_plot$cell_type, 
                              levels = rev(c("GABA", "NmCp_glut", "NpCm_glut", "NPC")))


ggplot(df_to_plot,
       aes(x = as.factor(source),
           y = as.factor(cell_type),
           size = pct.1 * 100,
           fill = 2 ^ avg_log2FC * pct.1 / pct.2)) +
  geom_point(shape = 21) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "red", "darkred")) +
  theme_grey(base_size = 10) +
  ggtitle("PPP1R13B")



# three cell types ####
#ATP2A2
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "glut", "NPC"),
                             avg_log2FC = c(NA, NA, NA),
                             pct.1 = c(NA, NA, NA),
                             pct.2 = c(NA, NA, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "glut", "NPC"),
                             avg_log2FC = c(0.027857631, 0.02150795, NA),
                             pct.1 = c(0.009, 0.008, NA),
                             pct.2 = c(0.005, 0.005, NA))
		
#ZSWIM6
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "glut", "NPC"),
                             avg_log2FC = c(0.069003373, 0.094513104, NA),
                             pct.1 = c(0.019, 0.043, NA),
                             pct.2 = c(0.008, 0.023, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "glut", "NPC"),
                             avg_log2FC = c(0.061571643, 0.058661602, NA),
                             pct.1 = c(0.019, 0.035, NA),
                             pct.2 = c(0.008, 0.023, NA))

#ZSWIM6_2
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "glut", "NPC"),
                             avg_log2FC = c(0.091410967, 0.077766124, NA),
                             pct.1 = c(0.083, 0.1, NA),
                             pct.2 = c(0.055, 0.075, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "glut", "NPC"),
                             avg_log2FC = c(0.05514496, 0.044930443, NA),
                             pct.1 = c(0.075, 0.09, NA),
                             pct.2 = c(0.055, 0.075, NA))


#BDNF
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "glut", "NPC"),
                             avg_log2FC = c(NA, 0.068396389, NA),
                             pct.1 = c(NA, 0.031, NA),
                             pct.2 = c(NA, 0.017, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "glut", "NPC"),
                             avg_log2FC = c(0.02675838, 0.084749955, NA),
                             pct.1 = c(0.022, 0.034, NA),
                             pct.2 = c(0.017, 0.017, NA))

		
#EGR1
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "glut", "NPC"),
                             avg_log2FC = c(0.195929101, 0.144720968, 0.121791281),
                             pct.1 = c(0.059, 0.06, 0.092),
                             pct.2 = c(0.018, 0.027, 0.061))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "glut", "NPC"),
                             avg_log2FC = c(NA, NA, NA),
                             pct.1 = c(NA, NA, NA),
                             pct.2 = c(NA, NA, NA))

#FOS
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "glut", "NPC"),
                             avg_log2FC = c(0.056708944, 0.026597994, NA),
                             pct.1 = c(0.162, 0.17, NA),
                             pct.2 = c(0.14, 0.157, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "glut", "NPC"),
                             avg_log2FC = c(NA, NA, NA),
                             pct.1 = c(NA, NA, NA),
                             pct.2 = c(NA, NA, NA))
		

		

df_to_plot <- rbind(df_to_plot_1v0,
                    df_to_plot_6v0)
df_to_plot$source <- c(rep_len("1v0hr", length.out = length(df_to_plot_1v0$cell_type)),
                       rep_len("6v0hr", length.out = length(df_to_plot_6v0$cell_type)))
df_to_plot$cell_type <- factor(df_to_plot$cell_type, 
                               levels = rev(c("GABA", "glut", "NPC")))


ggplot(df_to_plot,
       aes(x = as.factor(source),
           y = as.factor(cell_type),
           size = pct.1 * 100,
           fill = 2 ^ avg_log2FC * pct.1 / pct.2)) +
  geom_point(shape = 21) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "red", "darkred")) +
  theme_grey() 


# Using MAST
#BDNF
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(NA, 0.086334061, 0.062924852, NA),
                             pct.1 = c(NA, 0.031, 0.031, NA),
                             pct.2 = c(NA, 0.016, 0.018, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(NA, 0.110311057, 0.077020821, NA),
                             pct.1 = c(NA, 0.037, 0.033, NA),
                             pct.2 = c(NA, 0.016, 0.018, NA))

#FOS
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.127177523, 0.280515313, 0.206142961, 0.098411203),
                             pct.1 = c(0.072, 0.108, 0.098, 0.075),
                             pct.2 = c(0.04, 0.033, 0.043, 0.048))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.10232409, 0.172720994, 0.161107388, 0.082047247),
                             pct.1 = c(0.068, 0.08, 0.087, 0.076),
                             pct.2 = c(0.04, 0.033, 0.043, 0.048))

#NPAS4
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.23125368, 0.219222734, 0.222985304, NA),
                             pct.1 = c(0.143, 0.145, 0.147, NA),
                             pct.2 = c(0.071, 0.075, 0.071, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.220993672, 0.169692108, 0.184128655, 0.035536597),
                             pct.1 = c(0.143, 0.13, 0.135, 0.064),
                             pct.2 = c(0.071, 0.075, 0.071, 0.047))

#EGR1
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.195929101, 0.1808981, 0.137601656, 0.121791281),
                             pct.1 = c(0.059, 0.072, 0.057, 0.092),
                             pct.2 = c(0.018, 0.031, 0.025, 0.061))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.011246587, 0.02567091, NA, 0.001519996),
                             pct.1 = c(0.016, 0.047, NA, 0.092),
                             pct.2 = c(0.018, 0.031, NA, 0.061))

#NR4A1
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.230019339,0.27466292, 0.205806565,0.223923737),
                             pct.1 = c(0.136,0.151,0.132,0.172),
                             pct.2 = c(0.067,0.062,0.064, 0.097))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.08166616,
                                            0.109687278,
                                            0.098331862,
                                            0.103014975
                             ),
                             pct.1 = c(0.093,
                                       0.101,
                                       0.097,
                                       0.139
                             ),
                             pct.2 = c(0.067,
                                       0.062,
                                       0.064,
                                       0.097
                             ))

#MMP1
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.027347173, NA, 0.039549921, NA),
                             pct.1 = c(0.007, NA, 0.021, NA),
                             pct.2 = c(0.003, NA, 0.013, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.049629803, 0.079399049, 0.049511271, 0.066833924),
                             pct.1 = c(0.011, 0.022, 0.023, 0.034),
                             pct.2 = c(0.003, 0.009, 0.013, 0.021))


#RGS2
df_to_plot_1v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(NA, 0.094060914, NA, NA),
                             pct.1 = c(NA, 0.126, NA, NA),
                             pct.2 = c(NA, 0.094, NA, NA))
df_to_plot_6v0 <- data.frame(cell_type = c("GABA", "NmCp_glut", "NpCm_glut", "NPC"),
                             avg_log2FC = c(0.041126712, 0.055278912, NA, NA),
                             pct.1 = c(0.125, 0.117, NA, NA),
                             pct.2 = c(0.106, 0.094, NA, NA))

# make a bubble plot
df_to_plot <- rbind(df_to_plot_1v0,
                    df_to_plot_6v0)
df_to_plot$source <- c(rep_len("1v0hr", length.out = length(df_to_plot_1v0$cell_type)),
                       rep_len("6v0hr", length.out = length(df_to_plot_6v0$cell_type)))
df_to_plot$cell_type <- factor(df_to_plot$cell_type, 
                               levels = rev(c("GABA", "NmCp_glut", "NpCm_glut", "NPC")))


ggplot(df_to_plot,
       aes(x = as.factor(source),
           y = as.factor(cell_type),
           size = pct.1 * 100,
           fill = 2 ^ avg_log2FC * pct.1 / pct.2)) +
  geom_point(shape = 21) +
  ylab("cluster") +
  xlab("time points") +
  scale_fill_gradientn(colours = c("white", "red", "darkred")) +
  theme_grey(base_size = 10) +
  ggtitle("BDNF")
