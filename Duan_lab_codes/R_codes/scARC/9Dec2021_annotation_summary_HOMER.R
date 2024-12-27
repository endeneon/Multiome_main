# Chuxuan Li 12/9/2021
# plot bar graphs using the annotation summary produced by HOMER for peaks in 
#each cell type at 0hr and the differentiallly accessible peaks for each cell
#type comparing 1hr and 0hr to show that the compositions of these sets of peaks
#changed

# init
library(ggplot2)

# test
GABA_1v0 <- read.csv("../new_peak_set_plots/annotation_summary_bargraphs/GABA_1v0.csv",
                     header = T, sep = ",")
#rownames(GABA_1v0) <- GABA_1v0$Annotation
ggplot(GABA_1v0, aes(x = Annotation, y = log(Number_of_peaks))) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  ggtitle("GABA 1hr vs 0hr")


# read in the data
file_names <- c("GABA_1v0.csv", "GABA_all.csv", 
                "NmCp_1v0.csv", "NmCp_all.csv", 
                "NpCm_1v0.csv", "NpCm_all.csv")

for (i in 1:length(file_names)){
  n <- file_names[i]
  file_name <- paste0("../new_peak_set_plots/annotation_summary_bargraphs/", 
                      n)
  data <- read.csv(file = file_name, 
                   header = T, sep = ",")
  
  if (grepl(pattern = "all", x = n, fixed = T)){
    
    if (grepl(pattern = "GABA", x = n, fixed = T)){
      
      plot_title <- paste0("annotation of all peaks in ",
                           gsub(pattern = "_all.csv",
                                replacement = "",
                                x = n))
      
    } else {
      
      plot_title <- paste0("annotation of all peaks in ",
                           gsub(pattern = "_all.csv",
                                replacement = "",
                                x = n),
                           " glut")
      
    }
    
    
  } else {
    
    if (grepl(pattern = "GABA", x = n, fixed = T)){
      
      plot_title <- paste0("annotation of differentially accessible peaks in ",
                           gsub(pattern = "_1v0.csv",
                                replacement = "",
                                x = n),
                           " 1hr vs. 0hr")
      
    } else {
      plot_title <- paste0("annotation of differentially accessible peaks in ",
                           gsub(pattern = "_1v0.csv",
                                replacement = "",
                                x = n),
                           " glut 1hr vs. 0hr")
    }
    
  }
  
  print(plot_title)

  pdf(file = paste0(gsub(pattern = ".csv",
                         replacement = "",
                         x = file_name),
                    "_annotation_summary_graph.pdf"))

  p <- ggplot(data, aes(x = Annotation, y = log(Number_of_peaks))) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    ggtitle(plot_title)

  print(p)

  dev.off()
  
}

