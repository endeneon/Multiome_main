# 27 May 2020 Siwei
# recapitulate the scanpy stacked violin plot function in R
# https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/

# init
library(Seurat)
library(future)
library(readr)
library(sctransform)
library(stringr)
library(ggplot2)
library(patchwork)


### main functions #####
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot <- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
         axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(0.8), hjust = rel(0.1),
                                      vjust = rel(0.5), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

modify_vlnplot <- function(obj,
                           feature,
                           pt.size = 0,
                           # plot.margin = unit(c(0,0,0,0), "pt"),
                           plot.margin = unit(c(-3, 0, -3, 0), "pt"),
                           ...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    scale_y_continuous(trans = 'log2') +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(0.8), hjust = rel(0.1),
                                      vjust = rel(0.5), angle = 0),
          plot.margin = plot.margin,
          plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = rel(1)))
  return(p)
}

## extract the max value of the y axis
extract_max <- function(p){
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot <- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs <- purrr::map_dbl(plot_list, extract_max)
  plot_list <- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y)) #+ 
                            # scale_x_discrete(limits = c("0", "1", "6", "7", "12", "19",
                            #                             "3", "4", "9", "18", 
                            #                             "2", "5", "11",
                            #                             "8", "10", "13", "14", "15", "16", "17", "20")))
  
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

