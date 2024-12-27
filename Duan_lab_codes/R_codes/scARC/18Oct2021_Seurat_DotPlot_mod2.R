# Chuxuan Li 10/18/2021
# rewrite DotPlot function to adapt to split.by != NULL
library(stringr)
library(fields)

PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

DPnew2 <- 
  function (object, assay = NULL, features, cols = c("lightgrey", 
                                                     "blue"), 
            col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
            idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
            scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
  {
    assay <- assay %||% DefaultAssay(object = object) # 
    DefaultAssay(object = object) <- assay # set default assay
    split.colors <- !is.null(x = split.by) && !any(cols %in%  # if split by something, then TRUE
                                                     rownames(x = brewer.pal.info)) # if cols not in brewer, T
    scale.func <- switch(EXPR = scale.by, size = scale_size, 
                         radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    feature.groups <- NULL # initialize feature list
    if (is.list(features) | any(!is.na(names(features)))) { # if features input is a list or has at least 1 non-NA
      feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                          FUN = function(x) {
                                            return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                          }))
      if (any(is.na(x = feature.groups))) {
        warning("Some feature groups are unnamed.", call. = FALSE, 
                immediate. = TRUE)
      }
      features <- unlist(x = features)
      names(x = feature.groups) <- features
    }
    # Get cell names grouped by identity class
    cells <- unlist(x = CellsByIdentities(object = object, idents = idents)) 
    data.features <- FetchData(object = object, vars = features, 
                               cells = cells)
    data.features$id <- if (is.null(x = group.by)) {
      Idents(object = object)[cells, drop = TRUE]
    }
    else {
      object[[group.by, drop = TRUE]][cells, drop = TRUE]
    }
    if (!is.factor(x = data.features$id)) {
      data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) { # if split by something, split the data
      splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
      if (split.colors) {
        if (length(x = unique(x = splits)) > length(x = cols)) {
          stop("Not enough colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
      }
      data.features$id <- paste(data.features$id, splits, 
                                sep = "_")
      unique.splits <- unique(x = splits)
      id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                          "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 
                                1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
        return(mean(x = expm1(x = x)))
      })
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                       threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    if (cluster.idents) {
      mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                                 FUN = unlist))
      mat <- scale(x = mat)
      id.levels <- id.levels[hclust(d = dist(x = mat))$order]
    }
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    
    #print(id.levels)
    # id. levels = cell type name + _ + xhr
    if (!is.null(x = id.levels)) {
      data.plot$id <- factor(x = data.plot$id, levels = id.levels) # data.plot$id has levels (unique ids)
    }
    #print(data.plot$id)
    if (length(x = levels(x = data.plot$id)) == 1) {
      scale <- FALSE
      warning("Only one identity present, the expression values will be not scaled", 
              call. = FALSE, immediate. = TRUE)
    }
    #print(unique(x = data.plot$features.plot))
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                             FUN = function(x) {
                               data.use <- data.plot[data.plot$features.plot == 
                                                       x, "avg.exp"]
                               if (scale) {
                                 data.use <- scale(x = data.use)
                                 data.use <- MinMax(data = data.use, min = col.min, 
                                                    max = col.max)
                               }
                               else {
                                 data.use <- log1p(x = data.use)
                               }
                               return(data.use)
                             })

    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (split.colors) {
      avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                           breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                      levels = features)
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (split.colors) {
      splits.use <- str_sub(string = as.character(x = data.plot$id), start = -3L)
      # splits.use <- vapply(X = as.character(x = data.plot$id), # input is all ids
      #                      FUN = gsub, # gsub uses regex to alter strings
      #                      FUN.VALUE = character(length = 1L),  
      #                      pattern = paste0("^((", # ^ is start, () is to group
      #                                       paste(sort(x = levels(x = object), 
      #                                                  decreasing = TRUE), 
      #                                             collapse = "|"), 
      #                                       ")_)"), 
      #                      replacement = "", 
      #                      USE.NAMES = FALSE)
      print(splits.use)
      data.plot$colors <- mapply(FUN = function(color, value) {
        return(colorRampPalette(colors = c("grey", color))(20)[value])
      }, color = cols[splits.use], value = avg.exp.scaled)
      #print(data.plot$colors) # grey mostly
      #print(cols[splits.use]) # correct color names
      #print(data.plot$avg.exp.scaled) # why mainly 1s?
      
    }
    # there are 2 factors to arrange the data in the following vector: group_by and split_by
    print(avg.exp.scaled) # this vector is written in the order of group1_split1, group2_split1, ..., group1_split2, group2_split2,...
    color.by <- ifelse(test = split.colors, yes = "colors", 
                       no = "avg.exp.scaled")
    if (!is.na(x = scale.min)) {
      data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
      data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    if (!is.null(x = feature.groups)) {
      data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                         levels = unique(x = feature.groups))
    }
    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                          y = "id")) + 
      geom_point(mapping = aes_string(size = "pct.exp", 
                                      color = color.by)) + 
      scale.func(range = c(0, dot.scale), 
                 limits = c(scale.min, scale.max)) + 
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank()) +
      
      guides(size = guide_legend(title = "Percent Expressed")) + 
      labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                      yes = "Identity", no = "Split Identity")) + 
      cowplot::theme_cowplot()
    if (!is.null(x = feature.groups)) {
      plot <- plot + facet_grid(facets = ~feature.groups, 
                                scales = "free_x", space = "free_x", switch = "y") + 
        theme(panel.spacing = unit(x = 1, units = "lines"), 
              strip.background = element_blank())
    }
    if (split.colors) {
      plot <- plot + 
        scale_color_identity() 
    }
    else if (length(x = cols) == 1) {
      plot <- plot + scale_color_distiller(palette = cols)
    }
    else {
      plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (!split.colors) {
      plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    }
    return(plot)
  }
