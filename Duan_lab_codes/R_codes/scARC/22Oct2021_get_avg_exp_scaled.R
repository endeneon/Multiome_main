# Chuxuan Li Siwei Zhang 10/22/2021
# rewrite DotPlot function to adapt to split.by != NULL
library(stringr)
library(fields)

PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

GetAvgExpScaled <- 
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
    print(avg.exp.scaled)
    avg_exp_scaled_df <- data.frame(avg.exp.scaled, 
                                 stringsAsFactors = F)
    rownames(avg_exp_scaled_df) <- c("inh_0hr", "ext_0hr",
                                     "inh_1hr", "ext_1hr",
                                     "inh_6hr", "ext_6hr")
    return(avg_exp_scaled_df)

  }
