markerBarPlot <- function(data_lsdt, geneCol_v = "Gene", y_v = "logFC", fillCol_v = "Class", colors_v = NULL, 
                          markers_dt, markerName_v, expName_v = NULL, title_v = NULL) {
  #' Marker Bar Plot
  #' @description make a barplot displaying expression values for specific markers
  #' @param data_dt list of data.tables of proteomics DE results (can also just be one table)
  #' @param geneCol_v column in data_dt containing gene names
  #' @param y_v column in data_dt that has the y variable
  #' @param fillCol_v variable to color bars with
  #' @param markers_dt data.table of markers to check for in data_dt. must have geneCol_v
  #' @param markerName_v name of markers to include in plot
  #' @param expName_v name of phenotype to include in plot. Only used if data_lsdt is a single data.table, then it is required.
  #' @param title_v optional character vector to prepend to the title.
  #' @return ggplot
  #' @export
  #' 
  
  ### If only one data.table is provided, turn it into a list
  if (!is.logical(all.equal(class(data_lsdt), "list"))) {
    tmp <- copy(data_lsdt)
    data_lsdt <- list()
    data_lsdt[[expName_v]] <- tmp
  } # fi not list
  
  ### Wrangle data
    ### Merge with markers
      ### subset for marker genes only
      ### add fill column if present
    ### Factorize genes to keep in correct order
    ### Clean and melt
  plot_lsdt <- sapply(names(data_lsdt), function(x) {
    xx <- data_lsdt[[x]]
    xx <- merge(markers_dt, xx, by = geneCol_v, sort = F)
    xx <- melt(xx[,mget(c(geneCol_v, y_v, fillCol_v))], measure.vars = y_v)
    xx[[geneCol_v]] <- factor(xx[[geneCol_v]], levels = xx[[geneCol_v]])
    xx$variable <- NULL
    colnames(xx)[ncol(xx)] <- y_v
    xx$Type <- x
    return(xx)}, simplify = F, USE.NAMES = T)
  
  plot_dt <- do.call(rbind, plot_lsdt)
  
  ### Make Plot
  plot_gg <- ggplot(data = plot_dt, aes(x = !!sym(geneCol_v), y = !!sym(y_v), fill = !!sym(fillCol_v))) +
    geom_bar(stat = "identity") +
    my_theme() + angle_x() +
    ggtitle(paste0(title_v, " Murray ", markerName_v, " Expression"))
  
  if (length(plot_lsdt) > 1) {
    plot_gg <- plot_gg + facet_wrap(~Type, ncol = 2, nrow = ceiling(length(plot_lsdt)/2))
  } else {
    plot_gg <- plot_gg + ggtitle(paste0(title_v, " Murray ", markerName_v, " Expression", names(data_lsdt)[1]))
  }
  
  if (!is.null(fillCol_v)) {
    plot_gg <- plot_gg + scale_fill_manual(values = colors_v, breaks = names(colors_v))
  }
  
  return(plot_gg)
  
} # markerBarPlot