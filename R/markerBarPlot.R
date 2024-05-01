markerBarPlot <- function(data_dt, geneCol_v = "Gene", y_v = "logFC", fillCol_v = "Class", colors_v = NULL, markers_dt, markerName_v, expName_v) {
  #' Marker Bar Plot
  #' @description make a barplot displaying expression values for specific markers
  #' @param data_dt data.table of proteomics DE results
  #' @param geneCol_v column in data_dt containing gene names
  #' @param y_v column in data_dt that has the y variable
  #' @param fillCol_v variable to color bars with
  #' @param markers_dt data.table of markers to check for in data_dt. must have geneCol_v
  #' @param markerName_v name of markers to include in plot
  #' @param expName_v name of phenotype to include in plot
  #' @return ggplot
  #' @export
  
  ### Subset the data
  data_dt <- data_dt[get(geneCol_v) %in% markers_dt[[geneCol_v]],]
  markers_dt <- markers_dt[get(geneCol_v) %in% data_dt[[geneCol_v]]]
  
  ### Merge with markers to get color class
  if (!is.null(fillCol_v)) {
    data_dt <- merge(markers_dt[,mget(c(geneCol_v, fillCol_v))], data_dt, by = geneCol_v, sort = F)
  } # fi
  
  ### Sort genes
  data_dt[[geneCol_v]] <- factor(data_dt[[geneCol_v]], levels = markers_dt[[geneCol_v]])
  
  ### Make Plot
  plot_gg <- ggplot(data = data_dt, aes(x = !!sym(geneCol_v), y = !!sym(y_v), fill = !!sym(fillCol_v))) +
    geom_bar(stat = "identity") +
    my_theme() + angle_x() +
    ggtitle(paste0("Murray ", markerName_v, " Expression in\n", expName_v))
  
  if (!is.null(fillCol_v)) {
    plot_gg <- plot_gg + scale_fill_manual(values = colors_v, breaks = names(colors_v))
  }
  
  return(plot_gg)
  
} # markerBarPlot