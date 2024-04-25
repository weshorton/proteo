proteoVolcano <- function(data_dt, geneCol_v = "Gene", lfcCol_v, pvalCol_v,
                          lfc_v = 0.5, ident1_v, labelAll_v = F, labelTop_v = 20,
                          labelDir_v = "both", labelSize_v = 1, title_v) {
  #' ProteoVolcano
  #' @description Volcano plot for proteomics data (TMT and Silac)
  #' low sig: 0.1 > p >= 0.05
  #' medium sig: 0.05 > p >= 0.01
  #' high sig: 0.01 > p
  #' @param data_dt data.table of protein expression
  #' @param geneCol_v column in data_dt identifying the gene
  #' @param lfcCol_v column in data_dt indicating the fold change results (logFC for TMT; Z-score for Silac)
  #' @param pvalCol_v column in data_dt indicating the significance value (FDR for TMT; BH_adjusted for Silac)
  #' @param lfc_v value for cut-off for fold change. Will be included in plot as vertical dotted line
  #' @param ident1_v name of the reference treatment. Genes on right of plot are higher in this treatment than the other
  #' @param labelAll_v logical indication 
  #' @param labelTop_v optional, either NULL or a number indicating top N genes to label.
  #' @param labelDir_v can be "up", "down", or "both", indicating with side(s) of the volcano to label.
  #' @param labelSize_v pch size
  #' @param title_v what to name the plot
  #' @details Fill this in later
  #' @return volcano ggplot
  #' @export
  #' 
  
  ### Set local variables
  dirVarCheck_lsv <- list("up" = c("up", "Up", "u", "U"), 
                          "down" = c("down", "Down", "D", "d"), 
                          "both" = c("both", "Both", "B", "b"))
  multiGroupVolcanoColors_lsv <- list("NO" = "grey", 
                                      "DOWN" = c("downHigh" = "darkblue", "downMed" = "blue", "downLow" = "lightblue"),
                                      "UP" = c("upHigh" = "darkred", "upMed" = "red", "upLow" = "#FFB09C"))
  
  ### Add a column indicating up/dn direction
  if (!("diffExp" %in% colnames(data_dt))) {
    data_dt <- classifyProteo(data_dt = data_dt, lfcCol_v = lfcCol_v, pvalCol_v = pvalCol_v,
                              lfc_v = lfc_v, pval_v = c(0.01, 0.05), newName_v = "diffExp")
  }
  
  ### Get only significant ones and sort
  sig_lsdt <- list("UP" = data_dt[diffExp %in% names(multiGroupVolcanoColors_lsv$UP),][order(get(pvalCol_v))],
                   "DOWN" = data_dt[diffExp %in% names(multiGroupVolcanoColors_lsv$DOWN),][order(get(pvalCol_v))])
  
  ### Determine genes to label
  if (labelAll_v) {
    
    labelGenes_lsv <- lapply(sig_lsdt, function(x) x[[geneCol_v]])
    
  } else {
    
    if (!is.null(labelTop_v)) {
      
      labelGenes_lsv <- lapply(sig_lsdt, function(x) {
        y <- x[1:min(x[,.N], labelTop_v), get(geneCol_v)]
        return(y)
      })
      
    } # fi
    
  } # fi
  
  ### Add labels
  data_dt$DElabel <- ""
  
  if (labelDir_v %in% dirVarCheck_lsv$up | labelDir_v %in% dirVarCheck_lsv$both) {
    data_dt[get(geneCol_v) %in% setdiff(labelGenes_lsv$UP, "na"), DElabel := get(geneCol_v)]
  } # fi up
  
  if (labelDir_v %in% dirVarCheck_lsv$dn | labelDir_v %in% dirVarCheck_lsv$both) {
    data_dt[get(geneCol_v) %in% setdiff(labelGenes_lsv$DOWN, "na"), DElabel := get(geneCol_v)]
  } # fi dn
  
  plotColors_v <- unlist(multiGroupVolcanoColors_lsv)
  names(plotColors_v) <- gsub("DOWN\\.|UP\\.", "", names(plotColors_v))
  
  ### Make plot
  plot_gg <- ggplot2::ggplot(data = data_dt, aes(x = !!sym(lfcCol_v), 
                                                 y = -log10(!!sym(pvalCol_v)), 
                                                 col = diffExp,
                                                 label = DElabel)) +
    geom_point() + my_theme() +
    theme(legend.position = "right") +
    ggrepel::geom_label_repel(size = labelSize_v, show.legend = F, max.overlaps = Inf, label.padding = 0.1) +
    geom_vline(xintercept=c(-lfc_v, lfc_v), col="black", linetype = "dashed") +
    geom_hline(yintercept=-log10(0.01), col="black", linetype = "dashed") +
    geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
    scale_color_manual(values=plotColors_v, labels = names(plotColors_v)) +
    guides(color = guide_legend(title = paste0(ident1_v, " Dir"))) +
    ggtitle(title_v)
  
  ### Return
  return(plot_gg)
  
} # proteoVolcano

