plotComparePhenotypesProteo <- function(data1_dt, data2_dt = NULL, names_v = NULL,
                                   lfcCol1_v = "logFC", pValCol1_v = "FDR",
                                   lfcCol2_v = "logFC", pValCol2_v = "FDR",
                                   lfc_v = 0.5, pval_v = c(0.01, 0.05)) {
  #' Plot Compare Phenotypes Proteo
  #' @description
  #' Make various plots to compare the genes returned by two phenotypes
  #' @param data1_dt either a list of two data.tables, or a single data.table
  #' @param data2_dt if data1_dt is a list, this will be ignored. If it's a data.table, this must be a data.table
  #' @param names_v optional vector containing the experiment names for data1 and data2. Only used if data1_dt is NOT a list.
  #' @param lfcCol1_v column name in first data.table that corresponds to log fold change
  #' @param pValCol1_v column name in first data.table that corresponds to adjusted p-value
  #' @param lfcCol2_v column name in second data.table that corresponds to log fold change
  #' @param pValCol2_v column name in second data.table that corresponds to adjusted p-value
  #' @param lfc_v log-fold change cut-off value
  #' @param pval_v p-value cut-offs for grouping.
  #' @return volcano plot, venn diagrams, scatterplot
  #' @export
  
  ### Wrangle different input data types
  if (!is.logical(all.equal(class(data1_dt), "list"))) {
    data_lsdt <- list(data1_dt, data2_dt)
    names(data_lsdt) <- names_v
  } else {
    data_lsdt <- data1_dt
  } # fi
  
  ### Classify genes, if not already
  data_lsdt <- sapply(data_lsdt, function(x) {
    if (!"diffExp" %in% colnames(x)) {
      return(classifyProteo(x, lfcCol_v = lfcCol1_v, pvalCol_v = pValCol1_v,
                            lfc_v = lfc_v, pval_v = pval_v, newName_v = "diffExp"))
    } else {
      return(x)
    }
  }, simplify = F, USE.NAMES = T)
    
  ### Remove NA genes
  data_lsdt <- sapply(data_lsdt, function(x) x[Gene != "na",], simplify = F, USE.NAMES = T)
  
  ### Make volcano plot
  volcano1_gg <- proteoVolcano(data_dt = data_lsdt[[1]], lfcCol_v = lfcCol1_v, pvalCol_v = pValCol1_v,
                                 lfc_v = 1, ident1_v = names(data_lsdt)[1],
                                 labelAll_v = F, labelTop_v = 25, labelDir_v = "both",
                                 labelSize_v = 3, title_v = paste0("Genes expressed in ", names(data_lsdt)[1]))
  
  volcano2_gg <- proteoVolcano(data_dt = data_lsdt[[2]], lfcCol_v = lfcCol2_v, pvalCol_v = pValCol2_v,
                                   lfc_v = 1, ident1_v = names(data_lsdt)[1],
                                   labelAll_v = F, labelTop_v = 25, labelDir_v = "both",
                                   labelSize_v = 3, title_v = paste0("Genes expressed in ", names(data_lsdt)[1]))
  
  currVolcano_gg <- ggpubr::ggarrange(plotlist = list(volcano1_gg, volcano2_gg), ncol = 2)
  print(currVolcano_gg)
  
  ### Get up and down genes - all
  currUpGenesAll_lsv <- list("one" = data_lsdt[[1]][diffExp %in% c("upHigh", "upMed", "upLow"),Gene], 
                             "two" = data_lsdt[[2]][diffExp %in% c("upHigh", "upMed", "upLow"), Gene])
  currDnGenesAll_lsv <- list("one" = data_lsdt[[1]][diffExp %in% c("downHigh", "downMed", "downLow"),Gene],
                             "two" = data_lsdt[[2]][diffExp %in% c("downHigh", "downMed", "downLow"), Gene])
  names(currUpGenesAll_lsv) <- names(data_lsdt)
  names(currDnGenesAll_lsv) <- names(data_lsdt)
  
  ### Get up and down genes - high probability
  currUpGenesHigh_lsv <- list("one" = data_lsdt[[1]][diffExp == "upHigh",Gene], 
                              "two" = data_lsdt[[2]][diffExp == "upHigh", Gene])
  currDnGenesHigh_lsv <- list("one" = data_lsdt[[1]][diffExp == "downHigh",Gene],
                              "two" = data_lsdt[[2]][diffExp == "downHigh", Gene])
  
  ### All gene venn diagram
  currUpAll_venn <- invisible(venn.diagram(x = currUpGenesAll_lsv, filename = NULL, print.mode = c("raw", "percent"),
                                 main = paste0("Shared Up-reg genes All\n", paste(names(data_lsdt), collapse = " vs. "))))
  
  currDnAll_venn <- invisible(venn.diagram(x = currDnGenesAll_lsv, filename = NULL, print.mode = c("raw", "percent"),
                                 main = paste0("Shared Down-reg genes All\n", paste(names(data_lsdt), collapse = " vs. "))))
  
  ### High-Confidence gene venn diagram
  currUpHigh_venn <- invisible(venn.diagram(x = currUpGenesHigh_lsv, filename = NULL, print.mode = c("raw", "percent"),
                                  main = paste0("Shared Up-reg genes High-Conf\n", paste(names(data_lsdt), collapse = " vs. "))))
  
  currDnHigh_venn <- invisible(venn.diagram(x = currDnGenesHigh_lsv, filename = NULL, print.mode = c("raw", "percent"),
                                  main = paste0("Shared Down-reg genes High-Conf\n", paste(names(data_lsdt), collapse = " vs. "))))
  
  ### Arrange
  currVenn <- gridExtra::grid.arrange(currUpAll_venn, currDnAll_venn, currUpHigh_venn, currDnHigh_venn,
                                      nrow = 2)
  
  ### Compare LFC of shared genes
  sharedGenes_v <- c(intersect(currUpGenesAll_lsv[[1]], currUpGenesAll_lsv[[2]]),
                     intersect(currDnGenesAll_lsv[[1]], currDnGenesAll_lsv[[2]]))
  
  comp_dt <- merge(data_lsdt[[1]][,mget(c("Gene", lfcCol1_v, "diffExp"))],
                   data_lsdt[[2]][,mget(c("Gene", lfcCol2_v, "diffExp"))],
                   by = "Gene", all = F, sort = F, suffixes = paste0("_", names(data_lsdt)))
  
  tmpCol_v <- paste0("diffExp_", names(data_lsdt)[1])
  up_dt <- comp_dt[grep("up", get(tmpCol_v)),]; 
  up_dt[[tmpCol_v]] <- gsub("High|Med|Low", "", as.character(up_dt[[tmpCol_v]]))
  dn_dt <- comp_dt[grep("down", get(tmpCol_v)),]
  dn_dt[[tmpCol_v]] <- gsub("High|Med|Low", "", as.character(dn_dt[[tmpCol_v]]))
  
  sig_dt <- rbind(up_dt, dn_dt)
  sig_dt[[tmpCol_v]] <- "foo"
  
  colors_v <- c("NO" = "grey", 
                "downHigh" = "darkblue", "downMed" = "blue", "downLow" = "lightblue",
                "upHigh" = "darkred", "upMed" = "red", "upLow" = "#FFB09C")
  comp_dt[[tmpCol_v]] <- factor(as.character(comp_dt[[tmpCol_v]]), levels = names(colors_v))
  xCol_v <- paste0(lfcCol1_v, "_", names(data_lsdt)[1])
  yCol_v <- paste0(lfcCol2_v, "_", names(data_lsdt)[2])
  ggplot(comp_dt, aes(x = !!sym(xCol_v), y = !!sym(yCol_v), color = !!sym(tmpCol_v))) +
    geom_point() + my_theme() +
    scale_color_manual(values = colors_v, breaks = names(colors_v)) +
    geom_smooth(data = sig_dt, method = 'lm', formula = y~x) +
    ggpubr::stat_cor(data = sig_dt, show.legend = F) +
    ggtitle(paste0("Compare ", names(data_lsdt)[1], " (x) and ", 
                   names(data_lsdt)[2],"  (y) fold changes\n", " R value excludes 'NO'"))
    
} # compareDataTypesProteo

