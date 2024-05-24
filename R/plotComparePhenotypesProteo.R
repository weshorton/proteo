plotComparePhenotypesProteo <- function(data1_dt, data2_dt = NULL, names_v = NULL,
                                   lfcCol1_v = "logFC", pValCol1_v = "FDR",
                                   lfcCol2_v = "logFC", pValCol2_v = "FDR",
                                   lfc_v = 0.5, pval_v = c(0.01, 0.05),
                                   geneLists_lslsv, print_v = F, outDir_v = NULL) {
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
  #' @param geneLists_lslsv list of lists (HALLMARK, GO), for GSEA analysis
  #' @param print_v logical indicating whether or not to print output to console.
  #' @param outDir_v optional directory to write output to.
  #' @return volcano plot, venn diagrams, scatterplot, pathway
  #' @export
  
  ###
  ### DATA WRANGLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  out_ls <- list()
  
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
  
  ###
  ### VOLCANO PLOT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Make volcano plot
  volcano1_gg <- proteoVolcano(data_dt = data_lsdt[[1]], lfcCol_v = lfcCol1_v, pvalCol_v = pValCol1_v,
                                 lfc_v = 1, ident1_v = names(data_lsdt)[1],
                                 labelAll_v = F, labelTop_v = 25, labelDir_v = "both",
                                 labelSize_v = 3, title_v = paste0("Genes expressed in ", names(data_lsdt)[1]))
  
  volcano2_gg <- proteoVolcano(data_dt = data_lsdt[[2]], lfcCol_v = lfcCol2_v, pvalCol_v = pValCol2_v,
                                   lfc_v = 1, ident1_v = names(data_lsdt)[1],
                                   labelAll_v = F, labelTop_v = 25, labelDir_v = "both",
                                   labelSize_v = 3, title_v = paste0("Genes expressed in ", names(data_lsdt)[2]))
  
  ### Combine
  currVolcano_gg <- ggpubr::ggarrange(plotlist = list(volcano1_gg, volcano2_gg), ncol = 2)
  
  ### Output
  out_ls[["volcano"]] <- currVolcano_gg
  if (print_v) print(currVolcano_gg)
  
  if (!is.null(outDir_v)) {
    pdf(file = file.path(outDir_v, paste0("volcano_compare_", paste0(names(data_lsdt), collapse = "-_-"), ".pdf")),
        width = 18, height = 10)
    print(currVolcano_gg)
    dev.off()
  } # fi !is.null(outDir_v)
  
  ###
  ### VENN DIAGRAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
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
  venn_gg <- ggpubr::ggarrange(plotlist = list(currUpAll_venn, currDnAll_venn, currUpHigh_venn, currDnHigh_venn), nrow = 2)
  
  ### Output
  out_ls[["venn"]] <- venn_gg
  if (print_v) print(venn_gg)
  
  if (!is.null(outDir_v)) {
    pdf(file = file.path(outDir_v, paste0("venn_compare_", paste0(names(data_lsdt), collapse = "-_-"), ".pdf")),
        width = 12, height = 12)
    print(venn_gg)
    dev.off()
  } # fi !is.null(outDir_v)
  
  ###
  ### SCATTER PLOT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
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
  scatter_gg <- ggplot(comp_dt, aes(x = !!sym(xCol_v), y = !!sym(yCol_v), color = !!sym(tmpCol_v))) +
    geom_point() + my_theme() + coord_equal() +
    scale_color_manual(values = colors_v, breaks = names(colors_v)) +
    geom_smooth(data = sig_dt, method = 'lm', formula = y~x) +
    ggpubr::stat_cor(data = sig_dt, show.legend = F) +
    ggtitle(paste0("Compare ", names(data_lsdt)[1], " (x) and ", 
                   names(data_lsdt)[2],"  (y) fold changes\n", " R value excludes 'NO'"))
  
  out_ls[['scatter']] <- scatter_gg
  if (print_v) print(scatter_gg)
  
  if (!is.null(outDir_v)) {
    pdf(file = file.path(outDir_v, paste0("scatterplot_compare_", paste0(names(data_lsdt), collapse = "-_-"), ".pdf")),
        width = 10, height = 10)
    print(scatter_gg)
    dev.off()
  } # fi !is.null(outDir_v)
  
  ###
  ### GSEA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  lfcCols_lsv <- list(lfcCol1_v, lfcCol2_v); names(lfcCols_lsv) <- names(data_lsdt)
  getCols_v <- c("pathway", "padj", "NES", "Dir")
  gseaFullOut_lslsls <- gseaCompareOut_lsls <- list()
  widths_v <- c(20, 20, 50); names(widths_v) <- names(geneLists_lslsv)
  
  for (i in 1:length(geneLists_lslsv)) {
    
    currListName_v <- names(geneLists_lslsv)[i]
    currList_lsv <- geneLists_lslsv[[i]]
    
    ### Waterfall plot
    gsea_lsls <- sapply(names(data_lsdt), function(x) {
      plotGSEA(data_dt = data_lsdt[[x]], list_lsv = currList_lsv, listName_v = currListName_v,
               title_v = paste0("Pathway Enrichment for ", x), name_v = x,
               lfcCol_v = lfcCols_lsv[[x]], le_v = F, shortenNames_v = T, 
               lfc_v = 0.5, pval_v = 0.05, labelSize_v = 1)
    }, simplify = F, USE.NAMES = T)
    
    plot_gg <- ggpubr::ggarrange(plotlist = list(gsea_lsls[[1]]$plot, gsea_lsls[[2]]$plot), ncol = 2)
    
    ### Scatter
    data_dt <- merge(gsea_lsls[[1]]$data[,mget(getCols_v)], gsea_lsls[[2]]$data[,mget(getCols_v)],
                     by = "pathway", sort = F, all = T, suffixes = paste0("_", names(gsea_lsls)))
    data_dt <- data_dt[,mget(c("pathway", grep("NES", colnames(data_dt), value = T),
                               grep("padj", colnames(data_dt), value = T), grep("Dir", colnames(data_dt), value = T)))]
    gseaScatter_gg <- ggplot(data = data_dt, aes(x = !!sym(paste0("NES_", names(data_lsdt)[1])), y = !!sym(paste0("NES_", names(data_lsdt)[2])))) +
      geom_point() + my_theme() + coord_equal() +
      geom_smooth(data = data_dt, method = 'lm', formula = y~x) +
      ggpubr::stat_cor(data = data_dt, show.legend = F) +
      ggtitle(paste0("Compare ", names(data_lsdt)[1], " (x) and ", names(data_lsdt)[2], " (y)\n", currListName_v, " NES"))
    
    ### Venn diagram
    gseaUpVennList_lsv <- list(gsea_lsls[[1]]$data[NES > 0, pathway], gsea_lsls[[2]]$data[NES > 0, pathway]); names(gseaUpVennList_lsv) <- names(gsea_lsls)
    gseaDnVennList_lsv <- list(gsea_lsls[[1]]$data[NES < 0, pathway], gsea_lsls[[2]]$data[NES < 0, pathway]); names(gseaDnVennList_lsv) <- names(gsea_lsls)
    
    if (length(gseaUpVennList_lsv[[1]]) > 0 & length(gseaUpVennList_lsv[[2]]) > 0) {
      gseaUp_venn <- invisible(venn.diagram(x = gseaUpVennList_lsv, filename = NULL, print.mode = c("raw", "percent"),
                                            main = paste0("Shared Up-reg ", currListName_v, " Pathways\n", names(gsea_lsls)[[1]], " vs. ", names(gsea_lsls)[[2]])))
    } else {
      gseaUp_venn <- NULL
    }
    
    if (length(gseaDnVennList_lsv[[1]]) > 0 & length(gseaDnVennList_lsv[[2]]) > 0) {
      gseaDn_venn <- invisible(venn.diagram(x = gseaDnVennList_lsv, filename = NULL, print.mode = c("raw", "percent"),
                                            main = paste0("Shared Down-reg ", currListName_v, " Pathways\n", names(gsea_lsls)[[1]], " vs. ", names(gsea_lsls)[[2]])))
    } else {
      gseaDn_venn <- NULL
    }
    
    ### Arrange
    gseaVenn_gg <- ggpubr::ggarrange(plotlist = list(gseaUp_venn, gseaDn_venn), nol = 2)
    
    ### Output
    if (print_v) {
      print(plot_gg)
      print(DT::datatable(data_dt))
      print(gseaScatter_gg)
      print(gseaVenn_gg)
    } # fi print_v
    
    if (!is.null(outDir_v)) {
      pdf(file = file.path(outDir_v, paste0(currListName_v, "_compare_", paste0(names(data_lsdt), collapse = "-_-"), "_gsea.pdf")),
          width = widths_v[currListName_v], height = max(sapply(gsea_lsls, function(x) nrow(x$data))))
      print(plot_gg)
      dev.off()
      
      pdf(file = file.path(outDir_v, paste0(currListName_v, "_compare_", paste0(names(data_lsdt), collapse = "-_-"), "_gseaScatter.pdf")))
      print(gseaScatter_gg)
      dev.off()
      
      pdf(file = file.path(outDir_v, paste0(currListName_v, "_compare_", paste0(names(data_lsdt), collapse = "-_-"), "_gseaVenn.pdf")))
      print(gseaVenn_gg)
      dev.off()
    } # fi !is.null(outDir_v)
    
    ### Add to full list
    gseaFullOut_lslsls[[currListName_v]] <- gsea_lsls
    gseaCompareOut_lsls[[currListName_v]] <- list("data" = data_dt, "waterfall" = plot_gg, "scatter" = gseaScatter_gg, "venn" = gseaVenn_gg)
    
  } # for i
  
  out_ls[['fullGSEA']] <- gseaFullOut_lslsls
  out_ls[['compareGSEA']] <- gseaCompareOut_lsls
  
  return(out_ls)
  
  
  
} # compareDataTypesProteo

