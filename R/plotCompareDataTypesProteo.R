plotCompareDataTypesProteo <- function(tmt_dt, silac_dt, name_v,
                                   tmtLFCcol_v = "logFC", tmtPvalCol_v = "FDR",
                                   silacLFCcol_v = "Z-score", silacPvalCol_v = "BH_adjusted", silacVolcCol_v = "Log2(B/A)",
                                   lfc_v = 0.5, pval_v = c(0.01, 0.05),
                                   geneLists_lslsv, print_v = F, outDir_v = NULL) {
  #' Plot Compare Data Types Proteo
  #' @description
    #' Make various plots to compare the genes returned by each data type
  #' @param tmt_dt data.table of tmt results
  #' @param silac_dt data.table of silac results
  #' @param name_v Name of the experiment (e.g. M1E)
  #' @param tmtLFCcol_v column name in tmt_dt that corresponds to log fold change
  #' @param tmtPvalCol_v column name in tmt_dt that corresponds to adjusted p-value
  #' @param silacLFCcol_v column name in silac_dt that corresponds to log fold change
  #' @param silacPvalCol_v column name in silac_dt that corresponds to adjusted p-value
  #' @param silacVolCol_v not sure if lfcol or this col should be plotted for silac
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
  
  ### Classify genes (if not)
  if (!"diffExp" %in% colnames(tmt_dt)) {
    tmt_dt <- classifyProteo(tmt_dt, lfcCol_v = tmtLFCcol_v, pvalCol_v = tmtPvalCol_v,
                             lfc_v = lfc_v, pval_v = pval_v, newName_v = "diffExp")
  }
  if (!"diffExp" %in% colnames(silac_dt)) {
    silac_dt <- classifyProteo(silac_dt, lfcCol_v = silacLFCcol_v, pvalCol_v = silacPvalCol_v,
                               lfc_v = lfc_v, pval_v = pval_v, newName_v = "diffExp")
  }
    
  ### Remove NA genes
  tmt_dt <- tmt_dt[Gene != "na",]
  silac_dt <- silac_dt[Gene != "na",]
  
  ###
  ### VOLCANO PLOT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Make volcano plot
  tmtVolcano_gg <- proteoVolcano(data_dt = tmt_dt, lfcCol_v = tmtLFCcol_v, pvalCol_v = tmtPvalCol_v,
                                 lfc_v = 1, ident1_v = name_v,
                                 labelAll_v = F, labelTop_v = 25, labelDir_v = "both",
                                 labelSize_v = 3, title_v = paste0("Genes expressed in ", name_v, " \nTMT"))
  
  silacVolcano_gg <- proteoVolcano(data_dt = silac_dt, lfcCol_v = silacVolcCol_v, pvalCol_v = silacPvalCol_v,
                                   lfc_v = 1, ident1_v = name_v,
                                   labelAll_v = F, labelTop_v = 25, labelDir_v = "both",
                                   labelSize_v = 3, title_v = paste0("Genes expressed in ", name_v, " \nSilac"))
  
  ### Combine
  currVolcano_gg <- ggpubr::ggarrange(plotlist = list("tmt" = tmtVolcano_gg, "silac" = silacVolcano_gg), ncol = 2)
  
  ### Output
  out_ls[["volcano"]] <- currVolcano_gg
  if (print_v) print(currVolcano_gg)
  
  if (!is.null(outDir_v)) {
    pdf(file = file.path(outDir_v, paste0("volcano_compareDataTypes_", name_v, ".pdf")),
        width = 18, height = 10)
    print(currVolcano_gg)
    dev.off()
  } # fi !is.null(outDir_v)
  
  ###
  ### VENN DIAGRAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Get up and down genes - all
  currUpGenesAll_lsv <- list("tmt" = tmt_dt[diffExp %in% c("upHigh", "upMed", "upLow"),Gene], 
                             "silac" = silac_dt[diffExp %in% c("upHigh", "upMed", "upLow"), Gene])
  currDnGenesAll_lsv <- list("tmt" = tmt_dt[diffExp %in% c("downHigh", "downMed", "downLow"),Gene],
                             "silac" = silac_dt[diffExp %in% c("downHigh", "downMed", "downLow"), Gene])
  
  ### Get up and down genes - high probability
  currUpGenesHigh_lsv <- list("tmt" = tmt_dt[diffExp == "upHigh",Gene], 
                              "silac" = silac_dt[diffExp == "upHigh", Gene])
  currDnGenesHigh_lsv <- list("tmt" = tmt_dt[diffExp == "downHigh",Gene],
                              "silac" = silac_dt[diffExp == "downHigh", Gene])
  
  ### All gene venn diagram
  currUpAll_venn <- invisible(venn.diagram(x = currUpGenesAll_lsv, filename = NULL, print.mode = c("raw", "percent"),
                                 main = paste0("Shared Up-reg genes All\n", name_v)))
  
  currDnAll_venn <- invisible(venn.diagram(x = currDnGenesAll_lsv, filename = NULL, print.mode = c("raw", "percent"),
                                 main = paste0("Shared Down-reg genes All\n", name_v)))
  
  ### High-Confidence gene venn diagram
  currUpHigh_venn <- invisible(venn.diagram(x = currUpGenesHigh_lsv, filename = NULL, print.mode = c("raw", "percent"),
                                  main = paste0("Shared Up-reg genes High-Conf\n", name_v)))
  
  currDnHigh_venn <- invisible(venn.diagram(x = currDnGenesHigh_lsv, filename = NULL, print.mode = c("raw", "percent"),
                                  main = paste0("Shared Down-reg genes High-Conf\n", name_v)))
  
  ### Arrange
  venn_gg <- ggpubr::ggarrange(plotlist = list(currUpAll_venn, currDnAll_venn, currUpHigh_venn, currDnHigh_venn), nrow = 2)
  
  ### Output
  out_ls[["venn"]] <- venn_gg
  if (print_v) print(venn_gg)
  
  if (!is.null(outDir_v)) {
    pdf(file = file.path(outDir_v, paste0("venn_compareDataTypes_", name_v, ".pdf")),
        width = 12, height = 12)
    print(venn_gg)
    dev.off()
  } # fi !is.null(outDir_v)
  
  ###
  ### SCATTER PLOT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Compare LFC of shared genes
  sharedGenes_v <- c(intersect(currUpGenesAll_lsv$tmt, currUpGenesAll_lsv$silac),
                     intersect(currDnGenesAll_lsv$tmt, currDnGenesAll_lsv$silac))
  
  comp_dt <- merge(tmt_dt[,mget(c("Gene", tmtLFCcol_v, "diffExp"))],
                   silac_dt[,mget(c("Gene", silacLFCcol_v, "diffExp"))],
                   by = "Gene", all = F, sort = F, suffixes = c("_tmt", "_silac"))
  
  up_dt <- comp_dt[grep("up", diffExp_tmt),]; 
  up_dt$diffExp_tmt <- gsub("High|Med|Low", "", as.character(up_dt$diffExp_tmt))
  dn_dt <- comp_dt[grep("down", diffExp_tmt),]
  dn_dt$diffExp_tmt <- gsub("High|Med|Low", "", as.character(dn_dt$diffExp_tmt))
  
  sig_dt <- rbind(up_dt, dn_dt)
  sig_dt$diffExp_tmt <- "foo"
  
  colors_v <- c("NO" = "grey", 
                "downHigh" = "darkblue", "downMed" = "blue", "downLow" = "lightblue",
                "upHigh" = "darkred", "upMed" = "red", "upLow" = "#FFB09C")
  comp_dt$diffExp_tmt <- factor(as.character(comp_dt$diffExp_tmt), levels = names(colors_v))
  
  ### Make plot
  scatter_gg <- ggplot(comp_dt, aes(x = !!sym(tmtLFCcol_v), y = !!sym(silacLFCcol_v), color = diffExp_tmt)) +
    geom_point() + my_theme() +
    scale_color_manual(values = colors_v, breaks = names(colors_v)) +
    geom_smooth(data = sig_dt, method = 'lm', formula = y~x) +
    ggpubr::stat_cor(data = sig_dt, show.legend = F) +
    ggtitle(paste0("Compare TMT (x) and Silac (y) fold changes\n", name_v, " R value excludes 'NO'"))
  
  ### Output
  out_ls[['scatter']] <- scatter_gg
  if (print_v) print(scatter_gg)
  
  if (!is.null(outDir_v)) {
    pdf(file = file.path(outDir_v, paste0("scatterplot_compareDataTypes_", name_v, ".pdf")),
        width = 10, height = 10)
    print(scatter_gg)
    dev.off()
  } # fi !is.null(outDir_v)
  
  ###
  ### GSEA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  lfcCols_lsv <- list(tmtLFCcol_v, silacLFCcol_v); names(lfcCols_lsv) <- c("TMT", "Silac")
  getCols_v <- c("pathway", "padj", "NES", "Dir")
  gseaFullOut_lslsls <- gseaCompareOut_lsls <- list()
  widths_v <- c(20, 20, 50); names(widths_v) <- names(geneLists_lslsv)
  
  data_lsdt <- list("TMT" = tmt_dt, "Silac" = silac_dt)
  
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
    
    ### Scatterplot
    data_dt <- merge(gsea_lsls[[1]]$data[,mget(getCols_v)], gsea_lsls[[2]]$data[,mget(getCols_v)],
                     by = "pathway", sort = F, all = T, suffixes = paste0("_", names(gsea_lsls)))
    data_dt <- data_dt[,mget(c("pathway", grep("NES", colnames(data_dt), value = T),
                               grep("padj", colnames(data_dt), value = T), grep("Dir", colnames(data_dt), value = T)))]
    gseaScatter_gg <- ggplot(data = data_dt, aes(x = NES_TMT, y = NES_Silac)) +
      geom_point() + my_theme() +
      geom_smooth(data = data_dt, method = 'lm', formula = y~x) +
      ggpubr::stat_cor(data = data_dt, show.legend = F) +
      ggtitle(paste0("Compare TMT (x) and Silac (y) - ", name_v, "\n", currListName_v, " NES"))
    
    ### Venn diagram
    gseaUpVennList_lsv <- list("TMT" = gsea_lsls$TMT$data[NES > 0, pathway], "Silac" = gsea_lsls$Silac$data[NES > 0, pathway])
    gseaDnVennList_lsv <- list("TMT" = gsea_lsls$TMT$data[NES < 0, pathway], "Silac" = gsea_lsls$Silac$data[NES < 0, pathway])
    
    if (length(gseaUpVennList_lsv[[1]]) > 0 & length(gseaUpVennList_lsv[[2]]) > 0) {
      gseaUp_venn <- invisible(venn.diagram(x = gseaUpVennList_lsv, filename = NULL, print.mode = c("raw", "percent"),
                               main = paste0("Shared Up-reg ", currListName_v, " Pathways - ", name_v)))
    } else {
      gseaUp_venn <- NULL
    }
    
    if (length(gseaDnVennList_lsv[[1]]) > 0 & length(gseaDnVennList_lsv[[2]]) > 0) {
      gseaDn_venn <- invisible(venn.diagram(x = gseaDnVennList_lsv, filename = NULL, print.mode = c("raw", "percent"),
                               main = paste0("Shared Down-reg ", currListName_v, " Pathways - ", name_v)))
    } else {
      gseaDn_venn <- NULL
    }
    
    ### Arrange
    gseaVenn_gg <- ggpubr::ggarrange(plotlist = list(gseaUp_venn, gseaDn_venn), nol = 2)
    
    if (print_v) {
      print(plot_gg)
      print(DT::datatable(data_dt))
      print(gseaScatter_gg)
      print(gseaVenn_gg)
    } # fi print_v
    
    if (!is.null(outDir_v)) {
      pdf(file = file.path(outDir_v, paste0(currListName_v, "_compareDataTypes_", name_v, "_gsea.pdf")),
          width = widths_v[currListName_v], height = max(sapply(gsea_lsls, function(x) nrow(x$data))))
      print(plot_gg)
      dev.off()
      
      pdf(file = file.path(outDir_v, paste0(currListName_v, "_compareDataTypes_", name_v, "_gseaScatter.pdf")))
      print(gseaScatter_gg)
      dev.off()
      
      pdf(file = file.path(outDir_v, paste0(currListName_v, "_compareDataTypes_", name_v, "_gseaVenn.pdf")))
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

