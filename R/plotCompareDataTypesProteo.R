plotCompareDataTypesProteo <- function(tmt_dt, silac_dt, name_v,
                                   tmtLFCcol_v = "logFC", tmtPvalCol_v = "FDR",
                                   silacLFCcol_v = "Z-score", silacPvalCol_v = "BH_adjusted", silacVolcCol_v = "Log2(B/A)",
                                   lfc_v = 0.5, pval_v = c(0.01, 0.05)) {
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
  #' @return volcano plot, venn diagrams, scatterplot
  #' @export
  
  ### Classify genes
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
  
  ### Make volcano plot
  tmtVolcano_gg <- proteoVolcano(data_dt = tmt_dt, lfcCol_v = tmtLFCcol_v, pvalCol_v = tmtPvalCol_v,
                                 lfc_v = 1, ident1_v = name_v,
                                 labelAll_v = F, labelTop_v = 25, labelDir_v = "both",
                                 labelSize_v = 3, title_v = paste0("Genes expressed in ", name_v, " \nTMT"))
  
  silacVolcano_gg <- proteoVolcano(data_dt = silac_dt, lfcCol_v = silacVolcCol_v, pvalCol_v = silacPvalCol_v,
                                   lfc_v = 1, ident1_v = name_v,
                                   labelAll_v = F, labelTop_v = 25, labelDir_v = "both",
                                   labelSize_v = 3, title_v = paste0("Genes expressed in ", name_v, " \nSilac"))
  
  currVolcano_gg <- ggpubr::ggarrange(plotlist = list("tmt" = tmtVolcano_gg, "silac" = silacVolcano_gg), ncol = 2)
  print(currVolcano_gg)
  
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
  currVenn <- gridExtra::grid.arrange(currUpAll_venn, currDnAll_venn, currUpHigh_venn, currDnHigh_venn,
                                      nrow = 2)
  
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
  
  ggplot(comp_dt, aes(x = !!sym(tmtLFCcol_v), y = !!sym(silacLFCcol_v), color = diffExp_tmt)) +
    geom_point() + my_theme() +
    scale_color_manual(values = colors_v, breaks = names(colors_v)) +
    geom_smooth(data = sig_dt, method = 'lm', formula = y~x) +
    ggpubr::stat_cor(data = sig_dt, show.legend = F) +
    ggtitle(paste0("Compare TMT (x) and Silac (y) fold changes\n", name_v, " R value excludes 'NO'"))
    
} # compareDataTypesProteo

