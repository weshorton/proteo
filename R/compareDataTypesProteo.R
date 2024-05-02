compareDataTypesProteo <- function(tmt_lsdt, silac_lsdt, geneCol_v = "Gene", col_v = "diffExp",
                          levels_lsv = list("No" = "NO",
                                            "Up" = c("upHigh", "upMed", "upLow"),
                                            "Down" = c("downHigh", "downMed", "downLow"))) {
  #' Compare Proteo
  #' @description
    #' Compare the significance groups of proteomics results 
    #' between tmt and silac
  #' @param tmt_lsdt list of tmt results. Must be run through classifyProteo first.
  #' @param silac_lsdt list of silac results. Must be run through classifyProteo first.
  #' @param geneCol_v column name in both data.tables that refers to the genes
  #' @param col_v column name in both data.tables that refers to differential expression "levels", from classifyProteo
  #' @param levels_lsv levels created from classifyProteo. This shouldn't change, but could theoretically.
  #' @details
    #' Want to get an idea of what to expect in each experiment. How many genes significant/shared/etc.
  #' @return list
  #' @export
  
  ### Get count/pct summary of up/dn genes in each level
  tmtSummary_lsdt <- summarizeProteo(data_lsdt = tmt_lsdt, geneCol_v = geneCol_v, col_v = col_v,
                                     levels_lsv = levels_lsv)
  silacSummary_lsdt <- summarizeProteo(data_lsdt = silac_lsdt, geneCol_v = geneCol_v, col_v = col_v,
                                       levels_lsv = levels_lsv)
  summary_lsdt <- sapply(names(tmtSummary_lsdt), function(x) {
    merge(tmtSummary_lsdt[[x]], silacSummary_lsdt[[x]], by = "Category", sort = F, suffixes = c("_tmt", "_silac"))
  }, simplify = F, USE.NAMES = T)
  
  ### Extract the actual genes in each up/dn level for each data type.
  tmtGroups_lslslsv <- sapply(tmt_lsdt, function(x) getGeneGroups(data_dt = x, levels_lsv = levels_lsv, removeNA_v = T), simplify = F, USE.NAMES = T)
  silacGroups_lslslsv <- sapply(silac_lsdt, function(x) getGeneGroups(data_dt = x, levels_lsv = levels_lsv, removeNA_v = T), simplify = F, USE.NAMES = T)
  
  ### Get comparisons
  groupNames_mat <- listNames(tmtGroups_lslslsv)
  groupNames_dt <- suppressWarnings(as.data.table(melt(groupNames_mat)[,"value",drop=F]))
  groupNames_dt$Phenotype <- sapply(groupNames_dt$value, function(x) strsplit(x, split = "_")[[1]][2])
  groupNames_dt$Dir <- sapply(groupNames_dt$value, function(x) strsplit(x, split = "_")[[1]][3])
  groupNames_dt$Level <- sapply(groupNames_dt$value, function(x) strsplit(x, split = "_")[[1]][4])
  groupNames_dt$value <- NULL
  groupNames_dt$tmtGenes <- numeric()
  groupNames_dt$silacGenes <- numeric()
  groupNames_dt$sharedGenes <- numeric()
  
  for (i in 1:groupNames_dt[,.N]) {
    
    ### Get info
    currName_v <- groupNames_dt$Phenotype[i]
    currDir_v <- groupNames_dt$Dir[i]
    currLvl_v <- groupNames_dt$Level[i]
    
    ### Get genes
    currTMT_v <- tmtGroups_lslslsv[[currName_v]][[currDir_v]][[currLvl_v]]
    currSilac_v <- silacGroups_lslslsv[[currName_v]][[currDir_v]][[currLvl_v]]
    
    ### Compare
    currSharedGenes_v <- intersect(currTMT_v, currSilac_v)
    currTmtGenes_v <- setdiff(currTMT_v, currSilac_v)
    currSilacGenes_v <- setdiff(currSilac_v, currTMT_v)
    
    groupNames_dt[i, sharedGenes := length(currSharedGenes_v)]
    groupNames_dt[i, tmtGenes := length(currTmtGenes_v)]
    groupNames_dt[i, silacGenes := length(currSilacGenes_v)]
    
  } # for i
  
  ### Add percents
  groupNames_dt[, pctTmtShared := round(sharedGenes / (sharedGenes + tmtGenes) * 100, digits = 3)]
  groupNames_dt[, pctSilacShared := round(sharedGenes / (sharedGenes + silacGenes) * 100, digits = 3)]
  groupNames_dt[, pctTmtUnique := round(tmtGenes / (sharedGenes + tmtGenes) * 100, digits = 3)]
  groupNames_dt[, pctSilacUnique := round(silacGenes / (sharedGenes + silacGenes) * 100, digits = 3)]
  
  ### split to list
  compareSummary_lsdt <- sapply(names(summary_lsdt), function(x) {
    return(groupNames_dt[Phenotype == x,])}, simplify = F, USE.NAMES = T)
  
  out_lsdt <- list("indSummary" = summary_lsdt,
                   "compareSummary" = compareSummary_lsdt,
                   "geneGroups" = list("tmt" = tmtGroups_lslslsv,
                                       "silac" = silacGroups_lslslsv))
  
} # compareProteo