### This has not been changed from the datatypes version yet
comparePhenotypesProteo <- function(data1_dt, data2_dt = NULL, names_v = NULL, geneCol_v = "Gene", col_v = "diffExp",
                          levels_lsv = list("No" = "NO",
                                            "Up" = c("upHigh", "upMed", "upLow"),
                                            "Down" = c("downHigh", "downMed", "downLow"))) {
  #' Compare Proteo
  #' @description
    #' Compare the significance groups of proteomics results 
    #' between two different experiments from the same data type
  #' @param data1_dt either a list of two data.tables, or a single data.table
  #' @param data2_dt if data1_dt is a list, this will be ignored. If it's a data.table, this must be a data.table
  #' @param names_v optional vector containing the experiment names for data1 and data2. Only used if data1_dt is NOT a list.
  #' @param geneCol_v column name in both data.tables that refers to the genes
  #' @param col_v column name in both data.tables that refers to differential expression "levels", from classifyProteo
  #' @param levels_lsv levels created from classifyProteo. This shouldn't change, but could theoretically.
  #' @details
    #' Want to get an idea of what to expect in each experiment. How many genes significant/shared/etc.
  #' @return list
  #' @export
  
  ### Wrangle different input data types
  if (!is.logical(all.equal(class(data1_dt), "list"))) {
    data_lsdt <- list(data1_dt, data2_dt)
    names(data_lsdt) <- names_v
  } else {
    data_lsdt <- data1_dt
  } # fi
  
  ### Get count/pct summary of up/dn genes in each level
  summary_lsdt <- summarizeProteo(data_lsdt = data_lsdt, geneCol_v = geneCol_v, col_v = col_v,
                                     levels_lsv = levels_lsv)
  
  summary_dt <- merge(summary_lsdt[[1]], summary_lsdt[[2]], by = "Category",
                      sort = F, suffixes = paste0("_", names(summary_lsdt)))
  
  
  ### Extract the actual genes in each up/dn level for each data type.
  geneGroups_lslslsv <- sapply(data_lsdt, function(x) getGeneGroups(data_dt = x, levels_lsv = levels_lsv, removeNA_v = T), simplify = F, USE.NAMES = T)
  
  ### Get comparisons
  groupNames_dt <- data.table("value" = listNames(geneGroups_lslslsv[[1]]))
  groupNames_dt$Dir <- sapply(groupNames_dt$value, function(x) strsplit(x, split = "_")[[1]][2])
  groupNames_dt$Level <- sapply(groupNames_dt$value, function(x) strsplit(x, split = "_")[[1]][3])
  groupNames_dt$value <- NULL
  groupNames_dt$genes1 <- numeric()
  groupNames_dt$genes2 <- numeric()
  groupNames_dt$sharedGenes <- numeric()
  
  for (i in 1:groupNames_dt[,.N]) {
    
    ### Get info
    currDir_v <- groupNames_dt$Dir[i]
    currLvl_v <- groupNames_dt$Level[i]
    
    ### Get genes
    curr1_v <- geneGroups_lslslsv[[1]][[currDir_v]][[currLvl_v]]
    curr2_v <- geneGroups_lslslsv[[2]][[currDir_v]][[currLvl_v]]
    
    ### Compare
    currSharedGenes_v <- intersect(curr1_v, curr2_v)
    curr1Genes_v <- setdiff(curr1_v, curr2_v)
    curr2Genes_v <- setdiff(curr2_v, curr1_v)
    
    groupNames_dt[i, sharedGenes := length(currSharedGenes_v)]
    groupNames_dt[i, genes1 := length(curr1Genes_v)]
    groupNames_dt[i, genes2 := length(curr1Genes_v)]
    
  } # for i
  
  ### Add percents
  groupNames_dt[, pct1Shared := round(sharedGenes / (sharedGenes + genes1) * 100, digits = 3)]
  groupNames_dt[, pct2Shared := round(sharedGenes / (sharedGenes + genes2) * 100, digits = 3)]
  groupNames_dt[, pct1Unique := round(genes1 / (sharedGenes + genes1) * 100, digits = 3)]
  groupNames_dt[, pct2Unique := round(genes2 / (sharedGenes + genes2) * 100, digits = 3)]
  
  colnames(groupNames_dt) <- gsub("1", names(summary_lsdt)[1], colnames(groupNames_dt))
  colnames(groupNames_dt) <- gsub("2", names(summary_lsdt)[2], colnames(groupNames_dt))
  
  out_lsdt <- list("indSummary" = summary_dt,
                   "compareSummary" = groupNames_dt,
                   "geneGroups" = geneGroups_lslslsv)
  
} # compareProteo
