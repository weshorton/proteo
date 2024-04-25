summarizeProteo <- function(data_lsdt, geneCol_v = "Gene", col_v = "diffExp",
                            levels_lsv = list("No" = "NO",
                                              "Up" = c("upHigh", "upMed", "upLow"),
                                              "Down" = c("downHigh", "downMed", "downLow"))) {
  #' Summarize Proteo
  #' @description summarize proteomics results
  #' @param data_lsdt list of tmt or silac datasets to summarize
  #' @param geneCol_v column name in the data.tables that indicates the gene
  #' @param col_v column to use to summarize
  #' @param levels_lsv the up and down levels that are in col_v
  #' @details Summarize number of up/down genes by confidence group.
  #' Both input lists need to have been run through classifyProteo() prior.
  #' @return summary table
  #' @export
  
  ### Prepare output list
  out_lsdt <- list()
  
  ### Populate
  for (i in 1:length(data_lsdt)) {
    
    ### Get info
    currName_v <- names(data_lsdt)[i]
    currData_dt <- data_lsdt[[currName_v]]
    
    ### Get gene groups
    # geneGroups_lslsv <- sapply(levels_lsv, function(x) {
    #   sapply(x, function(y) setdiff(currData_dt[get(col_v) == y, get(geneCol_v)], "na"), simplify = F, USE.NAMES = T)
    # }, simplify = F, USE.NAMES = T)
    geneGroups_lslsv <- getGeneGroups(currData_dt, levels_lsv = levels_lsv, removeNA_v = T)
    
    ### Summarize up/dn genes
    summary_dt <- parseProteoList(geneGroups_lslsv)
    
    ### Add
    out_lsdt[[currName_v]] <- summary_dt
  } # for i
  
  ### Return
  return(out_lsdt)
  
} # summarize Proteo

