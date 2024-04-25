getGeneGroups <- function(data_dt, levels_lsv, removeNA_v = T, col_v = "diffExp", geneCol_v = "Gene") {
  #' Get Gene Groups
  #' @description
  #' Get all the genes in different signficance groups
  #' @param data_dt tmt or silace results table
  #' @param levels_lsv the significance levels
  #' @param removeNA_v logical. Should 'na' genes be removed?
  #' @return list of lists
  #' @export
  
  groups_lslsv <- sapply(levels_lsv, function(x) {
    sapply(x, function(y) {
      genes_v <- data_dt[get(col_v) == y, get(geneCol_v)]
      if (removeNA_v) setdiff(genes_v, "na")
      return(genes_v)
    }, simplify = F, USE.NAMES = T)
  }, simplify = F, USE.NAMES = T)
  
  return(groups_lslsv)
} # getGeneGroups
