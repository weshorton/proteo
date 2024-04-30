classifyProteo <- function(data_dt, lfcCol_v = "logFC", pvalCol_v = "FDR", lfc_v = 0.5, pval_v = c(0.01, 0.05), newName_v = "diffExp") {
  #' Classify Proteomics Data
  #' @description
  #' Group data_dt into 3 levels each for up/down expression:
  #' low sig: 0.1 > p >= 0.05
  #' medium sig: 0.05 > p >= 0.01
  #' high sig: 0.01 > p
  #' @param data_dt proteomics data.table with differential expression results
  #' @param lfcCol_v column in data_dt indicating the fold change results (logFC for TMT; Z-score for Silac)
  #' @param pvalCol_v column in data_dt indicating the significance value (FDR for TMT; BH_adjusted for Silac)
  #' @param lfc_v value for cut-off for fold change. Will be included in plot as vertical dotted line
  #' @param pval_v vector of length 2! first value indicates the threshold for 'high' significance, second for 'low'
  #' @param newName_v name of new column
  #' @import data.table
  #' @details
  #' Add a new column in data_dt summarizing the direction and level of each entry.
  #' @return data_dt with one new column
  #' @export
  
  ### Low
  data_dt[get(lfcCol_v) > lfc_v & get(pvalCol_v) >= pval_v[2] & get(pvalCol_v) < 0.1, (newName_v) := "upLow"]
  data_dt[get(lfcCol_v) < -lfc_v & get(pvalCol_v) >= pval_v[2] & get(pvalCol_v) < 0.1, (newName_v) := "downLow"]
  
  ### Med
  data_dt[get(lfcCol_v) > lfc_v & get(pvalCol_v) >= pval_v[1] & get(pvalCol_v) < pval_v[2], (newName_v) := "upMed"]
  data_dt[get(lfcCol_v) < -lfc_v & get(pvalCol_v) >= pval_v[1] & get(pvalCol_v) < pval_v[2], (newName_v) := "downMed"]
  
  ### High
  data_dt[get(lfcCol_v) > lfc_v & get(pvalCol_v) >= pval_v[1] & get(pvalCol_v) < pval_v[2], (newName_v) := "upHigh"]
  data_dt[get(lfcCol_v) < -lfc_v & get(pvalCol_v) >= pval_v[1] & get(pvalCol_v) < pval_v[2], (newName_v) := "downHigh"]
  
  ### Others
  data_dt[is.na(get(newName_v)), (newName_v) := "NO"]
  
  ### Factorize
  data_dt[[newName_v]] <- factor(data_dt[[newName_v]],
                            levels = c("NO", 
                                       "downLow", "downMed", "downHigh", 
                                       "upHigh", "upMed", "upLow"))
  
  ### Return
  return(data_dt)
  
}