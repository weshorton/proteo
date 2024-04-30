classifyProteo <- function(data_dt, lfcCol_v, pvalCol_v, lfc_v = 0.5, pval_v = c(0.01, 0.05), newName_v = "diffExp") {
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
  #' @details
  #' Add a new column in data_dt summarizing the direction and level of each entry.
  #' @return data_dt with one new column
  #' @export
  
  ### Add a column indicating up/dn direction
  #data_dt[[newName_v]] <- "NO"
  setDT(data_dt)
  
  ### Low
  print("one")
  #print(data_dt[((lfcCol_v) > lfc_v & (pvalCol_v) >= pval_v[2] & (pvalCol_v) < 0.1), ])
  #data_dt[, (newName_v) := ifelse((.SD[[1]] > lfc_v & .SD[[2]] >= pval_v[2] & .SD[[2]] < 0.1), "upLow", get(newName_v))]
  #data_dt[(get(lfcCol_v) > lfc_v & get(pvalCol_v) >= pval_v[2] & get(pvalCol_v) < 0.1), (newName_v) := "upLow"]
  print(class(data_dt))
  print(head(data_dt))
  print(data_dt[get(lfcCol_v) > lfc_v,])
  data_dt[get(lfcCol_v) > lfc_v, newName_v := "upLow"]
  print("oneb")
  data_dt[(get(lfcCol_v) < -lfc_v & get(pvalCol_v) >= pval_v[2] & get(pvalCol_v) < 0.1), (newName_v) := "downLow"]
  
  ### Med
  print("two")
  data_dt[(get(lfcCol_v) > lfc_v & get(pvalCol_v) >= pval_v[1] & get(pvalCol_v) < pval_v[2]), (newName_v) := "upMed"]
  print("twob")
  data_dt[(get(lfcCol_v) < -lfc_v & get(pvalCol_v) >= pval_v[1] & get(pvalCol_v) < pval_v[2]), (newName_v) := "downMed"]
  
  ### High
  print("three")
  data_dt[(get(lfcCol_v) > lfc_v & get(pvalCol_v) < pval_v[1]), (newName_v) := "upHigh"]
  print("threeb")
  data_dt[(get(lfcCol_v) < -lfc_v & get(pvalCol_v) < pval_v[1]), (newName_v) := "downHigh"]
  
  ### Factorize
  data_dt$diffExp <- factor(data_dt$diffExp,
                            levels = c("NO", 
                                       "downLow", "downMed", "downHigh", 
                                       "upHigh", "upMed", "upLow"))
  
  ### Return
  return(data_dt)
  
}