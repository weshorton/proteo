parseProteoList <- function(list_lslsv) {
  #' Parse proteo list
  #' @description get summary stats
  #' @param list_lslsv gene groups of tmt or silac results
  #' @details list_lslsv must have 3 elements. Each is a list:
  #' 1. "No" - 'NO'
  #' 2. "Up" - 'upHigh', 'upMed', 'upLow'
  #' 3. "Down" - 'downHigh', 'downMed', 'downLow'
  #' @return summary stats
  #' @export
  
  ### Get totals
  N_nonSig_v <- length(list_lslsv$No$NO)
  N_up_v <- sapply(list_lslsv$Up, length)
  N_dn_v <- sapply(list_lslsv$Down, length)
  N_sig_v <- sum(N_up_v, N_dn_v)
  N_tot_v <- sum(N_nonSig_v, N_sig_v)
  
  ### Get percentages
  pctOfTotal_allSig_v <- N_sig_v / N_tot_v * 100
  pctOfTotal_upSig_v <- sum(N_up_v) / N_tot_v * 100
  pctOfTotal_dnSig_v <- sum(N_dn_v) / N_tot_v * 100
  
  pctOfSig_up_v <- sum(N_up_v) / N_sig_v * 100
  pctOfSig_dn_v <- sum(N_dn_v) / N_sig_v * 100
  
  pctOfSigUp_lvl_v <- N_up_v / sum(N_up_v) * 100
  pctOfSigDn_lvl_v <- N_dn_v / sum(N_dn_v) * 100
  
  ### Output
  out_dt <- data.table("Category" = c("Pct of Total Genes Sig", 
                                      "Pct of Total Genes Sig Up",
                                      "Pct of Total Genes Sig Down",
                                      "Pct of Sig Genes Up",
                                      "Pct of Sig Up Genes High Conf.",
                                      "Pct of Sig Up Genes Med Conf.",
                                      "Pct of Sig Up Genes Low Conf.",
                                      "Pct of Sig Genes Dn",
                                      "Pct of Sig Dn Genes High Conf.",
                                      "Pct of Sig Dn Genes Med Conf.",
                                      "Pct of Sig Dn Genes Low Conf."),
                       "Amount" = c(pctOfTotal_allSig_v,
                                    pctOfTotal_upSig_v,
                                    pctOfTotal_dnSig_v,
                                    pctOfSig_up_v,
                                    pctOfSigUp_lvl_v,
                                    pctOfSig_dn_v,
                                    pctOfSigDn_lvl_v))
  
  return(out_dt)
  
} # parseProteoList