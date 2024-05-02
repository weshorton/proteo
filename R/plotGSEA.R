plotGSEA <- function(data_dt, list_v, listName_v, title_v = NULL, name_v, 
                     lfcCol_v = "logFC", le_v = F, shortenNames_v = T, 
                     lfc_v = 0.5, pval_v = 0.05, labelSize_v = 1) {
  #' Run GSEA
  #' @description
    #' Run GSEA and output results. Optionally output plots
  #' @param data_dt DEG table to run GSEA on
  #' @param list_v list of genes (gene set) used for GSEA
  #' @param listName_v name of the gene set
  #' @param title_v optional plot title
  #' @param name_v name of experiment (e.g. M1E)
  #' @param pvalCol_v column in data_dt that has p-values
  #' @param le_v logical indicating if leading edge should be run. DOESN'T WORK RIGHT NOW.
  #' @param plot_v logical indicating if plots should be made/saved or just the results calculated.
  #' @param shortenNames_v logical indicating if gene set names should be shortened.
  #' @param lfc_v log fold change - passed through for leading edge gene filter
  #' @param pval_v adjusted p-value - used to filter significant pathways and also passed through for leading edge filter
  #' @param labelSize_v how big to make leading edge volcano label text. Default is 1.
  #' @return not sure yet
  #' @export
  
  ### Make sure descending log fold change
  data_dt <- data_dt[order(-get(lfcCol_v))]
  
  ### Run GSEA calculations
  gsea_res <- wrh.scRNA::runFGSEA(data_df = data_dt, rankCol_v = lfcCol_v, pathways_v = list_v, minSize_v = 15, maxSize_v = 500, seed_v = 42)
  
  ### Skip if empty
  if (nrow(gsea_res) == 0) {
    
    if (!is.null(spop_v)) {
      cat(sprintf("No pathways found for: %s - %s\n", name_v, listName_v))
    } else {
      cat(sprintf("No pathways found for: %s - %s\n", name_v, listName_v))
    } # fi spop_v
    
  } else {
    
    ### Remove GENESET NAME from pathway name
    sub_v <- paste0("^", listName_v, "_")
    gsea_res$pathway <- gsub(sub_v, "", gsea_res$pathway)
    gsea_res$Parent <- gsub(sub_v, "", gsea_res$Parent)
    
    ### Shorten names
    if (shortenNames_v) {
      for (n_v in 1:nrow(gsea_res)) {
        path_v <- strsplit(gsea_res$pathway[n_v], split = "_")[[1]]
        if (length(path_v) < 4) {
          path_v <- paste(path_v, collapse = "_")
        } else {
          path_v <- paste(c(path_v[1:3], sapply(path_v[4:length(path_v)], function(x) strtrim(x, 4))), collapse = "_")
        } # fi length
      } # for n_v
    } # fi shortenNames_v
    
    ### Filter and skip if empty
    gsea_res <- gsea_res[padj <= pval_v,]
    
    if (nrow(gsea_res) == 0) {
      
      if (!is.null(spop_v)) {
        cat(sprintf("None of the pathways found for: %s - %s - %s - %s - %s were significant.\n", pop_v, spop_v, treat_v, otherTreat_v, listName_v))
      } else {
        cat(sprintf("None of the pathways found for: %s - %s - %s - %s were significant.\n", pop_v, treat_v, otherTreat_v, listName_v))
      } # fi spop
      
    } else {
      
      ### Add Color column
      gsea_res$Dir <- "UP"; gsea_res[NES < 0, Dir := "DOWN"]
      gsea_res$Dir <- factor(gsea_res$Dir, levels = c("UP", "DOWN"))
      
      ### Make plot
      plot_gg <- ggplot2::ggplot(data = gsea_res, aes(x = reorder(pathway, NES), y = NES, fill = Dir)) +
        geom_col() + coord_flip() + my_theme() +
        scale_fill_manual(values = gseaColors_v, labels = names(gseaColors_v)) +
        labs(x = "Pathway", y = "NES")
      
      ### Add title
      if (!is.null(title_v)) {
        plot_gg <- plot_gg + ggtitle(title_v)
      } # fi
      
      ### Output
      out_ls <- list("data" = gsea_res, "plot" = plot_gg)
      
      ### Volcano plot of Leading Edge genes
      if (le_v) {
        
        path_v <- gsea_res$pathway
        leVolcano_lsgg <- sapply(path_v, function(x) {
          wrh.scRNA::leadingEdgeVolcano(data_dt = data_dt, gsea_res = gsea_res, pathway_v = x, treat_v = name_v, lfc_v = lfc_v, 
                             pval_v = pval_v, title_v = title_v, labelSize_v = labelSize_v)
        }, simplify = F)
        
        out_ls[["le"]] <- leVolcano_lsgg
        
      } # fi le_v
      
      return(out_ls)
      
    } # fi nrow(gsea_res) == 0 (pval filter)
    
  } # fi nrow(gsea_res) == 0 (orig run)
  
} # runGSEA

