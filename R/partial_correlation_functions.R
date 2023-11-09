#### This is a script for the partial correlation analysis and visualization helper functions.


## multiMediators is a function for identifying cytokine-metabolite relationships with enrichments for genes in multiple gene sets. For exmple, there is evidence that genes in interferron gamma signaling and heme metabolism gene sets are mediators of IL-1RA-Glutamate.
## Input: matrix of NES, gene pathways X cytokine-met
##        vector of pathways to test          
## Output: printed list of cytokine-metabolite relationships and the total NES scores, variance of NES scores, statistic (total * variance)

multiMediators <- function(input, pathways) {
 
  input <- input %>%
    dplyr::filter(rownames(input) %in% pathways)
  input <- as.data.frame(t(input))
  input$total <- apply(input,1,sum)
  input$variance <- apply(input,1,var)
  input <- input %>%
    mutate(stat = total*(variance))
  input <- input %>%
    arrange(-stat) %>%
    select(total,variance,stat)
  
  print(head(input,20))
} 

## gseaFinder is a helper function for filtering the gsea results.
## Inputs: character cyt-met example
##         charcter key pathway to search
## Outputs: vector of adjusted p value, NES, leading edge genes

gseaFinder <- function(cyt_met_example, key_pathway) {

  cyt_met_labels <- names(gsea_results_formatted[[1]]) 
  cyt_met_labels <- cyt_met_labels[grepl("NES", cyt_met_labels)]
  cyt_met_labels <- gsub(" NES", "", cyt_met_labels)
  path_id <- which(cyt_met_example ==cyt_met_labels)
  
  tmp_gsea <- lapply(gsea_results, function(x)
    x[[path_id]])
  
  res <- lapply(tmp_gsea, function(x) {
    x %>% filter(pathway == key_pathway)
  })
  
 final_res <- bind_rows(res) %>% dplyr::select(c(padj,NES,leadingEdge)) %>% as.data.frame()
  final_res <- final_res[1,]
return(final_res)
  
}

tmp <- gseaFinder("IFN-gamma-kynurenine", "HALLMARK_INTERFERON_GAMMA_RESPONSE")
tmp$leadingEdge
gseaFinder("IFN-gamma-kynurenine", "HALLMARK_HEME_METABOLISM")
gseaFinder("IL-1RA-Glutamate", "HALLMARK_INTERFERON_GAMMA_RESPONSE")
tmp3 <- gseaFinder("IL-1RA-Glutamate", "HALLMARK_HEME_METABOLISM")
tmp3$leadingEdge

tmp2 <- gseaFinder("IFN-gamma-acyl-C14:1 (Tetradecenoyl Carnitine)", "HALLMARK_INTERFERON_GAMMA_RESPONSE")
tmp2$leadingEdge

gseaFinder("IL-29-Docosapentaenoic acid", "HALLMARK_HEME_METABOLISM")



