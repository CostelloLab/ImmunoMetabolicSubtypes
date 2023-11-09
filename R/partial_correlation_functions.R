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

