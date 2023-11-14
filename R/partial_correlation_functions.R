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







## Finding commonalities accross the topmost cyt-met relationships by pathway

# load the require libraries
library(ComplexHeatmap)
library(circlize)
library(fgsea)


### multi


multiHeatmap <- function(key_pathways, met_path, formatted_gsea) {
    
    met_path <- met_path %>%
        as.data.frame() %>%
        rownames_to_column("metabolite") %>%
        pivot_longer(cols = -metabolite, names_to = "pathway") %>%
        filter(value == 1)
    
  
  
  cyt_met_examples <- formatted_gsea %>%                                                  
    rownames_to_column("pathway") %>%
    filter(pathway %in% c("HALLMARK_INTERFERON_GAMMA_RESPONSE", 
                         "HALLMARK_HEME_METABOLISM",
                         "HALLMMARK_OXIDATIVE_PHOSPHORYLATION"))%>%
    pivot_longer(cols = contains("NES"), values_to = "NES", names_to = "cyt_met_NES") %>%
     pivot_longer(cols = contains("padj"), values_to = "padj", names_to = "cyt_met_padj") %>%
    filter(gsub(" NES", "", cyt_met_NES) == gsub(" padj", "", cyt_met_padj)) %>%
    mutate(cyt_met = gsub(" NES", "", cyt_met_NES)) %>%
    select(pathway, cyt_met, NES, padj) %>%
    filter(padj < .05 & NES > 0) %>%
    arrange(padj) %>%
    top_n(100) %>%
    distinct(cyt_met) %>%
    .$cyt_met
  
    #pathway_genes <- pathways_all[key_pathway][[1]]
  
  res <- lapply(cyt_met_examples, function(x) {
    tmp <- long_z[[1]] %>%
      filter(cyt_met == x) %>%
      arrange(desc(combined_Z)) %>%
      mutate(rank = row_number()) %>%
      select(gene,rank)
    names(tmp)[2] <- x
    return(tmp)
})

  
  rankings_df <- res %>% reduce(left_join, by = "gene") %>% column_to_rownames("gene")
  
  rankings_df <- as.matrix(rankings_df)
  
  tokeep <- apply(rankings_df, 1, function(x) sum(x < 100) > 2)
 toplot <- t(rankings_df[tokeep, ])
 
  
 IFN_gamma <-as.factor( ifelse(colnames(toplot) %in% pathways_all$HALLMARK_INTERFERON_GAMMA_RESPONSE, 1,0))
 names(IFN_gamma) <- colnames(toplot)
  Heme_Metabolism <-as.factor( ifelse(colnames(toplot) %in% pathways_all$HALLMARK_HEME_METABOLISM, 1,0))
 names(Heme_Metabolism) <- colnames(toplot)
  Oxidative_Phosphorylation <-as.factor( ifelse(colnames(toplot) %in% pathways_all$HALLMARK_OXIDATIVE_PHOSPHORYLATION, 1,0))
 names(Oxidative_Phosphorylation) <- colnames(toplot)
 
 

  


col_fun = colorRamp2(c(min(toplot), median(toplot), max(toplot) ), c("blue", "white", "gray"))

ha <- HeatmapAnnotation(IFN_gamma_response = IFN_gamma, 
                        Heme_Metabolism = Heme_Metabolism,
                        col = list(IFN_gamma_response = c("1"= "black", "0" = "white"),
                                   Heme_Metabolism =  c("1"= "red", "0" = "white")))

metabolites <-  full_results %>%
   filter(cyt_met %in% cyt_met_examples) %>%
   distinct(cyt_met, .keep_all = TRUE) %>%
   arrange(match(cyt_met, cyt_met_examples)) %>%
   .$metabolite


metabolites <- as.data.frame(metabolites)
names(metabolites) = "original"

met_pathways <- metabolites %>%
  left_join(met_path, by = c("original" = "metabolite"))

hr <- rowAnnotation(Metabolite_Class = met_pathways$pathway)
  

  p1 <- ComplexHeatmap::Heatmap(toplot,
                                name = "Rank",
                           row_names_gp = gpar(fontsize = 6),
                          column_names_gp = gpar(fontsize = 5),
                          show_column_names = FALSE,
                          show_row_names = FALSE,
                          Col = col_fun,
                           heatmap_legend_param = list(
                              legend_direction = "horizontal") ,
                          top_annotation = ha, 
                          right_annotation = hr)
  
   
  pdf(paste0("~/OneDrive - The University of Colorado Denver/Projects/subPhenoDS/results/partial_correlation_results/figures/",key_pathway, "_", "_mediation_genes_heatmap_", Sys.Date(),".pdf"), width = 8, height = 12)
  draw(p1, heatmap_legend_side = "bottom")
  dev.off()
  
  return(rankings_df)
}








