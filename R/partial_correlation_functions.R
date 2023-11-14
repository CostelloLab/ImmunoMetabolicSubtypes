#### This is a script for the partial correlation analysis and visualization helper functions.

## libraries for functions
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(fgsea)





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

### multiHeatmap is a function for visualizing the gene ranks and gsea results over the cyt-met relationships
## Inputs: key_pathways    character vector of pathways to label
##         met_path        table of metabolite pathway annotations
##         pathways_all    list of character vecotors of genes belonging to hallmark pathways
##         formatted_gsea  table of gsea results with NES and padj values for each cytokine-metabolite relationships
##         output_file     file where to output the figures
## Outputs: pdf file of figure


multiHeatmap <- function(key_pathways,pathways_all, met_path, formatted_gsea, output_file) {
    
    met_path <- met_path %>%
        as.data.frame() %>%
        rownames_to_column("metabolite") %>%
        pivot_longer(cols = -metabolite, names_to = "pathway") %>%
        filter(value == 1)
    
  cyt_met_examples <- formatted_gsea %>%                                                  
    rownames_to_column("pathway") %>%
    filter(pathway %in% key_pathways) %>%
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
  
  tokeep <- apply(rankings_df, 1, function(x) sum(x < 100) > 20)
    toplot <- t(rankings_df[tokeep, ])

    ## adding heatmap annotation for gene pathways
    col_anno <- lapply(key_pathways, function(x) {
        as.factor( ifelse(colnames(toplot) %in% pathways_all[[x]], 1,0))
    })
    col_anno <- col_anno %>%
        bind_rows()
    names(col_anno) <- colnames(toplot)
    rownames(col_anno) <- key_pathways
    
 ## IFN_gamma <-as.factor( ifelse(colnames(toplot) %in% pathways_all$HALLMARK_INTERFERON_GAMMA_RESPONSE, 1,0))
 ## names(IFN_gamma) <- colnames(toplot)
 ##  Heme_Metabolism <-as.factor( ifelse(colnames(toplot) %in% pathways_all$HALLMARK_HEME_METABOLISM, 1,0))
 ## names(Heme_Metabolism) <- colnames(toplot)
 ##  Oxidative_Phosphorylation <-as.factor( ifelse(colnames(toplot) %in% pathways_all$HALLMARK_OXIDATIVE_PHOSPHORYLATION, 1,0))
 ## names(Oxidative_Phosphorylation) <- colnames(toplot)
 
 

  


    col_fun <- colorRamp2(c(min(toplot), median(toplot), max(toplot) ), c("blue", "white", "gray"))

    ha <- HeatmapAnnotation(col_anno)
    

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
    
   
  pdf(output_file, width = 8, height = 12)
  draw(p1, heatmap_legend_side = "bottom")
  dev.off()
  
}








