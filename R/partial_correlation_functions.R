#### This is a script for the partial correlation analysis and visualization helper functions.

## libraries for functions
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(fgsea)
library(RColorBrewer)

source("~/OneDrive - The University of Colorado Denver/Projects/ImmunoMetabolicSubtypes/R/enrichment_functions.R")



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

gseaFinder <- function(cyt_met_example, key_pathway, tmp_gsea_results_formatted, tmp_gsea_results) {

  cyt_met_labels <- names(tmp_gsea_results_formatted) 
  cyt_met_labels <- cyt_met_labels[grepl("NES", cyt_met_labels)]
  cyt_met_labels <- gsub(" NES", "", cyt_met_labels)
  path_id <- which(cyt_met_example ==cyt_met_labels)
  
  tmp_gsea <- lapply(tmp_gsea_results, function(x)
    x[[path_id]])
  
  res <- lapply(tmp_gsea, function(x) {
    x %>% filter(pathway == key_pathway)
  })
  
 final_res <- bind_rows(res) %>% dplyr::select(c(padj,NES,leadingEdge)) %>% as.data.frame()
  final_res <- final_res[1,]
return(final_res)
  
}

## heatmapDataPrep is a function for preparing the data for visualization in multiHeatmap
## Inputs: key_pathways    character vector of pathways to label
##         met_path        table of metabolite pathway annotations
##         pathways_all    list of character vecotors of genes belonging to hallmark pathways
##         formatted_gsea  table of gsea results with NES and padj values for each cytokine-metabolite relationships
##         threshold       the rank threshold for a gene
##         frequency       the frequency of gene ranked below the threshold
## Outputs: pdf file of figure


heatmapDataPrep <- function(key_pathways,pathways_all, met_path, formatted_gsea, threshold = 200, frequency = 5, z_scores ) {
    
    
  cyt_met_examples <- formatted_gsea %>%                                                  
    rownames_to_column("pathway") %>%
    filter(pathway %in% key_pathways) %>%
    pivot_longer(cols = contains("NES"), values_to = "NES", names_to = "cyt_met_NES") %>%
     pivot_longer(cols = contains("padj"), values_to = "padj", names_to = "cyt_met_padj") %>%
    filter(gsub(" NES", "", cyt_met_NES) == gsub(" padj", "", cyt_met_padj)) %>%
    mutate(cyt_met = gsub(" NES", "", cyt_met_NES)) %>%
    select(pathway, cyt_met, NES, padj) %>%
      filter(padj < .05 & NES > 0)    %>%
      arrange(padj)    %>%
      top_n(-100)    %>%
    distinct(cyt_met) %>%
    .$cyt_met
  
  
  res <- lapply(cyt_met_examples, function(x) {
    tmp <- z_scores%>%
      filter(cyt_met == x) %>%
      arrange(desc(combined_Z)) %>%
      mutate(rank = row_number()) %>%
      select(gene,rank)
    names(tmp)[2] <- x
    return(tmp)
  })

  
  rankings_df <- res %>% reduce(left_join, by = "gene") %>% column_to_rownames("gene")
  
  rankings_df <- as.matrix(rankings_df)
  
  tokeep <- apply(rankings_df, 1, function(x) sum(x < threshold) > frequency)
  toplot <- t(rankings_df[tokeep, ])

    return(toplot)
}


## geneSetUnion is a function for preparing the data for visualization in multiHeatmap, filtering the genes by geneset in the key_pathways
## Inputs: key_pathways    character vector of pathways to label
##         met_path        table of metabolite pathway annotations
##         pathways_all    list of character vecotors of genes belonging to hallmark pathways
##         formatted_gsea  table of gsea results with NES and padj values for each cytokine-metabolite relationships
##         threshold       the rank threshold for a gene
##         frequency       the frequency of gene ranked below the threshold
##         num_cyt_met     the number of cyt_mets to include in the plot
## Outputs: pdf file of figure


geneSetUnion <- function(key_pathways,pathways_all, met_path, formatted_gsea, threshold = 200, frequency = 5, z_scores, num_cyt_met = 100, padj_threshold = .05 ) {
    
    

   cyt_met_examples <- formatted_gsea %>%                                                  
        rownames_to_column("pathway") %>%
        filter(pathway %in% key_pathways) %>%
        pivot_longer(cols = contains("NES"), values_to = "NES", names_to = "cyt_met_NES") %>%
        pivot_longer(cols = contains("padj"), values_to = "padj", names_to = "cyt_met_padj") %>%
        filter(gsub(" NES", "", cyt_met_NES) == gsub(" padj", "", cyt_met_padj)) %>%
        mutate(cyt_met = gsub(" NES", "", cyt_met_NES)) %>%
        dplyr::select(pathway, cyt_met, NES, padj) %>%
        filter(padj < padj_threshold & NES > 0) %>%
        arrange(padj) %>%
        top_n(n = -num_cyt_met) %>%
        distinct(cyt_met) %>%
        .$cyt_met
  
    res <- lapply(cyt_met_examples, function(x) {
        tmp <- z_scores%>%
            filter(cyt_met == x) %>%
            arrange(desc(combined_Z)) %>%
            mutate(rank = row_number()) %>%
            dplyr::select(gene,rank)
        names(tmp)[2] <- x
        return(tmp)
    })
  
    rankings_df <- res %>% reduce(left_join, by = "gene") %>% column_to_rownames("gene")
  
    rankings_df <- as.matrix(rankings_df)
    
    tokeep <- unique(unlist(pathways_all[key_pathways]))
    toplot <- rankings_df %>%
        t() %>%
        as.data.frame() 
    tokeep <- tokeep[tokeep %in% names(toplot)]
    toplot <- toplot %>%
        dplyr::select(tokeep)

    agg_rank_gene_set <- lapply(key_pathways, function(pathway) {
        gene_set <- tokeep[tokeep %in% pathways_all[[pathway]]]
        mean_ranks <- toplot %>%
            dplyr::select(gene_set)%>%
            summarise(across(everything(), mean)) %>%
            pivot_longer(cols = everything()) %>%
            arrange(value) %>%
            .$name
        return(mean_ranks)
    })

    agg_rank_gene_set <- unlist(agg_rank_gene_set)
    agg_rank_gene_set <- agg_rank_gene_set[!duplicated(agg_rank_gene_set)]
    toplot <- toplot %>%
        dplyr::select(agg_rank_gene_set)
    
    return(toplot)
}



### multiHeatmap is a function for visualizing the gene ranks and gsea results over the cyt-met relationships
## Inputs: toplot1         dataframe of cyt-met X gene ranks
##         toplot2         cluster results of gene set enrichment for cyt-met relationships 
##         key_pathways    character vector of pathways to label
##         met_path        table of metabolite pathway annotations
##         pathways_all    list of character vecotors of genes belonging to hallmark pathways
##         formatted_gsea  table of gsea results with NES and padj values for each cytokine-metabolite relationships
##         full_results    table of full partial correlation results
##         output_file     file where to output the figures
## Outputs: pdf file of figure

multiHeatmap <- function(toplot1, toplot2, key_pathways,pathways_all, met_path, formatted_gsea, full_results, output_file) {
    

    names(toplot2) <- gsub(" NES", "", names(toplot2))

    toplot1 <- toplot1 %>%
        as.data.frame() %>%
        filter(rownames(toplot1) %in% names(toplot2))

    metabolites <-  full_results %>%
        filter(cyt_met %in% rownames(toplot1)) %>%
        distinct(cyt_met, .keep_all = TRUE) %>%
        left_join(met_path, by = c('metabolite' = 'name'))    %>%
        arrange(Pathway)

    toplot1 <- toplot1[match( metabolites$cyt_met, rownames(toplot1)),]


    toplot2 <- toplot2 %>%
        dplyr::filter(rownames(toplot2) %in% key_pathways) %>%
        dplyr::select(rownames(toplot1)) %>%
        t() %>%
        as.data.frame() %>%
        arrange(match(rownames(toplot1), rownames(toplot2)))
    
  names(toplot2) <- gsub("HALLMARK_","", names(toplot2))
  names(toplot2) <- gsub("_"," ", names(toplot2))
  
 
    ## adding heatmap annotation for gene pathways
    ## col_anno <- lapply(key_pathways, function(x) {
    ##     as.data.frame(as.factor( ifelse(colnames(toplot1) %in% pathways_all[[x]], 1,0)))
    ## })
    ## col_anno <- col_anno %>%
    ##     bind_cols(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) 

    ## rownames(col_anno) <- colnames(toplot1)
    ## colnames(col_anno) <- key_pathways
    ## col_anno <- as.data.frame(col_anno)
    ## column_col <- lapply(key_pathways, function(x) {
    ##     c("1" = "black", "0" = "white")
    ## })
    ## names(column_col) <- key_pathways
    
    ## ha <- HeatmapAnnotation(df = col_anno, col = column_col, show_legend = rep(FALSE,length(key_pathways)))
    
    hr <- rowAnnotation(Class = metabolites$Pathway)

    col_fun <- colorRamp2(c(min(toplot1), median(apply(toplot1,2, median)), max(toplot1) ), c("blue", "white", "gray"))

    p1 <- ComplexHeatmap::Heatmap(as.matrix(toplot1),
                                  name = "Rank",
                                  row_names_gp = gpar(fontsize = 6),
                                  column_names_gp = gpar(fontsize = 5),
                                  show_column_names = FALSE,
                                  show_row_names = FALSE,
                                  show_row_dend = FALSE,
                                  row_order = rownames(toplot1),
                                  col = col_fun,
                                  heatmap_legend_param = list(
                                      legend_direction = "horizontal") ,
                                 # top_annotation = ha, 
                                  left_annotation = hr,
                                  width = unit(6,"in"))


  
  col_fun2 <- colorRamp2(c(-2,0, 2 ), c("white","gray", "red"))
  
    p2 <- Heatmap(name = "GSEA NES",
      as.matrix(toplot2),
      col = col_fun2,
      row_order = rownames(toplot2),
      show_row_names = FALSE,
      show_row_dend = FALSE,
      width = unit(4, "in"),
      column_names_side = "top",
      show_column_dend = FALSE)

    heatmaps <- p1+p2
   
  pdf(output_file, width = 12, height = 12)
  draw(heatmaps, heatmap_legend_side = "bottom")
    dev.off()

    return(toplot1)
  
}


### multiHeatmapOrderedGSEA is a function for visualizing the gene ranks and gsea results over the cyt-met relationships and ordering it by gene set enrichments
## Inputs: toplot1         dataframe of cyt-met X gene ranks
##         toplot2         cluster results of gene set enrichment for cyt-met relationships 
##         key_pathways    character vector of pathways to label
##         met_path        table of metabolite pathway annotations
##         pathways_all    list of character vecotors of genes belonging to hallmark pathways
##         formatted_gsea  table of gsea results with NES and padj values for each cytokine-metabolite relationships
##         full_results    table of full partial correlation results
##         output_file     file where to output the figures
##         order_pathway   character string of pathway to order the results by key pathay

## Outputs: pdf file of figure

multiHeatmapOrderedGSEA <- function(toplot1, toplot2, key_pathways,pathways_all, met_path, formatted_gsea, full_results, output_file, order_pathway,threshold = 1000, frequency = 5, col_vector,labels_vector) {
    
    names(toplot2) <- gsub(" NES", "", names(toplot2))

    toplot1 <- toplot1 %>%
        as.data.frame() %>%
        filter(rownames(toplot1) %in% names(toplot2))
    
   
      metabolites <-  full_results %>%
        filter(cyt_met %in% rownames(toplot1)) %>%
        distinct(cyt_met, .keep_all = TRUE) %>%
        left_join(met_path, by = c('metabolite' = 'name'))    %>%
        arrange(Pathway)


    toplot2 <- toplot2 %>%
        dplyr::filter(rownames(toplot2) %in% key_pathways) %>%
        dplyr::select(rownames(toplot1)) %>%
        t() %>%
        as.data.frame() %>%
        arrange(-!!sym(order_pathway), .keep_all = TRUE)

    toplot2 <- toplot2[, key_pathways]

    toplot1 <- toplot1[match(rownames(toplot2), rownames(toplot1)),]

    metabolites <- metabolites[match(rownames(toplot1), metabolites$cyt_met),]
    
    names(toplot2) <- gsub("HALLMARK_","", names(toplot2))
    names(toplot2) <- gsub("_"," ", names(toplot2))

    
    hr <- rowAnnotation(Class = metabolites$Pathway,
                        col = list(Class = col_vector),
                        show_annotation_name = FALSE,
                        annotation_legend_param = list(
                            Class = list(
                                nrow = 9))
                        
                        )

    ##adding heatmap annotation for gene pathways
    col_anno <- lapply(key_pathways, function(x) {
        as.data.frame(as.factor( ifelse(colnames(toplot1) %in% pathways_all[[x]], 1,0)))
    })
    col_anno <- col_anno %>%
        bind_cols(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) 
    rownames(col_anno) <- colnames(toplot1)
    colnames(col_anno) <- key_pathways
    col_anno <- as.data.frame(col_anno)
    column_col <- lapply(key_pathways, function(x) {
        c("1" = "black", "0" = "white")
    })
    names(column_col) <- key_pathways
    
    ha <- HeatmapAnnotation(df = col_anno,
                            col = column_col,
                            show_legend = rep(FALSE,length(key_pathways)),
                            annotation_name_side = "left",
                            annotation_name_gp = gpar(fontsize = 6),
                            show_annotation_name = FALSE
                            )

    ## col_fun <- colorRamp2(c(min(toplot1), 2000, max(toplot1) ), c("blue", "white", "gray"))
    ## col_fun <- colorRamp2(c(1, median(apply(toplot1,2,median)), 12624 ), c("purple", "white", "yellow"))
   col_fun <- colorRamp2(c(1, 4000, 12624 ), c("darkblue", "white", "gray"))

    p1 <- ComplexHeatmap::Heatmap(as.matrix(toplot1),
                                  name = "Rank",
                                  row_names_gp = gpar(fontsize = 6),
                                  column_names_gp = gpar(fontsize = 6),
                                  show_column_names = FALSE,
                                  show_row_names = FALSE,
                                  row_labels = row_names,
                                  row_names_side = "left",
                                  show_row_dend = FALSE,
                                  show_column_dend = FALSE,
                                  column_dend_reorder = FALSE,
                                  column_order = names(toplot1),
                                  row_order = rownames(toplot1),
                                  col = col_fun,
                                  top_annotation = ha,
                                  heatmap_legend_param = list(
                                      legend_direction = "horizontal",) ,
                                  ##  left_annotation = hr,
                                  ## width = unit(6,"in"))
                                  )
  
  col_fun2 <- colorRamp2(c(-3,0, 3 ), c("white","gray", "red"))
    row_names <- ifelse(rownames(toplot1) %in% labels_vector, rownames(toplot1), "")  
    p2 <- Heatmap(name = "GSEA NES",
      as.matrix(toplot2),
      col = col_fun2,
      row_order = rownames(toplot2),
      show_row_names = TRUE,
      row_labels = row_names,
      row_names_side = "right",
      row_names_gp = gpar(fontsize = 6),
      show_row_dend = FALSE,
      column_order = names(toplot2),
      width = unit(.5, "in"),
      column_names_side = "top",
      show_column_dend = FALSE,
      column_names_gp = gpar(fontsize = 8),
      column_names_rot = 45,
      show_column_names = FALSE,
      heatmap_legend_param = list(legend_direction = "horizontal")
      
      )

    heatmaps <- p1+p2
    
  pdf(output_file, width = 4.5, height = 4.5)
  draw(heatmaps, heatmap_legend_side = "bottom")
  dev.off()

  
}


## OrderedGSEA is a function for returning the top ranking genes and cytokine-metabolite relationships given a key pathway.
## Inputs:
## toplot1      dataframe of gene mediation rankings across cyt-met relationships
## toplot2      GSEA NES scores across pathways and cyt-met relationships within a cluster
## 

OrderedGSEA <- function(formatted_gsea, key_pathway, p_thresh ) {

    cyt_met_examples <- formatted_gsea %>%                                                  
        rownames_to_column("pathway") %>%
        filter(pathway %in% key_pathway) %>%
        pivot_longer(cols = contains("NES"), values_to = "NES", names_to = "cyt_met_NES") %>%
        pivot_longer(cols = contains("padj"), values_to = "padj", names_to = "cyt_met_padj") %>%
        filter(gsub(" NES", "", cyt_met_NES) == gsub(" padj", "", cyt_met_padj)) %>%
        mutate(cyt_met = gsub(" NES", "", cyt_met_NES)) %>%
        dplyr::select(pathway, cyt_met, NES, padj) %>%
        filter(padj < p_thresh & NES > 0)    %>%
        arrange(padj)



    return(cyt_met_examples)
    
    ## names(toplot2) <- gsub(" NES", "", names(toplot2))

    ## toplot1 <- toplot1 %>%
    ##     as.data.frame() %>%
    ##     filter(rownames(toplot1) %in% names(toplot2))
    
   
    ## toplot2 <- toplot2 %>%
    ##     filter(rownames(toplot2) == order_pathway) %>%
    ##     dplyr::select(rownames(toplot1)) %>%
    ##     t() %>%
    ##     as.data.frame() %>%
    ##     arrange(-!!sym(order_pathway), .keep_all = TRUE)

  
    ## toplot1 <- toplot1[match(rownames(toplot2), rownames(toplot1)),]

       
    ## names(toplot2) <- gsub("HALLMARK_","", names(toplot2))
    ## names(toplot2) <- gsub("_"," ", names(toplot2))

    ## return(toplot2)
    

}







filterCyt_Mets <- function(toplot, formatted_gsea, key_pathway, gene_sets, num_cyt_met) {

    gene_set <- gene_sets[[key_pathway]]
    tmp_gsea <- formatted_gsea %>%
        dplyr::select(contains(rownames(toplot1))) %>%
        t() %>%
        as.data.frame() %>%
        arrange(-!!sym(key_pathway))

    cyt_mets <- rownames(tmp_gsea)[1:num_cyt_met]

    cyt_mets <- gsub(" NES", "", cyt_mets)
    
    return(cyt_mets)

}


filterCyt_Mets_unbiased <- function(sig_list, formatted_gsea, key_pathway, gene_sets, num_cyt_met, key_cluster, all_results) {

    cluster_sig <- all_results %>%
        filter(cluster == key_cluster) %>%
        filter(cor_cyt_met_p < .05) %>%
        distinct(cyt_met)

    sig <- intersect(sig_list, cluster_sig$cyt_met)
    tmp_gsea <- formatted_gsea %>%
        dplyr::select(contains(sig)) %>%
        t()    %>%
         as.data.frame()    %>%
        arrange(-!!sym(key_pathway))

    cyt_mets <- rownames(tmp_gsea)[grepl("NES", rownames(tmp_gsea))]

    cyt_mets <- gsub(" NES", "", cyt_mets)


    return(cyt_mets)

}


clusterGeneRanks <- function(cyt_mets, long_z_list, key_pathway, gene_sets ) {

    clusters <- names(long_z_list)

    gene_set <- gene_sets[[key_pathway]]
    
    results <- lapply(clusters, function(y) {

        res <- lapply(cyt_mets, function(x) {
            tmp <- long_z_list[[y]]            %>%
                filter(cyt_met == x)            %>%
                arrange(desc(combined_Z)) %>%
                mutate(rank = row_number()) %>%
                dplyr::select(gene,rank)
            names(tmp)[2] <- x
            return(tmp)
        })

  
        rankings_df <- res %>% reduce(left_join, by = "gene") %>% column_to_rownames("gene")
        
        rankings_df <- as.matrix(rankings_df)
                
        toplot <- t(rankings_df[rownames(rankings_df) %in% gene_set, ])

        tokeep <- gene_set[ gene_set %in% colnames(toplot)]

        toplot <- toplot %>%
            as.data.frame() %>%
            dplyr::select(tokeep)

        
        
        return(toplot)
    })

    mean_ranks <- results[[1]] %>%
            summarise(across(everything(), mean)) %>%
            pivot_longer(cols = everything()) %>%
            arrange(value) %>%
            .$name

    results <- lapply(results, function(x) {
        tmp <- x %>%
            dplyr::select(mean_ranks)
        return(tmp)
    })
    names(results) <- names(long_z_list)
    return(results)

}



clusterHeatmap <- function(cluster_rank, name, cluster_correlations) {



    correlation_data <- cluster_correlations %>%
        filter(cyt_met %in% rownames(cluster_rank)) %>%
        arrange(match(cyt_met, rownames(cluster_rank)))

    row_col_fun <- colorRamp2(c(-.6,0,.6), c("blue", "white", "red"))
    ha <- rowAnnotation(cor = correlation_data$r,
                        col = list(cor = row_col_fun),
                        show_legend = F,
                        show_annotation_name = F
                        )


        
    col_fun <- colorRamp2(c(1, 2000, 12624 ), c("darkblue", "white", "gray"))
    p1 <- Heatmap(as.matrix(cluster_rank),
                  name = name,
            row_order = rownames(cluster_rank),
            show_row_dend = FALSE,
            show_column_names = TRUE,
            column_names_gp = gpar(fontsize = 2),
            column_names_side = "top",
            row_names_side = "left",
            show_column_dend = FALSE,
            column_order = names(cluster_rank),
            col = col_fun,
            row_names_gp = gpar(fontsize =4),
            show_heatmap_legend = FALSE,
            row_title = name,
            left_annotation = ha
            )

    return(p1)

}


clusterHeatmapWrapper <- function(cluster_ranks, name,cluster_correlations, output_file) {

    clusters <- c("T21", "D21", "1"  , "2"  , "3"  , "4"  )
    heatmaps <- lapply(clusters , function(cluster) {
        clusterHeatmap(cluster_ranks[[cluster]], cluster, cluster_correlations[[cluster]])
    })
    all_heatmaps <- heatmaps[[1]]%v%heatmaps[[2]]%v%heatmaps[[3]]%v%heatmaps[[4]]%v%heatmaps[[5]]%v%heatmaps[[6]]
    
    pdf(output_file)
    draw(all_heatmaps, column_title = name)
    dev.off()
}


### Cluster Diff Plot








#### fullHeatmap is a function for creating a heatmap of the GSEA results of the partial correlation analysis over all the cyt-met relationships

fullHeatmap <- function(input, title) {
 col_fun = colorRamp2(c(min(input),0, max(input) ), c("white","gray", "red"))  
  var_order <- apply(input, 1, function(x) var(x))
  var_order <- var_order[order(var_order, decreasing = TRUE)]
  
  p1 <- Heatmap(
    as.matrix(input), 
    row_names_gp = gpar(fontsize = 8), 
    row_dend_reorder = FALSE, 
    row_order = names(var_order),
    row_title = "",
    column_names_gp = gpar(fontsize = 6), 
    column_names_rot = 60,
    show_column_names = FALSE,
    name = "NES", 
    column_title = title, 
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      legend_direction = "horizontal"), 
    col = col_fun
  )
  pdf(paste0("~/OneDrive - The University of Colorado Denver/Projects/subPhenoDS/results/partial_correlation_results/figures/",title,"_full_heatmap.pdf"))
  draw(p1,heatmap_legend_side="bottom")
  dev.off()
}



mediationSignatures <- function(key_pathway, tmp_gsea_formatted, tmp_toplot_clusters) {
  
  tmp_input <-  tmp_gsea_results_formatted %>%                                              
    filter(rownames(tmp_toplot_clusters) == key_pathway) %>%
    rownames_to_column("pathway") %>%
    pivot_longer(cols = contains("NES"), values_to = "NES", names_to = "cyt_met_NES") %>%
    pivot_longer(cols = contains("padj"), values_to = "padj", names_to = "cyt_met_padj") %>% 
    filter(gsub(" NES", "", cyt_met_NES) == gsub(" padj", "", cyt_met_padj)) %>%
    mutate(cyt_met = gsub(" NES", "", cyt_met_NES)) %>%
    dplyr::select(pathway, cyt_met, NES, padj) %>%
    filter(padj < .05 & NES > 0 & pathway == key_pathway) %>%
    arrange(padj)
  
  signature <- tmp_input$cyt_met
  
  return(signature)
  
} 


clusterPathwayEnrichment <- function(gsea_results) {
  
    output <-  lapply(gsea_results, function(x) {
        path_scores <- lapply(x, function(y) {
            y %>%
                mutate(one_sided_p = pval/2) %>%
                mutate(signed_p = ifelse(NES < 0, 1-one_sided_p, one_sided_p)) 

        })
        
        summary_scores <- path_scores %>% 
            bind_rows() %>%
            filter(!is.na(signed_p)) %>%
            group_by(pathway) %>%
            summarise(combined_chisq = -2 * sum(log(signed_p)),df = 2* n(), .groups = "keep") %>%
            mutate(combined_p = pchisq(combined_chisq,df,lower.tail = FALSE, log.p = FALSE))

        
        summary_scores$adj_combined_p <- p.adjust(summary_scores$combined_p, method = "fdr")
        
        return(summary_scores)        
    })

    return(output  )
} 


clusterPathwayRatios <- function(gsea_results) {
  
    output <-  lapply(gsea_results, function(x) {

        path_scores <- lapply(x, function(y) {
            y %>%
                mutate(sig_pos = ifelse(NES > 0 &  padj < .05, TRUE, FALSE)) %>%
                dplyr::select(pathway, padj,NES, sig_pos)

        })
        
        summary_scores <- path_scores %>% 
            bind_rows() %>%
            filter(!is.na(NES)) %>%
            group_by(pathway) %>%
            summarise(sig_cyt_met = sum(sig_pos), sig_cyt_met_perc = sum(sig_pos)/n(), .groups = "keep") 


        
                
        return(summary_scores)        
    })

    return(output  )
} 



clusterMetabolitePathwayEnrichment <- function(gsea_results) {

    output <-  lapply(gsea_results, function(x) {

        path_scores <- lapply(x, function(y) {
            y %>%
                mutate(sig_pos = ifelse(NES > 0 &  padj < .05, TRUE, FALSE)) %>%
                dplyr::select(pathway, padj,NES, sig_pos)

        })
        
        summary_scores <- path_scores %>% 
            bind_rows() %>%
            filter(!is.na(NES)) %>%
            group_by(pathway) %>%
            summarise(sig_cyt_met = sum(sig_pos), sig_cyt_met_perc = sum(sig_pos)/n(), .groups = "keep") 

        return(summary_scores)        
    })

    return(output  )
  
}

vector_match <- function(x, y) {
    grepl(x, y)
}



clusterCytMetNES <- function(key_pathway,cluster,gsea_results, num_cyt_mets = 20,met_path,met_class) {
    tmp_gsea <- gsea_results[c("T21","D21", cluster)]
    pathway_filtered_gsea <- lapply(tmp_gsea, function(x) {
        path_scores <- lapply(x, function(y) {
            y %>%
                filter(pathway == key_pathway ) 
        })
        path_scores <- path_scores %>%
            bind_rows() %>%
            mutate(cyt_met = names(tmp_gsea[[1]])) %>%
            dplyr::select(padj,NES,cyt_met)
    })
    all_pathway_gsea <- pathway_filtered_gsea %>% reduce(left_join, by = "cyt_met")

    res <- lapply(met_path[met_path$Pathway == met_class,"name"], function(x) {
        all_pathway_gsea %>% filter(grepl(x,cyt_met))
    })
    all_pathway_gsea <- bind_rows(res)
    names(all_pathway_gsea) <- c("T21.padj", "T21.NES", "cyt_met", "D21.padj", "D21.NES", "cluster.padj", "cluster.NES")
    all_pathway_gsea <- all_pathway_gsea %>%
        filter(!if_any(everything(), is.na)) %>%
        mutate(diff = cluster.NES - T21.NES) %>%
        filter( cluster.padj < .1) %>%
        arrange(-cluster.NES) %>%
        top_n(num_cyt_mets,diff) %>%
        dplyr::select(cyt_met, T21.NES,  cluster.NES) %>%
        distinct(cyt_met, .keep_all =TRUE)
        names(all_pathway_gsea) <- c( "cyt_met", "T21",   paste0("cluster", cluster))
    all_pathway_gsea$cyt_met <- factor(all_pathway_gsea$cyt_met, levels = all_pathway_gsea$cyt_met)
    all_pathway_gsea <- all_pathway_gsea %>%
        pivot_longer(-cyt_met, names_to = "cluster", values_to = "NES")
    return(all_pathway_gsea)
}
    
clusterGeneZscores <- function(key_gene, long_z_list, cyt_mets) {
    gene_z <- lapply(long_z, function(x) {
        x %>% filter(gene %in% key_gene & cyt_met %in% cyt_mets) %>%
            dplyr::select(cyt_met,combined_Z) %>%
            arrange(-combined_Z)
            
    })
    gene_z <- gene_z %>% reduce(left_join,by = "cyt_met") %>% column_to_rownames("cyt_met")
    names(gene_z) <- names(long_z_list)

    return(gene_z)
}



clusterDiffPlot <- function(pathway, key_cluster = 1, gsea_results. = gsea_results,num_cyt_met) {

    cluster_gsea <- clusterCytMetNES(key_pathway = pathway,
                                     cluster = key_cluster,
                                     gsea_results = gsea_results.,
                                     num_cyt_met = num_cyt_met)

    p1<- ggplot(data = cluster_gsea, aes(x = cyt_met, y = NES, color = cluster )) +
        geom_point() +
        geom_linerange(inherit.aes = F,
                       data = spread(cluster_gsea, cluster, NES),
                       aes_string(x = "cyt_met", ymin = "T21", ymax =paste0("cluster",key_cluster) ), linetype = "dashed", color = "grey")+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 80, hjust = 1),
              plot.title = element_text(hjust = .5),
              legend.title = element_blank(),
              legend.position = c(.9,.8)
              ) +
        ylab("NES") +
        xlab("") +
        ggtitle(paste0("Cluster ",key_cluster,"     ", pathway))
    

    pdf(paste0("~/OneDrive - The University of Colorado Denver/Projects/ImmunoMetabolicSubtypes/results/cluster_NES_diff/",key_cluster,"_",pathway,".pdf"),
        height= 6,
        width = 12)
    print(p1)
    dev.off()
}


topMediatorsCluster <- function(long_z_list, cyt_mets, cluster, key_pathway, pathways_all, all_correlations) {
    ## sig_cyt_mets <- all_correlations[[as.character(cluster)]] %>%
    ##     filter(fdr < .1) %>%
    ##     mutate(cyt_met = paste0(A,"-", B)) %>%
    ##     .$cyt_met
    ## cyt_mets <- cyt_mets[cyt_mets %in% sig_cyt_mets]
    genes <- sapply(cyt_mets, function(x) {
        long_z_list[[as.character(cluster)]] %>%
            filter(cyt_met == x) %>%
           filter(gene %in% pathways_all[[key_pathway]]) %>%
            top_n(1,combined_Z) %>%
            .$gene
    })
            gene_z <- lapply(long_z_list, function(x) {
        tmp <- sapply(1:length(cyt_mets), function(y) {
            x %>% filter(gene== genes[y] & cyt_met ==cyt_mets[y]) %>%
                dplyr::select(gene,cyt_met,combined_Z) 
        })
        tmp %>% t() %>% as.data.frame()
    })
    gene_z <- gene_z %>% reduce(left_join,by = c("gene","cyt_met")) %>% column_to_rownames("cyt_met")
    names(gene_z)[2:7] <- names(long_z_list)
    return(gene_z)
}
