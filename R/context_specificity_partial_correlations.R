

########## context_specificity_partial_correlation
library(tidyverse)
library(igraph)
library(parallel)
library(limma)
library(metap)
library(fgsea)
library(fgsea)
library(purrr)
library(psych)
library(Hmisc)
library(ppcor)
library(openxlsx)

### create the gene X cyt-met matrix
createCytMetGeneMatrix <- function(input_results, values = "percent_diff") {
  gene_cytmet_matrix <- input_results %>%
    mutate(cyt_met = paste0(cytokine,"-",metabolite)) %>%
    group_by(gene) %>%
    mutate(row = row_number()) %>%
    pivot_wider( id_cols = c(row, cyt_met), names_from = gene, values_from = all_of(values)) %>%
    dplyr::select(-row) %>%
    distinct(cyt_met, .keep_all = TRUE) %>%
    column_to_rownames("cyt_met") 
  
  gene_cytmet_matrix[sapply(gene_cytmet_matrix, is.infinite)] <- NA
  
  gene_cytmet_matrix <- gene_cytmet_matrix %>%
    replace(is.na(.),0)
  
  return(gene_cytmet_matrix)
}



#### z score by each axis
### by genes
zScoreMatrix <- function(input_matrix) {
  gene_z <- apply(input_matrix,2, scale  )
  cytmet_z <- apply(input_matrix,1,scale)
  rownames(gene_z) <- colnames(cytmet_z)
  rownames(cytmet_z) <- colnames(gene_z)
  return(list(gene_z = gene_z, cytmet_z = cytmet_z))
}

### combine z scores with stouffer's method
stouffers <-function(z1,z2) {
  res <- (z1+z2)/sqrt(2)
  return(res)
}

### combine p-values with Fisher's method


fishers <- function(p1,p2) {
  if(any(is.na(c(p1,p2)) )) {
    return(NA)
  } else {
    return(sumlog(c(p1,p2))$p)
  }
}

# fishers <- function(p1,p2) {
#   chi <- -2*log(p1) + -2*log(p2)
#   res <- pchisq(q = chi, df = 2)
#   return(res)
# }

## create the z score matrix
combineZMatrices <- function(input_list, combination_method = "stouffers") {
  if(combination_method == "stouffers"){
    full_z <- sapply(rownames(input_list$gene_z), function(i) {
      sapply(colnames(input_list$gene_z), function(j) {
        z1 <- input_list$gene_z[i,j]
        z2 <- input_list$cytmet_z[j,i]
        res <- stouffers(z1,z2) 
        return(res)
      })
    })
  } else {
    full_z <- sapply(rownames(input_list$gene_z), function(i) {
      sapply(colnames(input_list$gene_z), function(j) {
        z1 <- input_list$gene_z[i,j]
        z2 <- input_list$cytmet_z[j,i]
        res <- fishers(z1,z2) 
        return(res)
      })
    })
    
  }
  return(full_z)
}

## wrapper function to calculate the full_z matrices over each cluster
fullZWrapper <- function(cluster_results, clusters = c("T21", "D21", 1,2,3,4), values ="percent_diff", combination_method = "stouffers") {
  results_list <- mclapply(clusters, function(x) {
    tmp_results <- cluster_results %>% 
      filter(cluster == x)
    cyt_met_gene_matrix <- createCytMetGeneMatrix(tmp_results, values  = values)
    z_score_matrix <- zScoreMatrix(cyt_met_gene_matrix)
    full_z <- combineZMatrices(z_score_matrix, combination_method = combination_method)
    return(full_z)
  }, mc.cores = 6)
  return(results_list)
}




# cluster_full_z_results <- fullZWrapper(cluster_results)
# save(cluster_full_z_results, file = "cluster_full_z_percent_r.RData")
# 
# T21_full_z_results <- fullZWrapper(T21_results)
# save(full_z, file = "T21_full_z_percent_r.RData")
# 
# T21_results$abs_diff_p <- abs(T21_results$diff_p)
# T21_results$percent_diff_p <- 1- T21_results$abs_diff_p/T21_results$cor_cyt_met_p
# T21_cyt_gene_matrix <- createCytMetGeneMatrix(T21_results, values = "percent_diff_r")

#### combining p values of z scores doesn't make sense since z scores are negative
# T21_full_z_p_results <- fullZWrapper(T21_results, values = "partial_cor_cyt_met_p", combination_method = "fishers")
# save(T21_full_z_p_results, file = "T21_full_z_p.RData")
# 
# clusters_full_z_p_results <- fullZWrapper(cluster_results, values = "partial_cor_cyt_met_p", combination_method = "fishers")
# save(clusters_full_z_p_results, file = "clusters_full_z_p.RData")

### Filter the z score matrix based on significance
filterZ <- function(full_z){
  bonf_correction <- .1/(ncol(full_z))
  critical_value <- qnorm(bonf_correction/2)
  sig_z <- ifelse(abs(full_z) > abs(critical_value), full_z, NA)
  return(sig_z)
}



### create a network from the filtered matrix
createEdgeList <- function(sig_z) {
  edge_list <- sig_z %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -c("gene") ,  names_to = "cytmet" ) %>%
    filter(!is.na(value))
  return(edge_list)
} 


createGraph <- function(edge_list) {
  graph_z <- graph_from_edgelist(as.matrix(edge_list[,1:2]), directed = F)
  E(graph_z)$weight <- edge_list$value
  return(graph_z)
}



### calculating the weighted centrality

weightedCentrality <- function(node,edge_list) {
  n <- length(c(unique(edge_list$gene, unique(edge_list$cytmet))))
  node_subgraph <- edge_list %>% filter(gene == node | cytmet == node) 
  node_degree <- nrow(node_subgraph)
  summed_node_weight <- sum(abs(node_subgraph$value))
  weighted_centrality <- node_degree * summed_node_weight / n
  return(weighted_centrality)
}


weightedCentralitiesWrapper <- function(edge_lists) {
  weighted_centralities_list <- mclapply(edge_lists, function(y) {
    
    weighted_centralities <- sapply(unique(y$gene), function(x) weightedCentrality( x,y))
    weighted_centralities <- weighted_centralities[order(weighted_centralities)]
    return(weighted_centralities)
  }, mc.cores = 6)
  # return(weighted_centralities_list)
}


withinClusterDEG <- function(omics_data, cluster, gene_list) {
  
  design = model.matrix(~0+as.factor(omics_data$clustering), data = omics_data)
  colnames(design) = paste0("c",sort(unique(omics_data$clustering)))
  c0 = sprintf("contrast =  makeContrasts(contrast= c%s, levels = design)",cluster)
  eval(parse(text = c0))	
  tmp <- t(omics_data[, gene_list])
  fit = eBayes(contrasts.fit(lmFit(tmp,design),contrast))
  all.siglist = topTable(fit, adjust = 'f', number = nrow(tmp)+1,
                         confint = T)
  sig.genes <- all.siglist
  return(sig.genes)
}

### wide to long for each matrix in the list
wideToLong <- function(full_z_results) {
  res <- lapply(full_z_results, function(x) {
    tmp <-  x %>% 
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(values_to = "combined_Z", names_to = "cyt_met", cols = where(is.numeric) )
  })
  return(res)
}



### Combine the list of results into a single matrix
addClusterID <- function(full_z_results) {
  clusters = c("T21", "D21", 1,2,3,4,5)
  for(i in 1:length(clusters)){
    full_z_results[[i]] = full_z_results[[i]] %>%
      as.data.frame() %>%
      mutate(cluster = clusters[i])
  }
  return(full_z_results)
}

#### fgsea function 


### gsea wrapper
# pathways1 <- gmtPathways("c2.cp.kegg.v2023.1.Hs.symbols.gmt")
# pathways2 <- gmtPathways("c3.tft.v2023.1.Hs.symbols.gmt")
# pathways_all <- c(pathways1,pathways2)
pathways_all <- gmtPathways("./pathways/h.all.v2023.1.Hs.symbols.gmt")

gseaWrapper <- function(long_results_matrix){
  
  cyt_met_rel <- unique(long_results_matrix$`cyt_met`)
  long_results_matrix <- long_results_matrix %>%
    arrange(-combined_Z)
  
  print("ranking genes")
  gene_ranks <- mclapply(cyt_met_rel, function(x) {
    ranks <- long_results_matrix %>%
      filter(`cyt_met` == x) 
    gene_names <- ranks$gene
    ranks <- ranks$combined_Z
    names(ranks) <- gene_names
    ranks <- ranks[!is.na(ranks)]
    return(ranks)
    
  }, mc.cores =18)
  
  failed <-  which(sapply(gene_ranks, function(x) class(x)) == "try-error")
  if(length(failed) > 0) {print(paste("failed with ", cyt_met_rel[failed]))}
  
  print("gsea")
  gsea_results <- list()
  n <- length(cyt_met_rel)
  gsea_results <- mclapply(1:n, function(i) {
    fgsea(pathways = pathways, stats= gene_ranks[[i]], minSize = 15, maxSize = 500, 
          nPermSimple = 1000, nproc = 1)
  }, mc.cores = 18)
  print("done")
  names(gsea_results) <- cyt_met_rel
  classes <- lapply(gsea_results, function(x) class(x))
  return(gsea_results)
}



# pathways1 <- gmtPathways("c2.cp.kegg.v2023.1.Hs.symbols.gmt")
# pathways2 <- gmtPathways("c3.tft.v2023.1.Hs.symbols.gmt")
# pathways<- c(pathways1,pathways2)
pathways <- gmtPathways("./pathways/h.all.v2023.1.Hs.symbols.gmt")
gseaWrapperGeneExpression <- function(diff_expression_results ){
  
  diff_expression_results <- diff_expression_results %>%
    arrange(-logFC) 
  gene_ranks <- diff_expression_results$logFC
  names(gene_ranks) <- rownames(diff_expression_results)
  
  gsea_results <- fgsea(pathways = pathways, stats= gene_ranks, minSize = 15, maxSize = 500, 
                        nPermSimple = 100000, nproc = 1,eps = 0 )
  
  return(gsea_results)
}



```