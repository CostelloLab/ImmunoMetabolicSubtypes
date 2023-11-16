### Necessary libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(fgsea)


### diffZ is a function that finds the differential mediation scores between karyotypes/clusters

diffZ <- function(z_scores_list, cluster1, cluster2) {
    diff_z <- z_scores_list[[cluster1]] - z_scores_list[[cluster2]]
    return(diff_z)
}


### NB: For comparison, need a function that finds difference in Z score for 1 cluster vs all other clusters

## ZandDEG is a function for combining the change in Z scores between karyotypes/clusters and DEG by karyotypes/clusters for a specific cytokine_metabolite relationship
ZandDEG <- function(delta_z, DEG) {

    tmp_diff <- DEG %>%
        rownames_to_column("gene")

    diff_and_z <- delta_z  %>%
        as.data.frame() %>%
        rownames_to_column("gene")%>%
        left_join(tmp_diff, by = "gene")

    return(diff_and_z)
}

### plotZandDEG is a function for plotting the relationship between the change in Z score and DEG by karyotype/cluster
plotZandDEG <- function(z_and_DEG, key_cyt_met) {

    tmp_z_and_DEG <- z_and_DEG %>%
        select(logFC, key_cyt_met)

    names(tmp_z_and_DEG)[2] <- "deltaZ"
    
   p <- ggplot(tmp_z_and_DEG, aes(x = deltaZ,y = logFC)) +
       geom_point() +
       ggtitle(key_cyt_met) +
       theme_classic()+
       stat_cor() +
       geom_smooth(method = "lm")
    return(p)
}

### corrZandDEG is a function for calculating the correlation coefficient between the deltaZ and logFC
corrZandDEG <- function(z_and_DEG) {
    corr <- cor.test(x = z_and_DEG[,1], z_and_DEG[,2])
    return(corr)
}


### corrZandDEGwrapper is a wrapper function that outputs the correlation coefficients across cytokine-metabolite relationships
corrZandDEGwrapper <- function(cyt_met_examples, diff_and_z) {
    res <- sapply(cyt_met_examples, function(x) {
        tmp_Z_DEG <- diff_and_z %>%
            select(logFC, x)
        names(tmp_Z_DEG)[2] <- "deltaZ"
        tmp_corr <- corrZandDEG(tmp_Z_DEG)
        return(c(tmp_corr$estimate, tmp_corr$p.value))
    })
    res <- as.data.frame(t(res))
    names(res)[2] <- "p.value"
    res <- res %>%
        arrange(-cor)
    
    return(res)
}


### GeneSetCorrZandDEG is a function for calculating the correlation coefficient between the deltaZ and logFC for a specific Gene Set

GeneSetCorrZandDEG <- function(z_and_DEG, gene_set) {
    tmp_z_and_DEG <- z_and_DEG %>%
        filter(gene %in% gene_set) %>%
        select(-gene)
    corr <- cor.test(x = tmp_z_and_DEG[,1], tmp_z_and_DEG[,2])
    return(corr)
}


### GeneSetCorrZandDEGwrapper is a wrapper function that outputs the correlation coefficients across cytokine-metabolite relationships across gene sets
GeneSetcorrZandDEGwrapper <- function(cyt_met_examples, diff_and_z, pathways_all) {
    pathways <- names(pathways_all)
    pathway_res <- lapply(pathways, function(pathway) {
        res <- sapply(cyt_met_examples, function(x) {
            tmp_Z_DEG <- diff_and_z %>%
                select(gene,logFC, x)
            tmp_corr <- GeneSetCorrZandDEG(tmp_Z_DEG, pathways_all[[pathway]])
            return(c(tmp_corr$estimate, tmp_corr$p.value))
        })
        res <- as.data.frame(t(res))
        names(res)[2] <- "p.value"
        res <- res %>%
            arrange(-cor)
    })

    names(pathway_res) <- pathways
    return(pathway_res)
}

### GeneSetplotZandDEG is a function for plotting the relationship between the change in Z score and DEG by karyotype/cluster for a specific gene set
GeneSetPlotZandDEG <- function(z_and_DEG, key_cyt_met, pathways_all, pathway) {

    gene_set <- pathways_all[[pathway]]
    
    tmp_z_and_DEG <- z_and_DEG %>%
        filter(gene %in% gene_set) %>%
        select(logFC, key_cyt_met)
    
    names(tmp_z_and_DEG)[2] <- "deltaZ"
    
   p <- ggplot(tmp_z_and_DEG, aes(x = deltaZ,y = logFC)) +
       geom_point() +
       ggtitle(paste(key_cyt_met,pathway)) +
       theme_classic()+
       stat_cor() +
       geom_smooth(method = "lm")
    return(p)
}
