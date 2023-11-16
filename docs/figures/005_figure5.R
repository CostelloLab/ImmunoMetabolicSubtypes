#### Figure 5

### Load the data
### Differential expression
load("~/OneDrive - The University of Colorado Denver/Projects/subPhenoDS/results/differential_expression/differential_expression_results_2023-09-08.RData")
### partial correlation z scores
load("~/OneDrive - The University of Colorado Denver/Projects/subPhenoDS/results/partial_correlation_results/full_z_percent_8_31.RData")


### Necessary libraries
library(tidyverse)

### Renmame the list list of z scores by karyotype/cluster
names(full_z_results) <- c("T21", "D21", 1,2,3,4,5)

### diffZ is a function that finds the differential mediation scores between karyotypes/clusters

diffZ <- function(z_scores_list, cluster1, cluster2) {
    diff_z <- z_scores_list[[cluster1]] - z_scores_list[[cluster2]]
    return(diff_z)
}


### NB: For comparison, need a function that finds difference in Z score for 1 cluster vs all other clusters

## ZandDEG is a function for combining the change in Z scores between karyotypes/clusters and DEG by karyotypes/clusters for a specific cytokine_metabolite relationship

ZandDEG <- function(key_cyt_met, delta_z, DEG) {

    tmp_diff <- DEG %>%
        rownames_to_column("gene")

    diff_and_z <- delta_z  %>%
        as.data.frame() %>%
        select(key_cyt_met) %>%
        rownames_to_column("gene")%>%
        left_join(tmp_diff, by = "gene") 

    return(diff_and_z)
}


plotZandDEG <- function(z_and_DEG) {
    

### Figure 5X
diff_T21_D21 <- diffZ(full_z_results, cluster1 = "T21", cluster2 = "D21")
Z_DEG_T21_D21 <- ZandDEG(key_cyt_met = "IFN-gamma-kynurenine",DEG =diff_genes[[1]],delta_z =  diff_T21_D21)


ggplot(diff_and_z, aes(x = combined_Z, y = logFC)) +
    geom_point()



