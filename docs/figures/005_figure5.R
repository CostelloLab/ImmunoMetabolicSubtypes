#### Figure 5

### Necessary libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(fgsea)

### Load the data
### Differential expression
load("~/OneDrive - The University of Colorado Denver/Projects/subPhenoDS/results/differential_expression/differential_expression_results_2023-09-08.RData")
### partial correlation z scores
load("~/OneDrive - The University of Colorado Denver/Projects/subPhenoDS/results/partial_correlation_results/full_z_percent_8_31.RData")
### Hallmark pathways
pathways_all <- gmtPathways("~/OneDrive - The University of Colorado Denver/Projects/subPhenoDS/data/pathway/h.all.v2023.1.Hs.symbols.gmt")



### Renmame the list list of z scores by karyotype/cluster
names(full_z_results) <- c("T21", "D21", 1,2,3,4,5)


diff_T21_D21 <- diffZ(full_z_results, cluster1 = "T21", cluster2 = "D21")
z_and_DEG_T21_D21 <- ZandDEG(DEG = diff_genes[[1]], delta_z = diff_T21_D21)

plotZandDEG(z_and_DEG = z_and_DEG_T21_D21, key_cyt_met ="MIP-3alpha-trans-4-Hydroxy-L-proline")

corr_Z_DEG <- corrZandDEGwrapper(colnames(full_z_results[[1]]), diff_and_z = z_and_DEG_T21_D21)

ggplot(corr_Z_DEG, aes(x = cor, y = -log10(p.value)))+
    geom_point()

## NB: What about by gene sets?

gene_set_corr_z_DEG <- GeneSetcorrZandDEGwrapper(colnames(full_z_results[[1]]), diff_and_z = z_and_DEG_T21_D21, pathways_all = pathways_all)
GeneSetPlotZandDEG(z_and_DEG_T21_D21,
                   key_cyt_met = "IL-27-Histidine",
                   gene_set = pathways_all$HALLMARK_HEME_METABOLISM)




### Figure 5X


