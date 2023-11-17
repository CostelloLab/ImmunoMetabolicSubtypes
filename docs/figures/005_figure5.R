#### Figure 5
## source functions
source("~/OneDrive - The University of Colorado Denver/Projects/ImmunoMetabolicSubtypes/R/diff_to_Z_functions.R")

### Load the data
### Differential expression
load("~/OneDrive - The University of Colorado Denver/Projects/subPhenoDS/results/differential_expression/differential_expression_results_2023-09-08.RData")
## The first data frame compare D21 to T21, need to get the inverse of logFC
diff_genes[[1]]$logFC <- diff_genes[[1]]$logFC*-1
### partial correlation z scores
load("~/OneDrive - The University of Colorado Denver/Projects/subPhenoDS/results/partial_correlation_results/full_z_percent_8_31.RData")
### Renmame the list list of z scores by karyotype/cluster
names(full_z_results) <- c("T21", "D21", 1,2,3,4,5)
### Hallmark pathways
pathways_all <- gmtPathways("~/OneDrive - The University of Colorado Denver/Projects/subPhenoDS/data/pathway/h.all.v2023.1.Hs.symbols.gmt")

##### T21 vs D21
diff_T21_D21 <- diffZ(full_z_results, cluster1 = "T21", cluster2 = "D21")
z_and_DEG_T21_D21 <- ZandDEG(DEG = diff_genes[[1]], delta_z = diff_T21_D21)


### Overall
corr_Z_DEG <- corrZandDEGwrapper(colnames(full_z_results[[1]]), diff_and_z = z_and_DEG_T21_D21)


ggplot(corr_Z_DEG, aes(x = cor, y = -log10(p.value)))+
    geom_point()

plotZandDEG(z_and_DEG = z_and_DEG_T21_D21,
            key_cyt_met = "MIP-3alpha-trans-4-Hydroxy-L-proline")

plotZandDEG(z_and_DEG = z_and_DEG_T21_D21,
            key_cyt_met = "IFN-gamma-kynurenine")

plotZandDEG(z_and_DEG = z_and_DEG_T21_D21,
            key_cyt_met = "SAA-Arginine")


## NB: What about by gene sets?

gene_set_corr_z_DEG <- GeneSetcorrZandDEGwrapper(colnames(full_z_results[[1]]), diff_and_z = z_and_DEG_T21_D21, pathways_all = pathways_all)

lapply(gene_set_corr_z_DEG, head)
lapply(gene_set_corr_z_DEG, tail)

GeneSetPlotZandDEG(z_and_DEG_T21_D21,
                   key_cyt_met = "MIP-3alpha-Lysine",
                   pathways_all = pathways_all,
                   pathway = "HALLMARK_INTERFERON_GAMMA_RESPONSE")

GeneSetPlotZandDEG(z_and_DEG_T21_D21,
                   key_cyt_met = "IFN-gamma-kynurenine",
                   pathways_all = pathways_all,
                   pathway = "HALLMARK_INTERFERON_GAMMA_RESPONSE")

GeneSetPlotZandDEG(z_and_DEG_T21_D21,
                   key_cyt_met = "MCP-1-Aspartate",
                   pathways_all = pathways_all,
                   pathway = "HALLMARK_HEME_METABOLISM")

GeneSetPlotZandDEG(z_and_DEG_T21_D21,
                   key_cyt_met = "IL-1RA-Glutamate",
                   pathways_all = pathways_all,
                   pathway = "HALLMARK_HEME_METABOLISM")


GeneSetPlotZandDEG(z_and_DEG_T21_D21,
                   key_cyt_met = "TNF-beta-5-Hydroxyindoleacetate",
                   pathways_all = pathways_all,
                   pathway = "HALLMARK_PANCREAS_BETA_CELLS")


GeneSetPlotZandDEG(z_and_DEG_T21_D21,
                   key_cyt_met = "IL-12/IL-23p40-propionyl-carnitine",
                   pathways_all = pathways_all,
                   pathway = "HALLMARK_OXIDATIVE_PHOSPHORYLATION")




source("~/OneDrive - The University of Colorado Denver/Projects/ImmunoMetabolicSubtypes/R/diff_to_Z_functions.R")
### By cluster
clusters <- c('1','2','3','4','5')
diff_Z_clusters <- lapply(clusters, function(cluster) diffZCluster(full_z_results, cluster1 = cluster))
names(diff_Z_clusters) <- clusters
Z_and_DEG_clusters <- lapply(clusters, function(cluster) ZandDEG(DEG = diff_genes[[cluster]], delta_z = diff_Z_clusters[[cluster]]))
names(Z_and_DEG_clusters) <- clusters

gene_set_corr_z_DEG_clusters <- lapply(clusters, function(cluster) {
    print(cluster)
    GeneSetcorrZandDEGwrapper(colnames(full_z_results[[1]]), diff_and_z = Z_and_DEG_clusters[[cluster]], pathways_all = pathways_all)
})




