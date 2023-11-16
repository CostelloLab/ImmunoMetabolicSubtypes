#### Figure 5

### Load the data
### Differential expression
load("~/OneDrive - The University of Colorado Denver/Projects/subPhenoDS/results/differential_expression/differential_expression_results_2023-09-08.RData")

## Fold Change correlation to z-score
key_cyt_met <- "IFN-gamma-kynurenine"

tmp_diff <- diff_genes[[1]] %>%
    rownames_to_column("gene")

diff_and_z <- long_z[[1]] %>%
    filter(cyt_met == key_cyt_met) %>%
    left_join(tmp_diff, by = "gene") 



ggplot(diff_and_z, aes(x = combined_Z, y = logFC)) +
    geom_point()



