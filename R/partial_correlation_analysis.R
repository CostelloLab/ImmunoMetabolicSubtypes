


load("../DS-transfer/data/HTP_transcription_FPKMs_wide_missing_variance_filtered.Rdata")

correlations <- read.xlsx("./results/tables/cluster_correlations_D21.xlsx")

exp <- expression_list$expression

## log transform the FPKMs
exp[, 2:ncol(exp)] <- log2(exp[, 2:ncol(exp)] + 1e-6)

## scale the log transformed FPKMs
exp[, 2:ncol(exp)] <- apply(exp[, 2:ncol(exp)], 2, scale)

#### 'toplot' is from the correlation comparisons above.
exp <- exp %>%
    filter(LabID %in% toplot$Row.names) %>%
    arrange(match(LabID, toplot$Row.names))

toplot <- toplot %>%
    filter(Row.names %in% exp$LabID) %>%
    arrange(match( Row.names, exp$LabID))
## There are only 290 overlapping individuals

table(toplot$clustering)


## create dataset with all the omic measurements
omics_data <- exp %>%
    inner_join(toplot, by = c("LabID" = "Row.names"))

## Create a karyotype field
omics_data$karyotype <- ifelse(omics_data$clustering == "D21", "D21", "T21")



## calculate the partial correlation overall for significant cytokine metabolite pairs and all genes
library(ppcor)

## partialCor is a function for calculating and reporting the total and partial correlations in a dataset
partialCor <- function(data, cyt, met, gene, cluster) {
   
    ## filter out all T21 subjects or specific cluster
    if(cluster == "T21") {
        tmp <- data %>%
            filter(karyotype == cluster)        %>%
            dplyr::select(c(cyt, met, gene)) %>%
            as.data.frame()
    } else {
         tmp <- data %>%
            filter(clustering == cluster) %>%
            dplyr::select(c(cyt, met, gene)) %>%
             as.data.frame()
    }
    
    ## Testing the total correlation in the data
    total_correlation <- cor.test(tmp[, cyt], tmp[,met])
    total_correlation_p <- total_correlation$p.value
    total_correlation_r <- total_correlation$estimate
    partial_correlation <- pcor(tmp)
    partial_correlation_r <- partial_correlation$estimate[1,2]
    partial_correlation_p <- partial_correlation$p.value[1,2]

    return(c(total_correlation_r, total_correlation_p, partial_correlation_r, partial_correlation_p))
}


partial_cor(omics_data, cyt, met, gene, cluster)


partialCorWrapper <- function(tmp_correlations, data, clusters, genes) {
    ## list of results for T21, D21, and each cluster
    partial_cor_results <- list()

    partial_cor_results <- lapply(clusters, function(x) {
        print(x)
        
        tmp_correlations <- tmp_correlations %>%
            dplyr::select("source", "target",paste0(x, "_p"), paste0(x, "_r"))

        ## filter to onlys significant correlations within the cluster
        tmp_correlations <- tmp_correlations[!is.na(tmp_correlations[, 3]), ]

        if(nrow(tmp_correlations) >=1) { 

        
            tmp_partial_cor_results <-  lapply(1:nrow(tmp_correlations), function(y) {

                print(y/nrow(tmp_correlations))

                cyt <- tmp_correlations[y, "source"]
                met <- tmp_correlations[y, "target"]
                
                tmp_results <- lapply(genes, function(z) {

                    return( c(cyt, met, z, partialCor(data, cyt, met, z, x)))
                    
                })
                tmp_results <- as.data.frame(do.call(rbind, tmp_results))
                return(tmp_results)
            })
            tmp_partial_cor_results <- as.data.frame(do.call(rbind, tmp_partial_cor_results))
            names(tmp_partial_cor_results) <- c("cytokine", "metabolite", "gene", "total_correlation_r", "total_correlation_p", "partial_correlation_r", "partial_correlation_p")
        } else {
            tmp_partial_cor_results <- NA
        }
        return(tmp_partial_cor_results)
    })
    
    names(partial_cor_results) <- clusters
    return(partial_cor_results)
}




cyt <- "IFN-gamma"
met <- "kynurenine"
gene <- "IDO1"
data <- omics_data
cluster <- "D21"




## setting the D21 p values to NA if not significant
## changing the names of the correlations 'All' field to 'T21'
names(correlations)[3:4] <- c("T21_p", "T21_r")
## removing the 'Cluster' from the names of correlations as well
names(correlations)[grepl("Cluster", names(correlations))] <- gsub("Cluster", "", names(correlations)[grepl("Cluster", names(correlations))])
correlations$D21_p <- ifelse(correlations$D21_p < .05, correlations$D21_p, NA)


tmp_correlations <- correlations
clusters <- c("T21", "D21", 1, 2, 3, 4, 5)
genes <- names(exp)[2:ncol(exp)]
partial_cor_results <- partialCorWrapper(tmp_correlations = tmp_correlations, data = omics_data, genes = genes, clusters = clusters)
save(partial_cor_results, file = "./results/partial_correlation_results.RData")

save(omics_data, file = "../data/processed/cyt_met_RNA_processed_2023-05-26.RData")

load('partial_correlation_results_T21.RData')
res <- partial_cor_results[[1]]
res[,4:7] <- apply(res[,4:7], 2, as.numeric)
res$diff <- res$partial_correlation_p-res$total_correlation_p

res_IFN <- res %>% filter(cytokine == "IFN-gamma") %>% arrange(-diff) %>% dplyr::select(gene,diff)
res_TNF <- res %>% filter(cytokine == "TNF-alpha") %>% arrange(-diff) %>% dplyr::select(gene,diff)

write.table(res_IFN, "res_IFN.rnk", row.names= FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(res_TNF, "res_TNF.rnk", row.names= FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



